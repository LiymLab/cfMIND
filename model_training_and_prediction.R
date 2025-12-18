rm(list = ls())
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(pROC)
  library(parallel)
  library(caret)
  library(Boruta)
  library(xgboost)
})

# ---------------- Arguments ----------------
args = commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop("Usage: Rscript model_training_and_prediction.R <feature_matrix_file> <prefix> <out_dir> <n_repeat> <n_fold> <n_thread>\n", call. = FALSE)
}

feature_matrix_file = args[1]
prefix   = args[2]
out_dir  = args[3]
n_repeat = as.numeric(args[4])
n_fold   = as.numeric(args[5])
n_thread = as.numeric(args[6])

# Ensure output directory exists
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
if (!grepl("/$", out_dir)) out_dir <- paste0(out_dir, "/")



# ---------------------- Load Data and Labels ----------------------
# Load feature matrix from tab file
if (!file.exists(feature_matrix_file)) {
  stop("Feature matrix file not found: ", feature_matrix_file)
}

# Read the feature matrix file
data_2d <- read.table(feature_matrix_file, header = TRUE, sep = "\t", 
                      row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
# 去除 tumor 列
if ("tumor" %in% colnames(data_2d)) {
  data_2d <- data_2d[, colnames(data_2d) != "tumor", drop = FALSE]
}
# Check if label column exists
if (!"label" %in% colnames(data_2d)) {
  stop("Label column not found in feature matrix file: ", feature_matrix_file)
}

# Extract labels and remove label column from feature matrix
sample_labels <- data_2d$label
data_2d <- data_2d[, !colnames(data_2d) %in% "label", drop = FALSE]

# Create labels_df for compatibility
labels_df <- data.frame(
  sample = rownames(data_2d),
  group = sample_labels,
  stringsAsFactors = FALSE
)

# Remove samples with missing labels
valid_samples <- !is.na(labels_df$group)
if (sum(valid_samples) == 0) {
  stop("No valid labels found in feature matrix file")
}

data_2d <- data_2d[valid_samples, , drop = FALSE]
labels_df <- labels_df[valid_samples, , drop = FALSE]
labels_df$sample <- trimws(labels_df$sample)
labels_df$group <- trimws(labels_df$group)

# Create sample to group mapping
sample_to_group <- setNames(trimws(labels_df$group), trimws(labels_df$sample))

get_true <- function(samples) {
  s <- trimws(samples)
  g <- unname(sample_to_group[s])
  if (anyNA(g)) warning("Some samples did not match labels in label file: ",
                        paste(unique(s[is.na(g)]), collapse = ", "))
  g
}

# Check that all samples in data_2d have corresponding labels
data_samples <- rownames(data_2d)
missing_labels <- setdiff(data_samples, labels_df$sample)
if (length(missing_labels) > 0) {
  stop("Missing labels for samples: ", paste(missing_labels, collapse = ", "))
}

# Filter data_2d to only include samples with labels
data_2d <- data_2d[data_samples %in% labels_df$sample, , drop = FALSE]

cat(sprintf("Loaded feature matrix: %d samples x %d features\n", nrow(data_2d), ncol(data_2d)))
cat(sprintf("Loaded labels for %d samples\n", nrow(labels_df)))

# If n_fold is 0, set it to the number of samples (Leave-One-Out CV)
if (n_fold == 0) {
  n_fold <- nrow(data_2d)
  cat(sprintf("Using Leave-One-Out CV with %d folds\n", n_fold))
}

# ---------------- Modeling & Prediction (Parallel) ----------------
# Get group labels for samples in the same order as data_2d
group_labels <- get_true(rownames(data_2d))
model_name   <- "XGBoost"
result_file  <- paste0(out_dir, prefix, "_", model_name, "_predict_probability.txt")
if (file.exists(result_file)) file.remove(result_file)

group_levels <- unique(group_labels)

# Generate folds for each repeat
folds_by_repeat <- vector("list", n_repeat)
for (repeat_idx in 1:n_repeat) {
  set.seed(repeat_idx)
  folds_by_repeat[[repeat_idx]] <- createFolds(group_labels, k = n_fold, list = TRUE, returnTrain = FALSE)
}

jobs <- expand.grid(repeat_idx = seq_len(n_repeat),
                    fold_idx   = seq_len(n_fold),
                    KEEP.OUT.ATTRS = FALSE,
                    stringsAsFactors = FALSE)

# Function for a single task: return output_data (data.frame) for the given (repeat, fold)

run_one_job <- function(job_i) {
  
  repeat_idx <- jobs$repeat_idx[job_i]
  fold_idx   <- jobs$fold_idx[job_i]
  folds      <- folds_by_repeat[[repeat_idx]]
  # --- progress logging (start) ---
  msg_start <- sprintf("[INFO] Starting repeat=%d fold=%d\n", repeat_idx, fold_idx)
  cat(msg_start); flush.console()

  test_index  <- folds[[fold_idx]]
  train_index <- setdiff(seq_len(nrow(data_2d)), test_index)
  train_data  <- data_2d[train_index, , drop = FALSE]
  test_data   <- data_2d[test_index,  , drop = FALSE]
  train_group <- group_labels[train_index]
  test_group  <- group_labels[test_index]

  # t-test filter (keep a feature if it is significant for any pairwise comparison)
  filter_by_t_test <- function(feature_data) {
    p_values <- c()
    for (i in 1:(length(group_levels) - 1)) {
      for (j in (i + 1):length(group_levels)) {
        g1 <- feature_data[train_group == group_levels[i]]
        g2 <- feature_data[train_group == group_levels[j]]
        test_result <- try(t.test(g1, g2), silent = TRUE)
        p_values <- c(p_values, ifelse(is(test_result, "try-error"), 1, test_result$p.value))
      }
    }
    any(p_values < 0.01)
  }
  filtered_feature <- names(which(apply(train_data, 2, filter_by_t_test)))
  cat(paste("Features retained after t-test screening (rep", repeat_idx, "fold", fold_idx, "):", length(filtered_feature), "\n"),
      file = paste0(out_dir, prefix, "_region_selection.txt"), append = TRUE)

  train_data <- train_data[, filtered_feature, drop = FALSE]
  test_data  <- test_data[,  filtered_feature, drop = FALSE]

  # Boruta feature selection
  set.seed(1)
  boruta_input <- train_data
  boruta_input$group <- as.factor(train_group)
  boruta_output <- Boruta(group ~ ., data = boruta_input, pValue = 0.01, mcAdj = TRUE, maxRuns = 500)
  boruta_result <- attStats(boruta_output)
  boruta_regions <- rownames(boruta_result[boruta_result$decision == 'Confirmed', ])
  
  cat(paste("Features retained after Boruta (rep", repeat_idx, "fold", fold_idx, "):", length(boruta_regions), "\n"),
      file = paste0(out_dir, prefix, "_region_selection.txt"), append = TRUE)
  line <- paste(paste0("rep", repeat_idx, "_fold", fold_idx),paste(rownames(train_data), collapse = ","),
                       paste(boruta_regions, collapse = ","),sep = "\t")
  cat(line, "\n",file = paste0(out_dir, prefix, "_region_markers.txt"),append = TRUE)
  
  # If Boruta selects nothing, fall back to the t-test filtered features
  if (length(boruta_regions) == 0L) {
    boruta_regions <- filtered_feature
  }
  train_data <- train_data[, boruta_regions, drop = FALSE]
  test_data  <- test_data[,  boruta_regions, drop = FALSE]

  # XGBoost
  train_data$group <- factor(train_group, levels = group_levels)
  test_data$group  <- factor(test_group,  levels = group_levels)
  x_train <- as.matrix(train_data[, -ncol(train_data), drop = FALSE])
  y_train <- as.numeric(train_data$group) - 1    ##### Encode factor labels to {0,1}: 0 = group_levels[1], 1 = group_levels[2];# XGBoost (binary:logistic) outputs P(label = 1).
  x_test  <- as.matrix(test_data[,  -ncol(test_data),  drop = FALSE])

  if (length(group_levels) == 2) {
    # ===== Binary classification =====
    set.seed(1)
    xgboost_model <- xgboost(
      data = x_train, label = y_train,
      eta = 0.01, max_depth = 3, gamma = 0, colsample_bytree = 0.8,
      nrounds = 200, objective = "binary:logistic",
      verbose = 0, nthread = 1
    )
    prob_pos <- predict(xgboost_model, newdata = x_test)  # NOTE: Under the 'binary:logistic' objective, predict() returns P(label = 1).
    # Since labels are encoded as {0,1} based on factor levels,
    # this corresponds to the probability of group_levels[2] (the positive class).
    prob_neg <- 1 - prob_pos
    negative_class <- group_levels[1]
    positive_class <- group_levels[2]
    pred_class <- ifelse(prob_pos >= 0.5, positive_class, negative_class)

    output_data <- data.frame(
      Repeat = repeat_idx,
      Fold   = fold_idx,
      Sample = rownames(test_data),
      Prediction = pred_class,
      check.names = FALSE
    )
    output_data[[paste0(positive_class, "_Prob")]] <- prob_pos
    output_data[[paste0(negative_class, "_Prob")]] <- prob_neg
    output_data$True_Label <- get_true(output_data$Sample)
    output_data <- output_data[, c("Repeat","Fold","Sample","True_Label","Prediction",
                                   paste0(group_levels, "_Prob"))]
    cat(sprintf("[INFO] Finished repeat=%d fold=%d (n=%d)\n", repeat_idx, fold_idx, nrow(output_data))); flush.console()
    return(output_data)

  } else {
    # ===== Multiclass classification =====
    set.seed(1)
    xgboost_model <- xgboost(
      data = x_train, label = y_train,
      objective = "multi:softprob", num_class = length(group_levels),
      eta = 0.01, max_depth = 3, gamma = 0, colsample_bytree = 0.8,
      nrounds = 200, verbose = 0, nthread = 1
    )
    prob_matrix <- matrix(predict(xgboost_model, newdata = x_test),
                          ncol = length(group_levels), byrow = TRUE)
    colnames(prob_matrix) <- group_levels
    pred_class <- colnames(prob_matrix)[max.col(prob_matrix)]

    output_data <- data.frame(
      Repeat = repeat_idx, Fold = fold_idx,
      Sample = rownames(test_data), Prediction = pred_class,
      check.names = FALSE
    )
    for (g in group_levels) output_data[[paste0(g,"_Prob")]] <- prob_matrix[, g]
    output_data$True_Label <- get_true(output_data$Sample)
    output_data <- output_data[, c("Repeat","Fold","Sample","True_Label","Prediction",
                                   paste0(group_levels, "_Prob"))]
    
    cat(sprintf("[INFO] Finished repeat=%d fold=%d (n=%d)\n", repeat_idx, fold_idx, nrow(output_data))); flush.console()
    return(output_data)
  }
}
# --- Parallel execution ---
all_outputs <- NULL
if (.Platform$OS.type == "windows") {
  # Windows: PSOCK
  cl <- parallel::makeCluster(n_thread)
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)

  parallel::clusterExport(
    cl,
    varlist = c("jobs","folds_by_repeat","data_2d","group_labels","group_levels",
                "out_dir","prefix","get_true","run_one_job"),
    envir = environment()
  )
  parallel::clusterEvalQ(cl, {
    suppressPackageStartupMessages({ library(Boruta); library(xgboost) })
    NULL
  })

  res_list <- tryCatch(
    parallel::parLapply(cl, seq_len(nrow(jobs)), run_one_job),
    error = function(e) {
      stop(sprintf(
        "[FATAL] Subtask failed (possibly due to resource limits). Please rerun with fewer threads by lowering n_thread (e.g., 4–8) and keeping xgboost nthread=1. Details: %s",
        conditionMessage(e)
      ))
    }
  )
  all_outputs <- dplyr::bind_rows(res_list)

} else {
  # Linux/macOS: fork
  res_list <- tryCatch(
    parallel::mclapply(seq_len(nrow(jobs)), run_one_job, mc.cores = n_thread),
    error = function(e) {
      stop(sprintf(
        "[FATAL] Subtask failed (possibly due to resource limits). Please rerun with fewer threads by lowering n_thread (e.g., 4–8) and keeping xgboost nthread=1. Details: %s",
        conditionMessage(e)
      ))
    }
  )
  all_outputs <- dplyr::bind_rows(res_list)
}

# Sort and write results
if (!is.null(all_outputs) && nrow(all_outputs) > 0) {
  all_outputs$Repeat <- as.integer(all_outputs$Repeat)
  all_outputs$Fold   <- as.integer(all_outputs$Fold)
  all_outputs <- dplyr::arrange(all_outputs, Repeat, Fold, Sample)
}
# ===== Integrity check (ensure every repeat × fold has results) =====
expected <- expand.grid(Repeat = seq_len(n_repeat),
                        Fold = seq_len(n_fold),
                        KEEP.OUT.ATTRS = FALSE,
                        stringsAsFactors = FALSE)
actual <- unique(all_outputs[, c("Repeat", "Fold")])
missing <- dplyr::anti_join(expected, actual, by = c("Repeat", "Fold"))

if (nrow(missing) > 0) {
  stop(sprintf(
    "[FATAL] Incomplete results detected! %d repeat/fold combinations missing.\nMissing combinations:\n%s",
    nrow(missing),
    paste(apply(missing, 1, function(x) paste0("Repeat=", x[1], ", Fold=", x[2])), collapse = "; ")
  ))
} else {
  cat("[INFO] All repeat/fold combinations processed successfully.\n")
}
write.table(all_outputs,file = result_file,sep = "	",row.names = FALSE,col.names = TRUE,quote = FALSE,append = FALSE)
cat(sprintf("[INFO] Prediction probabilities saved to: %s\n", result_file))
# ---------------- Evaluation ----------------
data <- read.table(result_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE,check.names = FALSE)
prob_cols <- grep("_Prob$", colnames(data), value = TRUE)
group_levels <- sub("_Prob$", "", prob_cols, perl = TRUE)
if (length(unique(group_labels)) == 2) {
  # ===== Binary evaluation =====
  pn <- list(pos = group_levels[2], neg = group_levels[1]) 
  # NOTE: Evaluation follows the same label convention as model training:
  # group_levels[1] is the negative/control class,
  # group_levels[2] is the positive disease class.
  # All binary metrics (Precision, Recall, F1, Sensitivity, etc.)
  # are computed with respect to the level-2 (positive) class.
  metrics_results <- data.frame(
    Repeat = integer(), AUC = numeric(), Accuracy = numeric(),
    Precision = numeric(), Recall = numeric(), F1 = numeric(),
    Sensitivity_at_90_Specificity = numeric()
  )

  for (rep in sort(unique(data$Repeat))) {
    df <- dplyr::filter(data, Repeat == rep)
    true_labels <- factor(df$True_Label, levels = c(pn$neg, pn$pos))
    pred_labels <- factor(df$Prediction, levels = c(pn$neg, pn$pos))
    prob_pos <- df[[paste0(pn$pos, "_Prob")]]
    roc_obj <- pROC::roc(true_labels, prob_pos, levels = c(pn$neg, pn$pos), direction = "<")
    auc_value <- as.numeric(pROC::auc(roc_obj))
    sens_90spec <- as.numeric(pROC::coords(roc_obj, x = 0.9, input = "specificity", ret = "sensitivity"))
    cm <- caret::confusionMatrix(pred_labels, true_labels, positive = pn$pos)
    acc <- as.numeric(cm$overall["Accuracy"])
    byc <- cm$byClass
    prec <- if ("Precision" %in% colnames(byc)) byc["Precision"] else byc["Pos Pred Value"]
    rec  <- if ("Recall"    %in% colnames(byc)) byc["Recall"]    else byc["Sensitivity"]
    f1   <- if ("F1" %in% colnames(byc)) byc["F1"] else if (prec + rec > 0) 2 * prec * rec / (prec + rec) else NA_real_

    metrics_results <- rbind(
      metrics_results,
      data.frame(Repeat = rep, AUC = auc_value, Accuracy = acc,
                 Precision = prec, Recall = rec, F1 = f1,
                 Sensitivity_at_90_Specificity = sens_90spec)
    )

  }

} else {
  # =====  Multiclass evaluation =====
  metrics_results <- data.frame(
    Repeat = integer(), Accuracy = numeric(), Precision = numeric(),
    Recall = numeric(), F1 = numeric()
  )

  for (rep in sort(unique(data$Repeat))) {
    df <- dplyr::filter(data, Repeat == rep)
    true_labels <- factor(df$True_Label)
    pred_labels <- factor(df$Prediction, levels = levels(true_labels))
    cm <- caret::confusionMatrix(pred_labels, true_labels)
    acc <- as.numeric(cm$overall["Accuracy"])
    precision <- mean(cm$byClass[, "Precision"], na.rm = TRUE)
    recall    <- mean(cm$byClass[, "Recall"], na.rm = TRUE)
    f1        <- mean(cm$byClass[, "F1"], na.rm = TRUE)

    metrics_results <- rbind(
      metrics_results,
      data.frame(Repeat = rep, Accuracy = acc, Precision = precision, Recall = recall, F1 = f1)
    )
  }
}
metrics_file <- paste0(out_dir, prefix, "_XGBoost_predict_metrics.txt")
write.table(metrics_results, file = metrics_file, sep = "\t",row.names = FALSE, quote = FALSE)
cat(sprintf("[INFO] Evaluation metrics saved to: %s\n", metrics_file))
