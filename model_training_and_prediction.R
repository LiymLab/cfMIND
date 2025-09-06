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
if (length(args) < 7) {
  stop("Usage: Rscript model_training_and_prediction_multiclass.R <manifest_file> <cut_off> <prefix> <out_dir> <n_repeat> <n_fold> <n_thread>\n", call. = FALSE)
}

manifest_file = args[1]
cut_off = args[2]
prefix   = args[3]
out_dir  = args[4]
n_repeat = as.numeric(args[5])
n_fold   = as.numeric(args[6])
n_thread = as.numeric(args[7])

# Ensure output directory exists
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
if (!grepl("/$", out_dir)) out_dir <- paste0(out_dir, "/")


# =====================================================================
# Function: data_process
# Purpose : Process input cvs files and build feature matrix
# =====================================================================
data_process = function(data_list,cut_off){
  #for (data in data_list){assign(data,drop_na(rename(read.csv(paste0(data_dir,data,".csv"),header = T,row.names = 1),M_0 = X0,M_0.25 = X0.25,M_0.5 = X0.5,M_0.75 =X0.75 ,M_1 = X1)))}
    for (data in data_list){assign(data,drop_na(rename(read.csv(paste0(data_dir,data,".csv"),header = T)[,-6],M_0 = region,M_0.25 = X0,M_0.5 = X0.25,M_0.75 = X0.5,M_1 = X0.75)))}
  #get the common regions from all data 
  regions = c()
  for (data in data_list){regions = append(regions,rownames(get(data)))}
  freq = as.data.frame(table(regions))
  regions = freq$regions[freq$Freq == length(data_list)]

  #get rid of the uncommon data
  regions = as.character(regions)
  for (data in data_list){assign(data,get(data)[regions,])}
  cat(paste("Common regions:", length(regions), "\n"),
      file = paste0(out_dir, prefix, "_region_selection.txt"), append = TRUE)  
  #get the average data from all data
  sum_data = data.frame(region = rownames(get(data)),M_0 = 0,M_0.25 = 0,M_0.5 = 0,M_0.75 = 0,M_1 = 0)
  rownames(sum_data) = sum_data$region
  sum_data = sum_data[,-1]
  for (data in data_list){
    sum_data$M_0 = sum_data$M_0 + get(data)$M_0
    sum_data$M_0.25 = sum_data$M_0.25 + get(data)$M_0.25
    sum_data$M_0.5 = sum_data$M_0.5 + get(data)$M_0.5
    sum_data$M_0.75 = sum_data$M_0.75 + get(data)$M_0.75
    sum_data$M_1 = sum_data$M_1 + get(data)$M_1
  }
  ave_data = sum_data/length(data_list)
  ave_data$SUM = apply(ave_data,1,sum)
  ave_data$SUM_ROUND = round(ave_data$SUM)
  
  #the distribution of region sum value from all data
  ggplot(data = ave_data,mapping = aes(x = SUM_ROUND,y = after_stat(prop),group = 1)) + 
    geom_bar() + 
    scale_x_continuous(limits = c(0, 100))
  
  #get rid of the low sum value data
  regions = rownames(ave_data[ave_data$SUM > cut_off,])
  for (data in data_list){assign(data,get(data)[regions,])}
  cat(paste("Regions retained after Coverage >", cut_off, "filter:", length(regions), "\n"),
      file = paste0(out_dir, prefix, "_region_selection.txt"), append = TRUE)   
  #normalize all data by their count per 10 million
  for (data in data_list){assign(data,get(data)/(sum(get(data)/10000000)))}
  
  #feature matrix
  data_2d = as.data.frame(matrix(nrow = length(data_list),ncol = length(regions)*5))
  turn_2d = function(region){paste0(region,c("_M_0","_M_0.25","_M_0.5","_M_0.75","_M_1"))}
  regions_2d = unlist(lapply(regions,turn_2d),)
  rownames(data_2d) = data_list
  colnames(data_2d) = regions_2d
  for (data in data_list){
    df = get(data)
    data_2d[data,] = c(t(df))
  }
  return(data_2d)
}
# ---------------------- Load Manifest ----------------------
manifest <- read.csv(manifest_file, header = TRUE, sep = "\t",
                     stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("path","group") %in% colnames(manifest))) {
  stop("Manifest must contain columns: path, group")
}
manifest$path  <- trimws(manifest$path)
manifest$group <- trimws(manifest$group)


manifest$Sample <- tools::file_path_sans_ext(basename(manifest$path))
sample_to_group <- setNames(trimws(manifest$group), trimws(manifest$Sample))

get_true <- function(samples) {
  s <- trimws(samples)
  g <- unname(sample_to_group[s])
  if (anyNA(g)) warning("Some samples did not match labels in manifest: ",
                        paste(unique(s[is.na(g)]), collapse = ", "))
  g
}
dirs <- unique(dirname(manifest$path))
if (length(dirs) != 1) {
  stop("All manifest$path must be in the same directory.")
}
data_dir  <- paste0(dirs, "/")
data_list <- tools::file_path_sans_ext(basename(manifest$path))
#run data_process
data_2d <- data_process(data_list = data_list, cut_off = cut_off)
Rdata <- file.path(out_dir, paste0(prefix, "_data_2d.Rdata"))
write.table(
  data_2d,
  file = paste0(out_dir,prefix, "_feature_matrix.txt"),
  sep = "\t",        
  row.names = TRUE,  
  col.names = NA,   
  quote = FALSE
)

# ---------------- Modeling & Prediction (Parallel) ----------------
group_labels <- manifest$group
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
  y_train <- as.numeric(train_data$group) - 1
  x_test  <- as.matrix(test_data[,  -ncol(test_data),  drop = FALSE])

  if (length(group_levels) == 2) {
    # ===== Binary classification =====
    negative_class <- group_levels[1]
    positive_class <- group_levels[2]
    xgboost_model <- xgboost(
      data = x_train, label = y_train,
      eta = 0.01, max_depth = 3, gamma = 0, colsample_bytree = 0.8,
      nrounds = 200, objective = "binary:logistic",
      verbose = 0, nthread = 1
    )
    prob_pos <- predict(xgboost_model, newdata = x_test)
    prob_neg <- 1 - prob_pos
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