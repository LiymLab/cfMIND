rm(list = ls())
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
})

# ---------------- Arguments ----------------
args = commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript feature_matrix.R <input_dir> <cut_off> <prefix> <out_dir> <manifest_file>\n", call. = FALSE)
}

input_dir = args[1]
cut_off = as.numeric(args[2])
prefix = args[3]
out_dir = args[4]
manifest_file = args[5]

# Ensure output directory exists
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
if (!grepl("/$", out_dir)) out_dir <- paste0(out_dir, "/")
if (!grepl("/$", input_dir)) input_dir <- paste0(input_dir, "/")

# =====================================================================
# Function: data_process
# Purpose : Process input csv files and build feature matrix
# =====================================================================
data_process = function(data_list, data_dir, cut_off, out_dir, prefix){
  for (data in data_list){
    assign(data, drop_na(rename(read.csv(paste0(data_dir, data, ".csv"), header = T, row.names = 1),
                                M_0 = X0, M_0.25 = X0.25, M_0.5 = X0.5, M_0.75 = X0.75, M_1 = X1)))
  }
  
  # Get the common regions from all data 
  regions = c()
  for (data in data_list){
    regions = append(regions, rownames(get(data)))
  }
  freq = as.data.frame(table(regions))
  regions = freq$regions[freq$Freq == length(data_list)]

  # Get rid of the uncommon data
  regions = as.character(regions)
  for (data in data_list){
    assign(data, get(data)[regions,])
  }
  cat(paste("Common regions:", length(regions), "\n"),
      file = paste0(out_dir, prefix, "_region_selection.log"), append = TRUE)  
  
  # Get the average data from all data
  sum_data = data.frame(region = rownames(get(data_list[1])), M_0 = 0, M_0.25 = 0, M_0.5 = 0, M_0.75 = 0, M_1 = 0)
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
  ave_data$SUM = apply(ave_data, 1, sum)
  ave_data$SUM_ROUND = round(ave_data$SUM)
  
  # The distribution of region sum value from all data
  p <- ggplot(data = ave_data, mapping = aes(x = SUM_ROUND, y = after_stat(prop), group = 1)) + 
    geom_bar() + 
    scale_x_continuous(limits = c(0, 100)) +
    labs(title = "Distribution of Region Coverage", x = "Coverage", y = "Proportion")
  ggsave(paste0(out_dir, prefix, "_coverage_distribution.png"), plot = p, width = 8, height = 6)
  
  # Get rid of the low sum value data
  regions = rownames(ave_data[ave_data$SUM > cut_off,])
  for (data in data_list){
    assign(data, get(data)[regions,])
  }
  cat(paste("Regions retained after Coverage >", cut_off, "filter:", length(regions), "\n"),
      file = paste0(out_dir, prefix, "_region_selection.log"), append = TRUE)   
  
  # Normalize all data by their count per 10 million
  for (data in data_list){
    assign(data, get(data)/(sum(get(data)/10000000)))
  }
  
  # Feature matrix
  data_2d = as.data.frame(matrix(nrow = length(data_list), ncol = length(regions)*5))
  turn_2d = function(region){paste0(region, c("_M_0", "_M_0.25", "_M_0.5", "_M_0.75", "_M_1"))}
  regions_2d = unlist(lapply(regions, turn_2d))
  rownames(data_2d) = data_list
  colnames(data_2d) = regions_2d
  for (data in data_list){
    df = get(data)
    data_2d[data,] = c(t(df))
  }
  return(data_2d)
}

# ---------------------- Main Processing ----------------------
# Get list of CSV files in input directory
csv_files <- list.files(input_dir, pattern = "\\.csv$", full.names = FALSE)
if (length(csv_files) == 0) {
  stop("No CSV files found in input directory: ", input_dir)
}

# Extract sample names (remove .csv extension)
data_list <- tools::file_path_sans_ext(csv_files)

cat(paste("Found", length(data_list), "CSV files in", input_dir, "\n"))
cat(paste("Processing with coverage cutoff:", cut_off, "\n"))

# Run data processing
data_2d <- data_process(data_list = data_list, data_dir = input_dir, 
                       cut_off = cut_off, out_dir = out_dir, prefix = prefix)

# ---------------------- Process Labels from Manifest ----------------------
# Read manifest file to get labels
if (file.exists(manifest_file)) {
  manifest <- read.table(manifest_file, header = TRUE, sep = "\t", 
                        stringsAsFactors = FALSE, check.names = FALSE)
  
  # Check if manifest has required columns including label
  required_cols <- c("bam", "cpg_ob", "cpg_ot", "prefix", "label")
  if (!all(required_cols %in% colnames(manifest))) {
    stop("Manifest file must contain columns: ", paste(required_cols, collapse = ", "))
  }
  
  # Create sample to label mapping
  sample_to_label <- setNames(manifest$label, manifest$prefix)
  
  # Get labels for samples in data_2d (in the same order)
  sample_labels <- sample_to_label[rownames(data_2d)]
  
  # Check for missing labels
  missing_labels <- is.na(sample_labels)
  if (any(missing_labels)) {
    warning("Missing labels for samples: ", 
            paste(rownames(data_2d)[missing_labels], collapse = ", "))
  }
  
  # Add labels to the data structure
  labels_df <- data.frame(
    sample = rownames(data_2d),
    group = sample_labels,
    stringsAsFactors = FALSE
  )
  
  cat(paste("Processed labels for", sum(!is.na(sample_labels)), "samples\n"))
} else {
  labels_df <- NULL
}

# Prepare final output with label column
if (!is.null(labels_df)) {
  # Add label column to feature matrix
  data_2d$label <- labels_df$group[match(rownames(data_2d), labels_df$sample)]
} else {
  # Add empty label column if no labels available
  data_2d$label <- NA
}

# Save results
write.table(
  data_2d,
  file = paste0(out_dir, prefix, "_feature_matrix.txt"),
  sep = "\t",        
  row.names = TRUE,  
  col.names = NA,   
  quote = FALSE
)

cat(paste("Data processing completed. Results saved with prefix:", prefix, "\n"))
cat(paste("Feature matrix saved to:", paste0(out_dir, prefix, "_feature_matrix.txt"), "\n"))
cat(paste("Final feature matrix dimensions:", nrow(data_2d), "samples x", ncol(data_2d), "features\n"))