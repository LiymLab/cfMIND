# =============================================
# install_R_packages.R
# Install required R packages for cfMIND_full
# =============================================

# ---- 1. Set CRAN mirror (Tsinghua University) ----
#Note: For faster package downloads, it is recommended to select a CRAN mirror geographically close to you.
#This script uses the Tsinghua University mirror in China, which is optimal for users in China. 
#You may replace it with another mirror from https://cran.r-project.org/mirrors.html if needed.
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
cat(">>> Using CRAN mirror:", getOption("repos"), "\n")
# ---- Clean environment & set Conda compilers ----
Sys.unsetenv(c("CPLUS_INCLUDE_PATH", "C_INCLUDE_PATH", "CPATH", "LIBRARY_PATH"))

Sys.setenv(
  CC  = file.path(Sys.getenv("CONDA_PREFIX"), "bin", "x86_64-conda-linux-gnu-cc"),
  CXX = file.path(Sys.getenv("CONDA_PREFIX"), "bin", "x86_64-conda-linux-gnu-c++"),
  FC  = file.path(Sys.getenv("CONDA_PREFIX"), "bin", "x86_64-conda-linux-gnu-gfortran")
)
# ---- 2. Install common CRAN packages ----
cran_packages <- c(
  "dplyr",      # 1.1.4
  "tidyr",      # 1.3.1
  "stringr",    # 1.5.1
  "pROC",       # 1.19.0.1
  "caret",      # 7.0.1
  "ranger",
  "jsonlite"
)

# Check and install missing packages
to_install <- setdiff(cran_packages, rownames(installed.packages()))
if (length(to_install) > 0) {
  cat(">>> Installing missing CRAN packages:", paste(to_install, collapse = ", "), "\n")
  install.packages(to_install)
} else {
  cat(">>> All required CRAN packages are already installed.\n")
}

# ---- 3. Install archived versions ----
install.packages("https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.6-1.1.tar.gz",
                 repos = NULL, type = "source")
install.packages("https://cran.r-project.org/src/contrib/Archive/Rcpp/Rcpp_1.0.11.tar.gz",
                 repos = NULL, type = "source")
# Boruta 8.0.0
if (!requireNamespace("Boruta", quietly = TRUE)) {
  cat(">>> Installing Boruta (v8.0.0) from CRAN archive...\n")
  install.packages(
    "https://cran.r-project.org/src/contrib/Archive/Boruta/Boruta_8.0.0.tar.gz",
    repos = NULL,
    type = "source"
  )
} else {
  cat(">>> Boruta is already installed.\n")
}

# xgboost 1.7.7.1
if (!requireNamespace("xgboost", quietly = TRUE)) {
  cat(">>> Installing xgboost (v1.7.7.1) from CRAN archive...\n")
  install.packages(
    "https://cran.r-project.org/src/contrib/Archive/xgboost/xgboost_1.7.7.1.tar.gz",
    repos = NULL,
    type = "source"
  )
} else {
  cat(">>> xgboost is already installed.\n")
}
cat(">>> DONE. All required packages installed successfully.\n")
