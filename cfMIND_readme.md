## cfMIND: A read-level methylation framework for accurate non-invasive disease detection using cell-free DNA
cfMIND  (<em><u>cf</u></em>DNA <em><u>M</u></em>ethylation signals of <em><u>IN</u></em>dividual read for disease <em><u>D</u></em>etection)is a **machine learning-based framework** that enables identification of stratified cfDNA methylation signals at the **individual read level** for sensitive and robust disease detection. 

---


## 🌟Installation 
cfMIND requires both **Python** and **R** environments. For stable performance, we recommend:

- **Python**: 3.7.12 
- **R**: 4.3.2

You can set up your conda environment and install the required R packages as follows:

```bash
git clone https://github.com/LiymLab/cfMIND.git
chmod a+x -R ./cfMIND
cd cfMIND
conda env create -f environment.yml
conda activate cfMIND
Rscript install_packages.R
```

---

## 📂Data Preparation

To get started with **cfMIND**, please prepare the following files:

1. **BAM files**
   - Coordinate-sorted and indexed (`.bam` with corresponding `.bai`).
2. **CpG_OB/T files** corresponding to each BAM file:
   - **CpG_OB\***: methylation information for CpGs on the *original bottom strand (OB)*.
   
     **CpG_OT\***: methylation information for CpGs on the *original top strand (OT)*.

**Tip:** Generate CpG_OB/OT with **Bismark methylation extractor**: ([Bismark GitHub](https://github.com/FelixKrueger/Bismark)):

```bash
bismark_methylation_extractor sample.nsorted.bam -o /path/to/output/ --bedGraph --counts --no_overlap --genome_folder /path/to/hg38
```
Once the installation and data preparation are complete, you are ready to start cfMIND. ✌️

## 🚀Quick Start

The main script **`cfMIND.sh`** supports **two major steps**:

1. **Feature extraction** – extract cfDNA methylation features from BAM and CpG files.

2. **Disease detection** – train and evaluate the model using cross-validation.

To help you verify your environment setup and quickly get started with cfMIND, we provide **example files** on Zenodo: https://doi.org/10.5281/zenodo.17067679.

You can run cfMIND with:

```bash
bash cfMIND.sh feature_extraction -i sample1.bam -r hg38 -b CpG_OB_sample1.txt -t CpG_OT_sample1.txt
bash cfMIND.sh disease_detection -m /path/to/manifest.txt
```

### 1. Feature Extraction
This step extracts stratified cfDNA methylation signals from BAM and CpG files. Example command:

```bash
bash cfMIND.sh feature_extraction -i sample.csorted.bam -r hg38 -b CpG_OB_sample.txt -t CpG_OT_sample.txt -p /csvdir/sample -@ 10
```
 Options:

```bash
-i: Input BAM file (sorted and indexed)
-r: BED file of genomic regions. Options: hg38/hg19, 
             or a custom BED file (4 tab-delimited columns, no header):
             Example: <chromosome> <start> <end> <region_id>
                       chr1       10000    10500  region_21
-b: CpG_OB file (methylation calls on the original bottom strand)
-t: CpG_OT file (methylation calls on the original top strand)
-p: Output prefix (path + filename prefix, default: sample).
    If no path is specified, outputs are written to the current directory.
-@: Threads for parallel processing (default: 1)
--cfTAPS: Optional flag if using cfTAPS sequencing data.
```
#### Output
cfMIND.sh feature_extraction produces `/csvdir/sample.csv`. **Columns:**

- **region** – Genomic region ID
- **0, 0.25, 0.5, 0.75, 1** – Read counts at each methylation level (0%, 25%, 50%, 75%, 100%)  

Example:

```txt
region,0,0.25,0.5,0.75,1
region_21,0,0,0,0,2
region_1734,0,0,0,1,2
region_1741,2,4,2,1,0
```
###  2. **Disease detection** 

This step uses the extracted features to train and evaluate the model with cross-validation. It requires a manifest file (e. g.,`manifest.txt`), which is a **tab-delimited text file** containing two columns:

- **`path`** – Full path to each sample’s `.csv` file produced by the feature_extraction step.
- **`group`** – Sample label (e. g., `CTR` for control, `Tumor` for cancer, or other custom labels).

```txt
path	group
/csvdir/sample1.csv	CTR
/csvdir/sample2.csv	CTR
/csvdir/sample3.csv	Tumor
/csvdir/sample4.csv	Tumor
```

Example command:

```bash
bash cfMIND.sh disease_detection -m /path/to/manifest.txt -c 20 -p test -o /outdir/  -n 10 -k 5 -@ 5
```

 Options:

```bash
-m: Manifest file (tab-delimited) with two columns:
   path – Full path to each .csv file
   group – sample label (e.g., CTR, tumor, or other multi-class tags)
-c: Coverage cutoff threshold for region (default: 20)
-p: Prefix for prediction results (default: test)
-o: Output directory for prediction results (default: current working directory)
-n: Number of repeats for cross-validation (default: 1)
-k: Number of folds for cross-validation (default: number of samples in manifest file, i.e. Leave-One-Out CV)
-@: Threads for parallel cross-validation (default: 1)
```

#### Output

- `<prefix>_XGBoost_predict_probability.txt` : Sample predictions with probabilities.

- `<prefix>_XGBoost_predict_metrics.txt` : Summary of evaluation metrics for each repeat.

  - For **binary classification**: AUC, accuracy, precision, recall, F1 score, and sensitivity at 90% specificity.
  - For **multi-class classification**: accuracy, precision, recall, and F1 score.

  Example:

  ```txt
  Repeat	AUC	Accuracy	Precision	Recall	F1	Sensitivity_at_90_Specificity
  1	0.888	0.853	0.894	0.790	0.839	0.790
  2	0.909	0.831	0.804	0.860	0.831	0.744
  ```

✅ With only two streamlined steps, **cfMIND** empowers users to capture read-level methylation features and achieve  sensitive, non-invasive disease detection.



