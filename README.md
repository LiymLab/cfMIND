# cfMIND: A read-level methylation framework for accurate non-invasive disease detection using cell-free DNA
cfMIND  (<em><u>cf</u></em>DNA <em><u>M</u></em>ethylation signals of <em><u>IN</u></em>dividual read for disease <em><u>D</u></em>etection)is a **machine learning-based framework** that enables identification of stratified cfDNA methylation signals at the **individual read level** for sensitive and robust disease detection. 

---


## 🌟Installation 
cfMIND requires both **Python** and **R** environments. For stable performance, we recommend the following versions:    

- **Python**: 3.7.12
- **R**: 4.3.2  

You can set up your conda environment and install the required R packages as follows:

```bash
git clone ***.git
chmod a+x -R ./cfMIND
cd cfMIND
conda env create -f environment.yml
conda activate cfMIND
Rscript install_R_packages.R
```

---

## 📂Data Preparation

To get started with **cfMIND**, please prepare the following files:

1. **BAM files**
   - Sorted by genomic coordinates and indexed (with`.bai`).
2. **CpG_OB/T files** corresponding to each BAM file:
   - **CpG_OB\***: methylation information for CpGs on the *original bottom strand (OB)*.
   
     **CpG_OT\***: methylation information for CpGs on the *original top strand (OT)*.

**Tip:** CpG_OB/T files can be generated using **Bismark methylation extractor** ([Bismark GitHub](https://github.com/FelixKrueger/Bismark)):

```bash
bismark_methylation_extractor sample.nsorted.bam -o /path/to/output/ --bedGraph --counts --no_overlap --genome_folder hg38
```
Once the installation and data preparation are complete, you are ready to start cfMIND. ✌️

## 🚀Quick Start

The main script **`cfMIND.sh`** supports **two major steps**:

1. **Feature extraction** – extract cfDNA methylation features from BAM and CpG files.
2. **Disease detection** – train and evaluate the model using leave-one-out (LOO) cross-validation.

### 1. Feature Extraction
This step extracts stratified cfDNA methylation signals from BAM and CpG files. Example command:

```bash
bash cfMIND.sh feature_extraction -i sample.csorted.bam -r hg38.500region3cpgs.bed -b CpG_OB_sample.txt -t CpG_OT_sample.txt -p /csvdir/sample -@ 8 
```
 Options:

```bash
-i: Input BAM file (sorted and indexed with .bai)
-r: BED file defining genomic regions of interest
-b: CpG_OB file (methylation calls on the original bottom strand)
-t: CpG_OT file (methylation calls on the original top strand)
-p: Output file prefix (normal samples must start with CTR)
-@: Number of threads (default: 8)
```
#### Output
cfMIND.sh feature_extraction produces `/csvdir/sample.csv`. Example content:
```txt
region,0,0.25,0.5,0.75,1
region_140,0,0,0,0,2
region_1734,0,0,0,1,2
region_1741,2,4,2,1,0
...
# Column definitions:
region: Genomic region ID
0: Number of reads with 0% methylation
0.25: Number of reads with ~25% methylation
0.5: Number of reads with ~50% methylation
0.75: Number of reads with ~75% methylation
1: Number of fully methylated reads (100%)
```
###  2. **Disease detection** 

This step uses the extracted features to train the model and perform disease detection using **leave-one-out (LOO) cross-validation**.(这里交叉验证的方式是否要修改)

```bash
bash cfMIND.sh disease_detection -d /csvdir/ -c 20 -p test -o /modeldir/ -@ 10
```

 Options:

```bash
-d: Directory containing feature files (*.csv) from feature extraction
-c: Coverage cutoff threshold for region (e.g., 20)
-p: Run prefix (used for output naming)
-o: Output directory for trained models and prediction results
-@: Number of threads for parallel LOO cross-validation
```

#### Output

- `*_predict_result.txt` : Sample predictions with probabilities.
- `*_predict_result.pdf` : ROC curve with evaluation metrics.
- `*_predict_metrics.txt` : Summary of AUC, accuracy, precision, recall, F1 score, and sensitivity at 90% specificity.
- `*region_length.txt` : Number of regions retained at each filtering step.

✅ With these two simple steps, you can efficiently extract read-level methylation features and perform **sensitive, non-invasive disease detection** using cfDNA.

