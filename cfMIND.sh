#!/bin/bash

# ===== Configuration =====
script_path="$(readlink -f "$0")"
script_dir="$(dirname "$script_path")"

# Use the Python / R interpreter from the current environment
PYTHON_BIN="$(command -v python)"
RSCRIPT_BIN="$(command -v Rscript)"

# Script paths (located in the same directory as cfMIND.sh)
SCRIPT_MLE_DEFAULT="$script_dir/methylation_levels_encoding.py"
SCRIPT_DATA_PROCESS="$script_dir/feature_matrix.R"
SCRIPT_MODEL="$script_dir/model_training_and_prediction.R"                
# Built-in regions location
BUILTIN_HG38="$script_dir/hg38.500region3cpgs.bed"
BUILTIN_HG19="$script_dir/hg19.500region3cpgs.bed"
usage() {
  cat <<'EOF'

=========================================================================================
                                cfMIND Usage
=========================================================================================

Usage:
  bash cfMIND.sh feature_extraction [options]
  bash cfMIND.sh disease_detection  [options]

-----------------------------------------------------------------------------------------
 Options for feature_extraction
-----------------------------------------------------------------------------------------
  -m <file>   Manifest file (tab-delimited, with header). Must contain columns:
                - bam    : Input BAM file (sorted & indexed)
                - cpg_ob : CpG_OB* file (bottom strand) from bismark methylation extractor
                - cpg_ot : CpG_OT* file (top strand) from bismark methylation extractor
                - prefix : Output prefix for each sample
                - label  : Sample label (optional, for disease detection)
  -r <file>   BED file of genomic regions. Can use 'hg38' or 'hg19' as shortcuts,
              or provide a custom BED file path (4 tab-delimited columns, no header):
                <chromosome> <start> <end> <region_id>
                chr1    10000    10500    region_21
  -c <num>    Coverage cutoff threshold (default: 20)
  -w <int>    Window size of genomic regions in bp (default: 500)
  -l <int>    Number of methylation levels (default: 5)
  -p <str>    Output prefix for processed data (default: test)
  -o <dir>    Output directory (default: current working directory)
  -@ <int>    Threads for parallel processing (default: 1)
  --cfTAPS    Use this option for cfTAPS sequencing data
----------------------------------------------------------------------------------------
 Options for disease_detection
----------------------------------------------------------------------------------------
  -d <file>   Feature matrix file from feature_extraction step (e.g., feature_matrix.txt)
              This file should contain both feature matrix and labels from manifest
  -p <str>    Output prefix (default: test)
  -o <dir>    Output directory for prediction results (default: current working directory)
  -n <int>    Number of repeats for cross-validation (default: 1)
  -k <int>    Number of folds for cross-validation
              (default: number of samples in data, i.e. Leave-One-Out CV)
  -@ <int>    Threads for parallel cross-validation (default: 1)
========================================================================================

EOF
}

# ------------------------
# Part 1: feature_extraction
# ------------------------
feature_extraction() {
  local manifest="" bed_file="" cut_off="20" threads=1 cfTAPS=false
  local output_prefix="test" out_dir="." window_size=500 n_levels=5
  local script_mle="$SCRIPT_MLE_DEFAULT"

  while [[ $# -gt 0 ]]; do
    case "$1" in
      -m) manifest="$2"; shift 2 ;;
      -r) bed_file="$2"; shift 2 ;;
      -c) cut_off="$2"; shift 2 ;;
      -w) window_size="$2"; shift 2 ;;
      -l) n_levels="$2"; shift 2 ;;
      -p) output_prefix="$2"; shift 2 ;;
      -o) out_dir="$2"; shift 2 ;;
      -@) threads="$2"; shift 2 ;;
      --cfTAPS) cfTAPS=true; shift 1 ;;
      -h|--help) usage; exit 0 ;;
      *) echo "Unknown option: $1"; usage; exit 2 ;;
    esac
  done

  # Required options
  [[ -n "${manifest:-}" ]] || { echo "ERROR: -m <manifest> is required"; exit 2; }
  [[ -f "$manifest" ]] || { echo "ERROR: Manifest not found: $manifest"; exit 2; }
  [[ -n "${bed_file:-}" ]] || { echo "ERROR: -r <bed_file> is required"; exit 2; }
  
  # Resolve built-in shorthand to full path and validate BED file
  case "$(basename "$bed_file")" in
    "hg38") bed_file="$BUILTIN_HG38" ;;
    "hg19") bed_file="$BUILTIN_HG19" ;;
    *)
      # custom user path; check format
      if [[ ! -f "$bed_file" ]]; then
        echo "ERROR: BED file not found: $bed_file"; exit 2
      fi
      awk -F'\t' 'NF && NF!=4 {
        print "ERROR: BED file must have exactly 4 TAB-delimited columns."
        print "Expect: <chr> <start> <end> <region_id>"
        print "Example: chr1\t1\t500\tregion_1"
        exit 2
      }' "$bed_file" || exit 2
    ;;
  esac
  [[ -f "$bed_file" ]] || { echo "ERROR: BED file not found: $bed_file"; exit 2; }
  
  # Check manifest format
  header="$(head -n 1 "$manifest" | tr -d "\r")"
  if ! echo "$header" | awk -F'\t' 'BEGIN{IGNORECASE=1}
    {for(i=1;i<=NF;i++){if($i=="bam")b=1;if($i=="cpg_ob")o=1;if($i=="cpg_ot")t=1;if($i=="prefix")p=1;if($i=="label")l=1}}
    END{exit !(b&&o&&t&&p&&l)}'
  then
    echo "ERROR: Manifest header must contain columns: bam, cpg_ob, cpg_ot, prefix, label"; exit 2
  fi

  # Validate inputs
  [[ "$cut_off" =~ ^-?[0-9]+(\.[0-9]+)?$ ]] || { echo "ERROR: -c must be numeric"; exit 2; }
  [[ "$window_size" =~ ^[0-9]+$ && "$window_size" -ge 1 ]] || { echo "ERROR: -w must be a positive integer"; exit 2; }
  [[ "$n_levels" =~ ^[0-9]+$ && "$n_levels" -ge 2 ]] || { echo "ERROR: -l must be an integer >= 2"; exit 2; }
  awk -F'\t' -v w="$window_size" 'NF && ($3-$2)!=w {bad++} END{if (bad > 0) printf "WARNING: %d region(s) in BED file do not match window size (-w %s)\n", bad, w}' "$bed_file"
  [[ -x "$PYTHON_BIN" ]] || { echo "ERROR: python not found in PATH"; exit 127; }
  [[ -x "$RSCRIPT_BIN" ]] || { echo "ERROR: Rscript not found in PATH"; exit 127; }
  [[ -r "$script_mle" ]] || { echo "ERROR: Python script not readable: $script_mle"; exit 1; }
  [[ -r "$SCRIPT_DATA_PROCESS" ]] || { echo "ERROR: R script not readable: $SCRIPT_DATA_PROCESS"; exit 1; }
  
  # Create output directory if it doesn't exist
  mkdir -p "$out_dir"

  echo "Processing manifest file: $manifest"
  echo "BED file: $bed_file"
  echo "Coverage cutoff threshold: $cut_off"
  echo "Window size: $window_size"
  echo "Methylation levels: $n_levels"
  echo "Output prefix: $output_prefix"
  echo "Output directory: $out_dir"
  # Process each line in manifest (skip header)
  sed '1d' "$manifest" | while IFS=$'\t' read -r bam ob ot prefix rest; do
    # Skip empty lines
    echo $bam
    [[ -n "$bam" ]] || continue
    echo "Processing sample: $prefix"
    
    # Validate inputs for this sample
    [[ -f "$bam" ]] || { echo "ERROR: BAM file not found: $bam"; exit 2; }
    [[ -f "${bam}.bai" ]] || { echo "ERROR: BAM index (.bai) missing for: $bam"; exit 2; }
    [[ -f "$ob" ]] || { echo "ERROR: CpG_OB file not found: $ob"; exit 2; }
    [[ -f "$ot" ]] || { echo "ERROR: CpG_OT file not found: $ot"; exit 2; }
    [[ "$ob" == *CpG_OB* ]] || { echo "ERROR: CpG_OB file must contain 'CpG_OB' in its name: $ob"; exit 2; }
    [[ "$ot" == *CpG_OT* ]] || { echo "ERROR: CpG_OT file must contain 'CpG_OT' in its name: $ot"; exit 2; }
    
    # Build and run command for this sample
    cmd=("$PYTHON_BIN" "$script_mle" \
         -i "$bam" -r "$bed_file" -b "$ob" -t "$ot" \
         -p "$out_dir/$prefix" -@ "$threads" -l "$n_levels")
    if $cfTAPS; then
      cmd+=(--cfTAPS)
    fi
    
    echo "Running feature extraction for $prefix: ${cmd[*]}"
    "${cmd[@]}" || { echo "ERROR: Feature extraction failed for $prefix"; exit 1; }
  done
  
  # After all feature extractions are complete, run data processing
  echo "Feature extraction completed. Running data processing..."
  
  # CSV files are now in the output directory
  csv_dir="$out_dir"
  
  echo "Running data processing:"
  echo "  $RSCRIPT_BIN $SCRIPT_DATA_PROCESS $csv_dir $cut_off $output_prefix $out_dir $manifest $n_levels"

  "$RSCRIPT_BIN" "$SCRIPT_DATA_PROCESS" \
    "$csv_dir" "$cut_off" "$output_prefix" "$out_dir" "$manifest" "$n_levels" || { echo "ERROR: Data processing failed"; exit 1; }
    
  echo "Feature extraction and data processing completed successfully!"
}

# ------------------------
# Part 2: disease_detection
# ------------------------
disease_detection() {
  local feature_matrix_file="" prefix="test" out_dir="." n_repeat=1 n_fold=0
  local n_thread="1"    
  while getopts ":d:p:o:n:k:@:h" opt; do
    case "$opt" in
      d) feature_matrix_file="$OPTARG" ;;
      p) prefix="$OPTARG" ;;
      o) out_dir="$OPTARG" ;;
      n) n_repeat="$OPTARG" ;;
      k) n_fold="$OPTARG" ;;
      @) n_thread="$OPTARG" ;;
      h) usage; exit 0 ;;
      \?) echo "Unknown option: -$OPTARG"; usage; exit 2 ;;
      :)  echo "Option -$OPTARG requires an argument."; usage; exit 2 ;;
    esac
  done

  [[ -n "$feature_matrix_file" ]] || { echo "ERROR: -d <feature_matrix_file> is required"; exit 2; }
  [[ -f "$feature_matrix_file" ]] || { echo "ERROR: Feature matrix file not found: $feature_matrix_file"; exit 2; }

  # If n_fold is not specified, it will be determined from the data in the R script
  if [[ -z "${n_fold:-}" || "$n_fold" -eq 0 ]]; then
    n_fold=0  # Let R script determine from data
  fi

  # Validate inputs
  [[ "$n_repeat" =~ ^[0-9]+$ ]] || { echo "ERROR: -n must be integer"; exit 2; }
  [[ "$n_fold"   =~ ^[0-9]+$ ]] || { echo "ERROR: -k must be integer"; exit 2; }
  [[ "$n_thread" =~ ^[0-9]+$ ]] || { echo "ERROR: -@ must be integer"; exit 2; }
  [[ "$n_thread" -ge 1       ]] || { echo "ERROR: -@ must be >= 1"; exit 2; }

  [[ -x "$RSCRIPT_BIN" ]]  || { echo "ERROR: Rscript not found in PATH"; exit 127; }
  [[ -r "$SCRIPT_MODEL" ]] || { echo "ERROR: R script not readable: $SCRIPT_MODEL"; exit 1; }


  export OMP_NUM_THREADS=1
  export MKL_NUM_THREADS=1
  export OPENBLAS_NUM_THREADS=1

  echo "Running disease detection:"
  echo "$RSCRIPT_BIN $SCRIPT_MODEL \\"
  if [[ "$n_fold" -eq 0 ]]; then
    # Get sample count from feature matrix file (excluding header)
    sample_count=$(sed '1d' "$feature_matrix_file" | wc -l)
    echo "  $feature_matrix_file $prefix $out_dir $n_repeat $sample_count $n_thread"
  else
    echo "  $feature_matrix_file $prefix $out_dir $n_repeat $n_fold $n_thread"
  fi

  "$RSCRIPT_BIN" "$SCRIPT_MODEL" \
    "$feature_matrix_file" "$prefix" "$out_dir" "$n_repeat" "$n_fold" "$n_thread" || { echo "ERROR: Disease detection failed"; exit 1; }
}

main() {
  local cmd="${1:-}"
  shift || true
  case "$cmd" in
    feature_extraction) feature_extraction "$@" ;;
    disease_detection)  disease_detection  "$@" ;;
    ""|-h|--help)       usage ;;
    *) echo "Unknown command: $cmd"; usage; exit 2 ;;
  esac
}

main "$@"