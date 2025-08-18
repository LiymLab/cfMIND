#!/bin/bash
set -euo pipefail

# ===== Configuration =====
script_path="$(readlink -f "$0")"
script_dir="$(dirname "$script_path")"

# Use the Python / R interpreter from the current environment
PYTHON_BIN="$(command -v python)"
RSCRIPT_BIN="$(command -v Rscript)"

# Script paths (located in the same directory as cfMIND.sh)
SCRIPT_5MLE="$script_dir/5_methylation_levels_encoding.py"   # Part 1: Feature extraction (Python script)
SCRIPT_MODEL="$script_dir/model_training_and_prediction.R"   # Part 2: Disease detection (R script)

usage() {
  cat <<'EOF'
Usage:
  bash cfMIND.sh feature_extraction [options]
  bash cfMIND.sh disease_detection  [options]

Options for feature_extraction:
  -i <file>   Input BAM file (sorted & indexed)
  -r <file>   BED file of genomic regions
  -b <file>   CpG_OB* file (bottom strand) from bismark methylation extractor
  -t <file>   CpG_OT* file (top strand)  from bismark methylation extractor
  -p <str>    Output file prefix
  -@ <int>    Threads (default: 8)

Example:
  bash cfMIND.sh feature_extraction \
    -i sample.bam \
    -r regions.bed \
    -b CpG_OB_sample.txt.gz \
    -t CpG_OT_sample.txt.gz \
    -p output/sample \
    -@ 8
Note:
  For normal samples, the prefix (-p) must start with 'CTR'.

Options for disease_detection:
  -d <dir>    Directory containing *.csv feature files (from feature_extraction)
  -c <num>    Coverage cutoff threshold (numeric, e.g., 20)
  -p <str>    Prefix for outputs (e.g., run1)
  -o <dir>    Output directory for disease prediction results
  -@ <int>    Threads for parallel leave-one-out (LOO) cross-validation

Example:
  bash cfMIND.sh disease_detection \
    -d /path/to/csv_dir \
    -c 50 \
    -p liver \
    -o /path/to/model_out \
    -@ 8
EOF
}

# ------------------------
# Part 1: feature_extraction
# ------------------------
feature_extraction() {
  local bam="" bed="" ob="" ot="" prefix="" threads=8
  while getopts ":i:r:b:t:p:@:h" opt; do
    case "$opt" in
      i) bam="$OPTARG" ;;
      r) bed="$OPTARG" ;;
      b) ob="$OPTARG" ;;
      t) ot="$OPTARG" ;;
      p) prefix="$OPTARG" ;;
      @) threads="$OPTARG" ;;
      h) usage; exit 0 ;;
      \?) echo "Unknown option: -$OPTARG"; usage; exit 2 ;;
      :)  echo "Option -$OPTARG requires an argument."; usage; exit 2 ;;
    esac
  done

  # Validate inputs
  [[ -n "$bam" && -f "$bam" ]]       || { echo "ERROR: BAM file is missing"; exit 2; }
  [[ -f "${bam}.bai" ]]              || { echo "ERROR: BAM index (.bai) is missing"; exit 2; }
  [[ -n "$bed" && -f "$bed" ]]       || { echo "ERROR: BED file is missing"; exit 2; }
  [[ -n "$ob" && -f "$ob" ]]         || { echo "ERROR: CpG_OB file is missing"; exit 2; }
  [[ -n "$ot" && -f "$ot" ]]         || { echo "ERROR: CpG_OT file is missing"; exit 2; }
  [[ "$ob" == *CpG_OB* ]]            || { echo "ERROR: -b file name must contain CpG_OB"; exit 2; }
  [[ "$ot" == *CpG_OT* ]]            || { echo "ERROR: -t file name must contain CpG_OT"; exit 2; }
  [[ -n "$prefix" ]]                 || { echo "ERROR: Output prefix is missing"; exit 2; }

  # Execute the Python feature extraction script
  "$PYTHON_BIN" "$SCRIPT_5MLE" \
    -i "$bam" \
    -r "$bed" \
    -b "$ob" \
    -t "$ot" \
    -p "$prefix" \
    -@ "$threads"
}

# ------------------------
# Part 2: disease_detection
# ------------------------
disease_detection() {
  local data_dir="" cut_off="" prefix="" out_dir="" threads=""
  while getopts ":d:c:p:o:@:h" opt; do
    case "$opt" in
      d) data_dir="$OPTARG" ;;
      c) cut_off="$OPTARG" ;;
      p) prefix="$OPTARG" ;;
      o) out_dir="$OPTARG" ;;
      @) threads="$OPTARG" ;;
      h) usage; exit 0 ;;
      \?) echo "Unknown option: -$OPTARG"; usage; exit 2 ;;
      :)  echo "Option -$OPTARG requires an argument."; usage; exit 2 ;;
    esac
  done

  # Validate inputs
  [[ -n "$data_dir" && -d "$data_dir" ]] || { echo "ERROR: -d <data_dir> directory not found"; exit 2; }
  [[ -n "$cut_off" ]]                    || { echo "ERROR: -c <cut_off> is missing"; exit 2; }
  [[ "$cut_off" =~ ^-?[0-9]+(\.[0-9]+)?$ ]] || { echo "ERROR: -c <cut_off> must be numeric"; exit 2; }
  [[ -n "$prefix" ]]                     || { echo "ERROR: -p <prefix> is missing"; exit 2; }
  [[ -n "$out_dir" ]]                    || { echo "ERROR: -o <out_dir> is missing"; exit 2; }
  [[ -n "$threads" && "$threads" =~ ^[0-9]+$ ]] || { echo "ERROR: -@ <threads> must be an integer"; exit 2; }
  [[ -x "$RSCRIPT_BIN" ]]                || { echo "ERROR: Rscript not found in PATH"; exit 127; }
  [[ -r "$SCRIPT_MODEL" ]]               || { echo "ERROR: R model script is not readable: $SCRIPT_MODEL"; exit 1; }

  "$RSCRIPT_BIN" "$SCRIPT_MODEL" \
    "$data_dir" \
    "$cut_off" \
    "$prefix" \
    "$out_dir" \
    "$threads"
}

main() {
  local cmd="${1:-}"
  shift || true
  case "$cmd" in
    feature_extraction) feature_extraction "$@" ;;
    disease_detection)  disease_detection "$@" ;;
    ""|-h|--help)       usage ;;
    *) echo "Unknown command: $cmd"; usage; exit 2 ;;
  esac
}

main "$@"
