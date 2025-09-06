#!/bin/bash
set -euo pipefail

# ===== Configuration =====
script_path="$(readlink -f "$0")"
script_dir="$(dirname "$script_path")"

# Use the Python / R interpreter from the current environment
PYTHON_BIN="$(command -v python)"
RSCRIPT_BIN="$(command -v Rscript)"

# Script paths (located in the same directory as cfMIND.sh)
SCRIPT_5MLE_DEFAULT="$script_dir/5_methylation_levels_encoding.py"      
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
  -i <file>   Input BAM file (sorted & indexed)                     (required)
  -r <file>   BED file of genomic regions. Options: hg38/hg19       (required)
              or the path to a custom BED file (4 tab-delimited columns, no header):
                <chromosome> <start> <end> <region_id>
                chr1    10000    10500    region_21
  -b <file>   CpG_OB* file (bottom strand) from bismark methylation extractor (required)
  -t <file>   CpG_OT* file (top strand)  from bismark methylation extractor   (required)
  -p <str>    Output prefix (path + filename prefix, default: sample)
              If no path is specified, outputs are written to the current directory
  -@ <int>    Threads for parallel processing (default: 1)
  --cfTAPS    Use this option for cfTAPS sequencing data
----------------------------------------------------------------------------------------
 Options for disease_detection
----------------------------------------------------------------------------------------
  -m <file>   Manifest file (tab-delimited, with header). Must contain columns:
                - path : .csv file from feature_extraction
                - group: sample label (e.g., CTR, tumor, or other multi-class tags) (required)
  -c <num>    Coverage cutoff threshold (default: 20)
  -p <str>    Output prefix (default: test)
  -o <dir>    Output directory for prediction results (default: current working directory)
  -n <int>    Number of repeats for cross-validation (default: 1)
  -k <int>    Number of folds for cross-validation
              (default: number of samples in manifest file, i.e. Leave-One-Out CV)
  -@ <int>    Threads for parallel cross-validation (default: 1)
========================================================================================

EOF
}


# ------------------------
# Part 1: feature_extraction
# ------------------------
feature_extraction() {
  local bam="" bed="" ob="" ot="" prefix="sample" threads=1 cfTAPS=false
  local script_5mle="$SCRIPT_5MLE_DEFAULT"

  while [[ $# -gt 0 ]]; do
    case "$1" in
      -i) bam="$2"; shift 2 ;;
      -r) bed="$2"; shift 2 ;;
      -b) ob="$2"; shift 2 ;;
      -t) ot="$2"; shift 2 ;;
      -p) prefix="$2"; shift 2 ;;
      -@) threads="$2"; shift 2 ;;
      --cfTAPS) cfTAPS=true; shift 1 ;;
      -h|--help) usage; exit 0 ;;
      *) echo "Unknown option: $1"; usage; exit 2 ;;
    esac
  done

  # Required options
  [[ -n "${bam:-}"    ]] || { echo "ERROR: -i <bam> is required"; exit 2; }
  [[ -n "${bed:-}"    ]] || { echo "ERROR: -r <bed> is required"; exit 2; }
  [[ -n "${ob:-}"     ]] || { echo "ERROR: -b <CpG_OB> is required"; exit 2; }
  [[ -n "${ot:-}"     ]] || { echo "ERROR: -t <CpG_OT> is required"; exit 2; }

  # Resolve built-in shorthand to full path (same directory as script)
  case "$(basename "$bed")" in
    "hg38") bed="$BUILTIN_HG38" ;;
    "hg19") bed="$BUILTIN_HG19" ;;
    *)
      # custom user path; leave as-is
      # check format: 4 tab-delimited columns, no header
      if [[ ! -f "$bed" ]]; then
        echo "ERROR: BED file not found: $bed"; exit 2
      fi
      awk -F'\t' 'NF && NF!=4 {
        print "ERROR: BED file must have exactly 4 TAB-delimited columns."
        print "Expect: <chr> <start> <end> <region_id>"
        print "Example: chr1\t1\t500\tregion_1"
        exit 2
      }' "$bed" || exit 2
    ;;
esac

  # Validate inputs
  [[ -f "$bam" ]]                    || { echo "ERROR: BAM file not found: $bam"; exit 2; }
  [[ -f "${bam}.bai" ]]              || { echo "ERROR: BAM index (.bai) missing for: $bam"; exit 2; }
  [[ -f "$bed" ]]                    || { echo "ERROR: BED file not found: $bed"; exit 2; }
  [[ -f "$ob"  ]]                    || { echo "ERROR: CpG_OB file not found: $ob"; exit 2; }
  [[ -f "$ot"  ]]                    || { echo "ERROR: CpG_OT file not found: $ot"; exit 2; }
  [[ "$ob" == *CpG_OB* ]]            || { echo "ERROR: -b file must contain 'CpG_OB' in its name"; exit 2; }
  [[ "$ot" == *CpG_OT* ]]            || { echo "ERROR: -t file must contain 'CpG_OT' in its name"; exit 2; }
  [[ -x "$PYTHON_BIN" ]]             || { echo "ERROR: python not found in PATH"; exit 127; }
  [[ -r "$script_5mle" ]]            || { echo "ERROR: Python script not readable: $script_5mle"; exit 1; }

  # Build and run command
  cmd=("$PYTHON_BIN" "$script_5mle" \
       -i "$bam" -r "$bed" -b "$ob" -t "$ot" \
       -p "$prefix" -@ "$threads")
  if $cfTAPS; then
    cmd+=(--cfTAPS)
  fi

  echo "Running feature extraction: ${cmd[*]}"
  "${cmd[@]}"
}

# ------------------------
# Part 2: disease_detection
# ------------------------
disease_detection() {
  local manifest="" cut_off="20" prefix="test" out_dir="." n_repeat=1 n_fold=0
  local n_thread="1"    
  while getopts ":m:c:p:o:n:k:@:h" opt; do
    case "$opt" in
      m) manifest="$OPTARG" ;;
      c) cut_off="$OPTARG" ;;
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

  [[ -n "$manifest" ]] || { echo "ERROR: -m <manifest> is required"; exit 2; }
  [[ -f "$manifest" ]] || { echo "ERROR: Manifest not found: $manifest"; exit 2; }
  # check format
  header="$(head -n 1 "$manifest" | tr -d "\r")"
  if ! echo "$header" | awk -F'\t' 'BEGIN{IGNORECASE=1}
    {for(i=1;i<=NF;i++){if($i=="path")p=1;if($i=="group")g=1}}
    END{exit !(p&&g)}'
  then
    echo "ERROR: Manifest header must contain columns: path, group"; exit 2
  fi

  if [[ -z "${n_fold:-}" || "$n_fold" -eq 0 ]]; then
    n_fold="$(tail -n +2 "$manifest" | grep -v '^[[:space:]]*$' | wc -l | tr -d ' ')"
  fi

  # Validate inputs
  [[ "$cut_off"  =~ ^-?[0-9]+(\.[0-9]+)?$ ]] || { echo "ERROR: -c must be numeric"; exit 2; }
  [[ "$n_repeat" =~ ^[0-9]+$               ]] || { echo "ERROR: -n must be integer"; exit 2; }
  [[ "$n_fold"   =~ ^[0-9]+$               ]] || { echo "ERROR: -k must be integer"; exit 2; }
  [[ "$n_thread" =~ ^[0-9]+$               ]] || { echo "ERROR: -@ must be integer"; exit 2; }
  [[ "$n_thread" -ge 1                     ]] || { echo "ERROR: -@ must be >= 1"; exit 2; }

  [[ -x "$RSCRIPT_BIN" ]]  || { echo "ERROR: Rscript not found in PATH"; exit 127; }
  [[ -r "$SCRIPT_MODEL" ]] || { echo "ERROR: R script not readable: $SCRIPT_MODEL"; exit 1; }


  export OMP_NUM_THREADS=1
  export MKL_NUM_THREADS=1
  export OPENBLAS_NUM_THREADS=1

  echo "Running disease detection:"
  echo "  $RSCRIPT_BIN $SCRIPT_MODEL \\"
  echo "    $manifest $cut_off $prefix $out_dir $n_repeat $n_fold $n_thread"

  "$RSCRIPT_BIN" "$SCRIPT_MODEL" \
    "$manifest" "$cut_off" "$prefix" "$out_dir" "$n_repeat" "$n_fold" "$n_thread"
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