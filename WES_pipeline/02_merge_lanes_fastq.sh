#!/usr/bin/env bash
# ==============================================================================
# 02_merge_lanes_fastq.sh
#
# Purpose
# -------
# Merge FASTQs across sequencing lanes (L001-L004) for each sample,
# separately for R1 and R2, producing one merged FASTQ per read:
#
#   Sample_X/..._L001_R1_001.fastq.gz
#   Sample_X/..._L002_R1_001.fastq.gz
#   ...
#   -> merged_fastq/Sample_X_R1.fastq.gz
#
# Same for R2.
#
# Inputs / assumptions
# --------------------
# - You are in a directory containing sample folders, e.g. Sample_*/...
# - FASTQs are gzipped (*.fastq.gz)
# - Lane is in filename as _L001_ _L002_ _L003_ _L004_
# - Read is in filename as _R1_ or _R2_
#
# Output
# ------
# - Creates: merged_fastq/
# - Writes merged FASTQs there:
#     <SampleDirName>_R1.fastq.gz
#     <SampleDirName>_R2.fastq.gz
#
# Safety / modes
# --------------
# DRY RUN by default (prints actions, does not execute).
# To execute merges:  RUN=1 bash 02_merge_lanes_fastq.sh
#
# Optional cleanup (off by default):
#   REMOVE_LANES=1  -> delete original lane-level fastqs after merge
#   REMOVE_EMPTY=1  -> remove empty leftover folders like *_qjbcuy, etc.
#
# Examples
# --------
# Preview:
#   bash WES_pipeline/02_merge_lanes_fastq.sh
#
# Run merges:
#   RUN=1 bash WES_pipeline/02_merge_lanes_fastq.sh
#
# Run merges + delete lane FASTQs:
#   RUN=1 REMOVE_LANES=1 bash WES_pipeline/02_merge_lanes_fastq.sh
#
# ==============================================================================

set -euo pipefail

RUN="${RUN:-0}"
REMOVE_LANES="${REMOVE_LANES:-0}"
REMOVE_EMPTY="${REMOVE_EMPTY:-0}"
OUTDIR="${OUTDIR:-merged_fastq}"

log() {
  printf "[%s] %s\n" "$(date '+%F %T')" "$*"
}

die() {
  log "ERROR: $*"
  exit 1
}

run_cmd() {
  if [[ "$RUN" -eq 1 ]]; then
    eval "$@"
  else
    echo "[DRYRUN] $*"
  fi
}

log "Merging lanes across FASTQs"
log "RUN=$RUN REMOVE_LANES=$REMOVE_LANES REMOVE_EMPTY=$REMOVE_EMPTY OUTDIR=$OUTDIR"

mkdir -p "$OUTDIR"

# ------------------------------------------------------------------------------
# Step 1: Identify sample directories that contain multi-lane FASTQs
# We detect samples with >= 2 lane files for R1 (or R2).
# ------------------------------------------------------------------------------
log "Scanning for sample directories with lane FASTQs..."

# Collect sample dirs that contain lane FASTQs
mapfile -t SAMPLE_DIRS < <(
  find . -maxdepth 1 -type d -name "Sample_*" -print |
  sort
)

[[ "${#SAMPLE_DIRS[@]}" -gt 0 ]] || die "No Sample_* directories found."

# Helper: gather lane fastqs for a sample + read
gather_fastqs() {
  local sample_dir="$1"
  local read="$2"   # R1 or R2

  # Search within sample_dir for lane fastqs for that read
  find "$sample_dir" -type f -name "*_L00[1-4]_${read}_*.fastq.gz" | sort
}

# ------------------------------------------------------------------------------
# Step 2: Merge per sample, per read
# ------------------------------------------------------------------------------
merged_any=0

for sd in "${SAMPLE_DIRS[@]}"; do
  sample_base="$(basename "$sd")"

  for read in R1 R2; do
    mapfile -t fq_list < <(gather_fastqs "$sd" "$read")

    # Need at least 2 lane FASTQs to merge
    if [[ "${#fq_list[@]}" -lt 2 ]]; then
      continue
    fi

    merged_any=1
    out_fq="${OUTDIR}/${sample_base}_${read}.fastq.gz"

    log "Will merge ${#fq_list[@]} files for ${sample_base} ${read} -> ${out_fq}"

    # Build a safe cat command
    # NOTE: cat gz files is valid; result is concatenated gzip stream.
    cmd="cat"
    for f in "${fq_list[@]}"; do
      cmd+=" \"${f}\""
    done
    cmd+=" > \"${out_fq}\""

    run_cmd "$cmd"
  done
done

if [[ "$merged_any" -eq 0 ]]; then
  log "No multi-lane samples detected (nothing to merge)."
  exit 0
fi

# ------------------------------------------------------------------------------
# Step 3 (Optional): Remove original lane-level FASTQs after successful merge
# ------------------------------------------------------------------------------
if [[ "$REMOVE_LANES" -eq 1 ]]; then
  log "REMOVE_LANES=1: removing original lane-level FASTQs (_L00[1-4]_R[12]_*.fastq.gz)"
  run_cmd "find . -type f -name '*_L00[1-4]_R[12]_*fastq.gz' -print -delete"
else
  log "REMOVE_LANES=0: keeping original lane-level FASTQs"
fi

# ------------------------------------------------------------------------------
# Step 4 (Optional): Remove empty leftover folders 
# command: rmdir *_[[:lower:]]*
# Here we apply it safely and only if empty.
# ------------------------------------------------------------------------------
if [[ "$REMOVE_EMPTY" -eq 1 ]]; then
  log "REMOVE_EMPTY=1: removing empty sample subfolders with lowercase suffixes"
  # This targets directories like Sample_..._qjbcuy or ..._bsfzgj if empty.
  # It will not remove non-empty directories.
  run_cmd "find . -maxdepth 1 -type d -name 'Sample_*_[a-z]*' -empty -print -delete"
else
  log "REMOVE_EMPTY=0: not removing any directories"
fi

log "Done."
log "Merged FASTQs are in: ${OUTDIR}/"
