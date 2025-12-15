#!/usr/bin/env bash
# ==============================================================================
# 01_fix_sample_names_and_paths.sh
#
# Purpose
# -------
# Normalize sample directory names after FASTQ download.
#
# Specifically:
#   1) Remove the substring "Sample_SV_" -> "Sample_"
#      to match naming conventions used by the rest of the cohort.
#
# Context
# -------
# Some sample directories were delivered with an extra "SV" tag
# (e.g. Sample_SV_10901987_PFC_BA10_375_202),
# while the rest of the dataset uses:
#   Sample_<ID>_<REGION>_<BA>_<I7>_<I5>
#
# Example:
#   Sample_SV_10901987_PFC_BA10_375_202_qjbcuy
#     -> Sample_10901987_PFC_BA10_375_202_qjbcuy
#
# This script ONLY renames directories; FASTQ filenames inside
# remain untouched at this stage.
#
# Usage
# -----
#   cd fastq_downloads
#   bash ../WES_pipeline/01_fix_sample_names_and_paths.sh
#
# Safety
# ------
# - Will NOT overwrite existing directories
# - Will print each rename before executing
#
# ==============================================================================

set -euo pipefail

log() {
  printf "[%s] %s\n" "$(date '+%F %T')" "$*"
}

log "Starting sample directory renaming (removing 'SV')"

# ------------------------------------------------------------------
# 1. Rename directories that contain 'Sample_SV_'
# ------------------------------------------------------------------
for dir in Sample_SV_*; do
  # Skip if glob doesn't match anything
  [[ -e "$dir" ]] || continue

  new_dir="${dir/Sample_SV_/Sample_}"

  if [[ -d "$new_dir" ]]; then
    log "WARNING: target directory already exists, skipping:"
    log "  $dir -> $new_dir"
    continue
  fi

  log "Renaming:"
  log "  $dir"
  log "  -> $new_dir"
  mv "$dir" "$new_dir"
done

log "Finished removing 'SV' from sample names"

# ------------------------------------------------------------------
# OPTIONAL STEP (commented out):
# Normalize accidental double underscores "__" -> "_"
#
# ------------------------------------------------------------------
#
# log "Normalizing double underscores '__' to '_'"
#
# for dir in Sample__*; do
#   [[ -e "$dir" ]] || continue
#
#   new_dir="$(echo "$dir" | sed 's/__/__/g; s/__/__/g')"
#   # If you actually want to collapse "__" -> "_", use:
#   # new_dir="$(echo "$dir" | sed 's/__*/_/g')"
#
#   if [[ "$dir" == "$new_dir" ]]; then
#     continue
#   fi
#
#   if [[ -d "$new_dir" ]]; then
#     log "WARNING: target directory already exists, skipping:"
#     log "  $dir -> $new_dir"
#     continue
#   fi
#
#   log "Renaming:"
#   log "  $dir"
#   log "  -> $new_dir"
#   mv "$dir" "$new_dir"
# done
#
# log "Finished optional underscore normalization"

log "All sample directory renaming complete."
