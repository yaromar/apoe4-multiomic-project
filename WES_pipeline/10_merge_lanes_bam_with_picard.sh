#!/bin/bash
# 10_merge_lanes_bam_with_picard.sh
#
# Purpose
# -------
# Merge multiple lane-level (already coordinate-sorted) BAMs within each sample folder
# into a single *_merged.bam using Picard MergeSamFiles, then delete the inputs + indices.
#
# Assumptions
# -----------
# - Input BAMs are already coordinate-sorted (e.g., produced by 09_sort_bam_with_picard.sh)
# - Each input BAM has a corresponding .bai in the same folder
# - Lanes live in the same sample folder under ~/scratch60/*
#
# DSQ usage
# ---------
# This script PRINTS one command per folder that needs merging.
#   bash 10_merge_lanes_bam_with_picard.sh > joblist.txt
#   module load dSQ
#   dsq --job-file joblist.txt -c 1 --mem-per-cpu=5GB -t 3:00:00 --max-jobs 200 --output DSQ_results
#
# Slurm parameters (already validated)
# -----------------------------------
# DSQ-generated .slurm should include something like:
#   #SBATCH --array 0-208%200
#   #SBATCH --output DSQ_results
#   #SBATCH --job-name dsq-joblist
#   #SBATCH -c 1 --mem-per-cpu=5GB -t 3:00:00
#
# Output naming
# -------------
# Output is placed in the same folder and named:
#   <prefix>_merged.bam
# where <prefix> is derived from the first BAM's basename, taking the first two
# underscore-delimited fields (kept consistent with your existing convention).

for folder in ~/scratch60/*; do
  # Count BAMs (ignore already-merged BAMs if present)
  bamCount=$(ls -1q "$folder"/*bam 2>/dev/null | grep -v '_merged\.bam$' | wc -l)

  # Only merge if >1 BAM exists
  if [[ "$bamCount" -ge 2 ]]; then
    # Collect BAM and BAI files (exclude merged outputs)
    mapfile -t bamFiles < <(ls -1 "$folder"/*.bam 2>/dev/null | grep -v '_merged\.bam$')
    mapfile -t baiFiles < <(ls -1 "$folder"/*.bai 2>/dev/null | grep -v '_merged\.bai$')

    # Derive output name (first BAM's basename -> first two underscore fields)
    base0=$(basename "${bamFiles[0]}")
    prefix=$(echo "$base0" | cut -d_ -f1,2)
    outName="$folder/${prefix}_merged.bam"

    # Build Picard INPUT=... list dynamically (works for 2, 3, 4, ... BAMs)
    inputs=""
    for b in "${bamFiles[@]}"; do
      inputs="${inputs} INPUT=${b}"
    done

    # Build rm list: BAMs + BAIs (only those that exist)
    rm_list="${bamFiles[*]}"
    if [[ "${#baiFiles[@]}" -gt 0 ]]; then
      rm_list="${rm_list} ${baiFiles[*]}"
    fi

    echo "module load picard; \
java -Xmx3G -jar \$EBROOTPICARD/picard.jar MergeSamFiles \
ASSUME_SORTED=true CREATE_INDEX=true \
${inputs} \
OUTPUT=$outName MERGE_SEQUENCE_DICTIONARIES=false SORT_ORDER=coordinate \
USE_THREADING=false VALIDATION_STRINGENCY=STRICT; \
rm ${rm_list}"
  fi
done
