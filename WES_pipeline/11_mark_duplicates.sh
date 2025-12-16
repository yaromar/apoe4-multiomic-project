#!/bin/bash
# 11_mark_duplicates.sh
#
# Purpose
# -------
# Run Picard MarkDuplicates on each (merged/sorted) BAM per sample folder, create index,
# write metrics, then remove the input BAM + its .bai to save space.
#
# Assumptions
# -----------
# - Each folder under ~/scratch60/* contains exactly one coordinate-sorted BAM to process
#   (typically *_merged.bam from step 10, or a single sorted BAM if only one lane existed).
# - The BAM has a matching .bai in the same folder.
#
# DSQ usage
# ---------
# This script PRINTS one command per folder.
#   bash 11_mark_duplicates.sh > joblist.txt
#   module load dSQ
#   dsq --job-file joblist.txt -c 2 --mem-per-cpu=15GB -t 5:00:00 --output DSQ_results
#
# Slurm parameters (already validated)
# -----------------------------------
# DSQ-generated .slurm should include something like:
#   #SBATCH --output DSQ_results
#   #SBATCH --job-name dsq-joblist
#   #SBATCH -c 2 --mem-per-cpu=15GB -t 5:00:00

for folder in ~/scratch60/*; do
  # Prefer merged BAM if present; otherwise take the only BAM in the folder.
  # Exclude already-marked BAMs to avoid reprocessing.
  bam_candidates=()
  if ls "$folder"/*_merged.bam >/dev/null 2>&1; then
    bam_candidates=("$folder"/*_merged.bam)
  else
    bam_candidates=("$folder"/*.bam)
  fi

  # Filter out *_marked.bam if glob picked it up
  bam_candidates=($(printf "%s\n" "${bam_candidates[@]}" 2>/dev/null | grep -v '_marked\.bam$'))

  # If there isn't exactly one BAM to run on, skip (safest behavior for DSQ joblists)
  if [[ "${#bam_candidates[@]}" -ne 1 ]]; then
    continue
  fi

  file="${bam_candidates[0]}"
  out="${file/.bam/_marked.bam}"
  metrics="${file/.bam/_marked_metrics.txt}"
  bai_in="${file/.bam/.bai}"

  echo "module load picard; \
java -Xmx12G -jar \$EBROOTPICARD/picard.jar MarkDuplicates \
MAX_RECORDS_IN_RAM=4000000 CREATE_INDEX=true \
INPUT=$file OUTPUT=$out \
ASSUME_SORT_ORDER=coordinate VALIDATION_STRINGENCY=STRICT \
METRICS_FILE=$metrics; \
rm $file ${bai_in}"
done
