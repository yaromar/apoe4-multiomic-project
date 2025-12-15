#!/bin/bash
# 09_sort_bam_with_picard.sh
#
# Purpose
# -------
# Sort each BAM by coordinate using Picard SortSam and create a BAM index.
# After a successful sort, delete the unsorted BAM.
#
# DSQ usage
# ---------
# This script PRINTS one command per BAM. Redirect to a joblist, then run via dSQ.
#   bash 09_sort_bam_with_picard.sh > joblist.txt
#   module load dSQ
#   dsq --job-file joblist.txt -c 1 --mem-per-cpu=6GB -t 10:00:00 --max-jobs 200 --output DSQ_results
#
# Slurm parameters (already validated)
# -----------------------------------
# The DSQ-generated .slurm should include something like:
#   #SBATCH --array 0-908%200
#   #SBATCH --output DSQ_results
#   #SBATCH --job-name dsq-joblist
#   #SBATCH -c 1 --mem-per-cpu=6GB -t 10:00:00
#
# Inputs
# ------
# - Unsorted BAMs in ~/scratch60/*/*.bam
#
# Outputs
# -------
# - Sorted BAM: *_sorted.bam
# - BAM index: *_sorted.bai (Picard uses .bai)
#
# Notes
# -----
# - This script only targets BAMs that DO NOT already end with "_sorted.bam"
# - Output naming preserves the original prefix and appends "_sorted.bam"
#   (more robust than substituting a specific substring like "_R1_001_val_1.bam")

for folder in ~/scratch60/*; do
  for file in "$folder"/*.bam; do
    # Skip if no BAMs in folder
    [ -e "$file" ] || continue

    # Skip already-sorted BAMs
    if [[ "$file" == *"_sorted.bam" ]]; then
      continue
    fi

    # Define output name robustly
    out="${file%.bam}_sorted.bam"

    # Emit Picard SortSam command
    echo "module load picard; \
java -Xmx4G -jar \$EBROOTPICARD/picard.jar SortSam \
MAX_RECORDS_IN_RAM=1000000 CREATE_INDEX=true \
INPUT=$file OUTPUT=$out SORT_ORDER=coordinate VALIDATION_STRINGENCY=STRICT; \
rm $file"
  done
done
