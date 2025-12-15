#!/bin/bash
# 06_trimgalore.sh
#
# Purpose:
#   Adapter and quality trimming of paired-end FASTQ files using Trim Galore.
#   This script GENERATES a DSQ joblist (one command per sample),
#   which is then executed as a Slurm array job.
#
# Notes:
#   - Assumes merged lanes (single R1/R2 per sample)
#   - Assumes files are named with R1 / R2 in the filename
#   - Original FASTQs are removed after successful trimming to save space
#
# Execution model:
#   1) Run this script to generate joblist.txt
#   2) Submit via dSQ using FIXED, PRE-VALIDATED resources (see Slurm header below)
#
# Example:
#   bash 06_trimgalore.sh > joblist.txt
#   module load dSQ
#   dsq --job-file joblist.txt --output DSQ_results
#
# Resources (already validated on largest samples):
#   CPUs: 15
#   Walltime: 5 hours
#   Max concurrent jobs: 13
#
# Corresponding Slurm parameters (auto-injected by dSQ):
#   #SBATCH --array 0-908%13
#   #SBATCH --output DSQ_results
#   #SBATCH --job-name dsq-joblist
#   #SBATCH -c 15
#   #SBATCH -t 5:00:00

set -euo pipefail

for folder in ~/scratch60/*; do
  # Loop over R1 files only (paired-end assumed)
  for r1 in "${folder}"/*R1*.fastq.gz; do
    # Skip if no R1 files exist
    [ -e "$r1" ] || continue

    # Infer R2 filename
    r2="${r1/R1/R2}"

    # Sanity check
    if [ ! -f "$r2" ]; then
      echo "# WARNING: Missing R2 for ${r1}" >&2
      continue
    fi

    # Emit Trim Galore command
    echo "module load Trim_Galore; \
trim_galore --paired --cores 8 -o ${folder} ${r1} ${r2} && \
rm ${r1} ${r2}"
  done
done
