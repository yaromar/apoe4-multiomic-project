#!/bin/bash
# 04_fastqc_raw.sh
#
# Purpose:
#   Generate a DSQ joblist to run FastQC on *raw* FASTQ(.gz) files under ~/scratch60/Sample*/
#
# Output:
#   FastQC outputs will be written to the working directory unless you set -o.
#   Recommended: write into a central folder like ~/scratch60/fastqc_raw/
#
# DSQ parameters (already validated by dry run):
#   #SBATCH --array 0-1817%200
#   #SBATCH --output DSQ_results
#   #SBATCH --job-name dsq-joblist_trim
#   #SBATCH -c 1
#   #SBATCH -t 3:00:00
#
# Usage:
#   bash 04_fastqc_raw.sh
#   module load dSQ
#   dsq --job-file joblist_fastqc_raw.txt -c 1 --mem-per-cpu=<MEM> -t 3:00:00 --max-jobs 200 --output DSQ_results
#   sbatch dsq-joblist_trim-*.slurm   # (filename depends on dsq output)
#
# Notes:
#   - FastQC is single-threaded by default here (-c 1).
#   - If you want a single output folder, set FASTQC_OUT below and keep -o.
#   - This script prints one FastQC command per input file into the joblist.

set -euo pipefail

JOBLIST="joblist_fastqc_raw.txt"
FASTQC_OUT="${FASTQC_OUT:-$HOME/scratch60/fastqc_raw}"

mkdir -p "${FASTQC_OUT}"
: > "${JOBLIST}"

for folder in "$HOME"/scratch60/*; do
  [ -d "$folder" ] || continue

  for file in "$folder"/*; do
    [ -f "$file" ] || continue

    # Optional: limit to FASTQ inputs (uncomment if needed)
    # [[ "$file" =~ \.fastq(\.gz)?$ ]] || continue

    echo "module load FastQC; fastqc -o ${FASTQC_OUT} ${file}" >> "${JOBLIST}"
  done
done

echo "Wrote joblist: ${JOBLIST}"
echo "FastQC output directory: ${FASTQC_OUT}"
