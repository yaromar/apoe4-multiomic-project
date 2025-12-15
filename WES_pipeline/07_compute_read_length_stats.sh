#!/bin/bash
# 07_compute_read_length_stats.sh
#
# Purpose:
#   Generate a DSQ joblist to compute mean read length for each trimmed FASTQ
#   (Trim Galore outputs *val_1.fq.gz / *val_2.fq.gz by default).
#
# What this script does:
#   - Prints one command per sample folder:
#       /home/ym362/project/AD_WES_yaro/MeanLength/get_mean.sh <folder>
#   - You redirect stdout to create joblist.txt
#
# How to run:
#   cd /gpfs/ysm/project/morgan_levine/ym362/AD_WES_yaro/MeanLength
#   bash 07_compute_read_length_stats.sh > joblist.txt
#
# Then submit with dSQ (parameters already validated):
#   module load dSQ
#   dsq --job-file joblist.txt --output DSQ_results
#
# Slurm / DSQ parameters (already validated):
#   #SBATCH --array 0-368%199
#   #SBATCH --output DSQ_results
#   #SBATCH --job-name dsq-joblist
#   #SBATCH -c 1
#   #SBATCH --mem-per-cpu=10G
#   #SBATCH -t 10:00

set -euo pipefail

for folder in ~/scratch60/*; do
  [ -d "$folder" ] || continue
  echo "/home/ym362/project/AD_WES_yaro/MeanLength/get_mean.sh $folder"
done
