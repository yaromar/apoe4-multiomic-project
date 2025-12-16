#!/bin/bash
# 13_base_recalibration.sh
#
# Purpose
# -------
# Generate BQSR recalibration tables for each (final) BAM using GATK BaseRecalibrator.
#
# Inputs
# ------
# - One BAM per sample folder under ~/scratch60/* (typically *_marked.bam)
# - Reference: decoy.fa (same reference used for alignment)
# - Known sites: dbSNP VCFs 
#
# Output
# ------
# - *_recal_data.table next to each BAM
#
# DSQ usage
# ---------
#   bash 13_base_recalibration.sh > joblist.txt
#   module load dSQ
#   dsq --job-file joblist.txt -c 1 --mem-per-cpu=1GB -t 4:00:00 --output DSQ_results
#
# Slurm parameters (already validated)
# -----------------------------------
#   #SBATCH --array 0-688%200
#   #SBATCH --output DSQ_results
#   #SBATCH --job-name dsq-joblist
#   #SBATCH -c 1 --mem-per-cpu=1GB -t 4:00:00

for folder in ~/scratch60/*; do
  # Prefer marked BAMs if present
  if ls "$folder"/*_marked.bam >/dev/null 2>&1; then
    bam_candidates=("$folder"/*_marked.bam)
  else
    bam_candidates=("$folder"/*.bam)
  fi

  # Require exactly one BAM to avoid accidental multi-BAM folders
  if [[ "${#bam_candidates[@]}" -ne 1 ]]; then
    continue
  fi

  file="${bam_candidates[0]}"

  echo "module load GATK; \
gatk BaseRecalibrator \
-R ~/project/AD_WES_yaro/reference_genome/decoy.fa \
--known-sites ~/project/AD_WES_yaro/reference_genome/dbSNP/00-All_chr.vcf.gz \
--known-sites ~/project/AD_WES_yaro/reference_genome/dbSNP/00-All_papu_chr.vcf.gz \
-I $file \
-O ${file/.bam/_recal_data.table}"
done
