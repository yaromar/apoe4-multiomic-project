#!/bin/bash
# 12_collect_alignment_metrics.sh
#
# Purpose
# -------
# Collect alignment and quality metrics on final BAM files using Picard
# CollectMultipleMetrics.
#
# Metrics collected:
#   - QualityScoreDistribution
#   - CollectAlignmentSummaryMetrics
#
# Assumptions
# -----------
# - Each sample folder under ~/scratch60/* contains exactly ONE final BAM
#   (typically *_marked.bam from step 11).
# - BAM files are coordinate-sorted and indexed.
#
# Note on alignment metrics
# Alignment metrics are collected against the same decoy.fa reference used 
# for read alignment (GRCh38 + common Illumina contaminants). As a result, 
# mapped-read statistics include reads aligning to decoy contigs. 
# This ensures internal consistency across the pipeline but may yield slightly 
# higher mapping rates compared to genome-only references.
#
# DSQ usage
# ---------
# This script PRINTS one command per sample folder.
#
#   bash 12_collect_alignment_metrics.sh > joblist.txt
#   module load dSQ
#   dsq --job-file joblist.txt -c 1 --mem-per-cpu=6GB -t 5:00:00 --output DSQ_results
#
# Slurm parameters (already validated)
# -----------------------------------
#   #SBATCH --array 0-688%200
#   #SBATCH --output DSQ_results
#   #SBATCH --job-name dsq-joblist
#   #SBATCH -c 1 --mem-per-cpu=6GB -t 5:00:00

for folder in ~/scratch60/*; do
  # Prefer marked BAMs
  bam_candidates=("$folder"/*_marked.bam)

  # Skip folders without exactly one marked BAM
  if [[ "${#bam_candidates[@]}" -ne 1 ]]; then
    continue
  fi

  file="${bam_candidates[0]}"
  out_prefix="${file/.bam/}"

  echo "module load picard; module load R; \
java -Xmx3G -jar \$EBROOTPICARD/picard.jar CollectMultipleMetrics \
INPUT=$file \
OUTPUT=$out_prefix \
REFERENCE_SEQUENCE=/home/ym362/project/AD_WES_yaro/reference_genome/decoy.fa \
VALIDATION_STRINGENCY=LENIENT \
PROGRAM=QualityScoreDistribution \
PROGRAM=CollectAlignmentSummaryMetrics"
done
