#!/bin/bash
# 15_get_pileups_and_contamination.sh
#
# Purpose:
#   Generate pileup summary tables for contamination estimation (GATK GetPileupSummaries).
#
# DSQ usage:
#   1) Generate joblist:
#        bash 15_get_pileups_and_contamination.sh > joblist.txt
#   2) Run DSQ (parameters already validated):
#        module load dSQ
#        dsq --job-file joblist.txt \
#            -c 1 --mem-per-cpu=8GB -t 4:00:00 \
#            --max-jobs 160 \
#            --output DSQ_results
#
# Expected output per BAM:
#   <sample>_pileups.table

REF=~/project/AD_WES_yaro/reference_genome/decoy.fa
GNOMAD=~/project/AD_WES_yaro/reference_genome/gnomAD/somatic-hg38_af-only-gnomad.hg38.vcf.gz
TARGETS=~/project/AD_WES_yaro/ycga_idt_target_regions/hg38_idt_and_spikein_regions.bed

for folder in ~/scratch60/*; do
  for bam in "$folder"/*.bam; do
    [[ -e "$bam" ]] || continue
    out="${bam%.bam}_pileups.table"

    echo "module load GATK; \
gatk GetPileupSummaries \
  -R \"$REF\" \
  -V \"$GNOMAD\" \
  -L \"$TARGETS\" \
  -I \"$bam\" \
  -O \"$out\""
  done
done
# SBATCH --job-name dsq-joblist
# SBATCH -c 1 --mem-per-cpu=8GB -t 4:00:00
