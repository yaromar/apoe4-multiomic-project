#!/bin/bash
# 16_mutect2_calling.sh
#
# Purpose:
#   Run GATK Mutect2 in tumor-only mode (no matched normal), producing:
#     - unfiltered Mutect2 VCF
#     - F1R2 counts tarball (for orientation bias modeling later)
#
# DSQ usage:
#   1) Generate joblist:
#        bash 16_mutect2_calling.sh > joblist.txt
#   2) Run DSQ (parameters already validated):
#        module load dSQ
#        dsq --job-file joblist.txt \
#            -c 1 --mem-per-cpu=3GB -t 15:00:00 \
#            --max-jobs 200 \
#            --output DSQ_results
#
# Notes:
#   - Mutect2 requires a sample name for -tumor. We derive it from BAM filename.

REF=~/project/AD_WES_yaro/reference_genome/decoy.fa
TARGETS=~/project/AD_WES_yaro/ycga_idt_target_regions/hg38_idt_and_spikein_regions.bed
GNOMAD=~/project/AD_WES_yaro/reference_genome/gnomAD/somatic-hg38_af-only-gnomad.hg38.vcf.gz
PON=~/project/AD_WES_yaro/reference_genome/1000g_pon/somatic-hg38_1000g_pon.hg38.vcf.gz

for folder in ~/scratch60/*; do
  for bam in "$folder"/*.bam; do
    [[ -e "$bam" ]] || continue

    vcf="${bam%.bam}_mt2.vcf"
    f1r2="${bam%.bam}_f1r2.tar.gz"

    # Derive tumor sample name from BAM filename (basename without extension).
    tumor="$(basename "$bam" .bam)"

    echo "module load GATK; \
gatk Mutect2 \
  -R \"$REF\" \
  -L \"$TARGETS\" \
  -I \"$bam\" \
  -O \"$vcf\" \
  -tumor \"$tumor\" \
  --af-of-alleles-not-in-resource 0.00003125 \
  --germline-resource \"$GNOMAD\" \
  -pon \"$PON\" \
  --f1r2-tar-gz \"$f1r2\""
  done
done

# SBATCH --job-name dsq-joblist
# SBATCH -c 1 --mem-per-cpu=3GB -t 15:00:00
