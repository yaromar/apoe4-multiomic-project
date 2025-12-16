#!/bin/bash
# 17_learn_read_orientation_and_sort_vcf.sh
#
# Purpose:
#   Per-sample post-Mutect2 steps:
#     1) LearnReadOrientationModel from F1R2 counts
#     2) SortVcf (and index)
#     3) FilterMutectCalls using:
#         - orientation-bias priors
#         - contamination table
#         - (optional) segmentation table
#
# DSQ usage:
#   1) Generate joblist:
#        bash 17_learn_read_orientation_and_sort_vcf.sh > joblist.txt
#   2) Run DSQ (parameters already validated):
#        module load dSQ
#        dsq --job-file joblist.txt \
#            -c 1 -t 5:00:00 \
#            --max-jobs 200 \
#            --output DSQ_results
#
# Notes:
#   - This script assumes these exist per BAM (from prior steps):
#       *_f1r2.tar.gz              (Mutect2)
#       *_mt2.vcf                  (Mutect2)
#       *_contamination.table      (CalculateContamination step)
#       *_segments.table           (optional; from tumor segmentation workflow)

REF=/home/ym362/project/AD_WES_yaro/reference_genome/decoy.fa
DICT=/home/ym362/project/AD_WES_yaro/reference_genome/decoy.dict
TARGETS=~/project/AD_WES_yaro/ycga_idt_target_regions/hg38_idt_and_spikein_regions.bed

for folder in ~/scratch60/*; do
  for bam in "$folder"/*.bam; do
    [[ -e "$bam" ]] || continue

    f1r2="${bam%.bam}_f1r2.tar.gz"
    rom="${bam%.bam}_read-orientation-model.tar.gz"

    mt2="${bam%.bam}_mt2.vcf"
    mt2_sorted="${bam%.bam}_mt2_sorted.vcf"

    contam="${bam%.bam}_contamination.table"
    seg="${bam%.bam}_segments.table"

    filtered="${bam%.bam}_mt2_contFiltered.vcf.gz"

    # The original stats filename was "${file/.bam/_mt2.vcf.stat}" which expands to *_mt2.vcf.stat
    stats="${mt2}.stat"

    echo "module load GATK; module load picard; \
gatk LearnReadOrientationModel -I \"$f1r2\" -O \"$rom\"; \
java -jar \$EBROOTPICARD/picard.jar SortVcf \
  SEQUENCE_DICTIONARY=\"$DICT\" \
  I=\"$mt2\" \
  OUTPUT=\"$mt2_sorted\" \
  CREATE_INDEX=true; \
gatk FilterMutectCalls \
  -R \"$REF\" \
  -L \"$TARGETS\" \
  -V \"$mt2_sorted\" \
  -O \"$filtered\" \
  --ob-priors \"$rom\" \
  --contamination-table \"$contam\" \
  --tumor-segmentation \"$seg\" \
  --stats \"$stats\""
  done
done
