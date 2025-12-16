#!/bin/bash
# 19_snpeff_annotation.sh
#
# Purpose:
#   Run snpEff annotation on each sample’s filtered VCF, producing:
#     1) an annotated VCF
#     2) an HTML snpEff summary report
#     3) a CSV stats file (via -csvStats)
#
# Usage:
#   1) Generate a DSQ joblist:
#        bash 19_snpeff_annotation.sh > joblist_19_snpeff.txt
#   2) Submit with dSQ (parameters already validated):
#        module load dSQ
#        dsq --job-file joblist_19_snpeff.txt -c 1 --mem-per-cpu=8G -t 02:00:00 --max-jobs 160 --output DSQ_results_filtered
#   3) sbatch dsq-joblist_filtered-*.slurm  (whatever dsq outputs)

SNP_EFF_JAR="$HOME/tools/snpEff/snpEff.jar"
SNP_EFF_DB="GRCh38.99"
INTERVALS="$HOME/project/AD_WES_yaro/ycga_idt_target_regions/hg38_idt_and_spikein_regions.bed"

for folder in ~/scratch60/*; do
  # Expecting uncompressed VCF from the Step 18 filtering, named *_filtered.vcf
  file=($folder/*_filtered.vcf)

  # If no match, skip
  [[ -e "${file[0]}" ]] || continue

  for vcf in "${file[@]}"; do
    # Output naming
    html="${vcf/_filtered.vcf/_snpeff.html}"
    csv="${vcf/_filtered.vcf/_snpeff.csv}"
    outvcf="${vcf/.vcf/_snpeff.vcf}"

    echo "module load Java; java -Xmx6g -jar ${SNP_EFF_JAR} -o gatk ${SNP_EFF_DB} \
-csvStats ${csv} -s ${html} -interval ${INTERVALS} ${vcf} > ${outvcf}"
  done
done


# SBATCH --job-name dsq-joblist_filtered
# SBATCH -c 1 --mem-per-cpu=8G -t 02:00:00
