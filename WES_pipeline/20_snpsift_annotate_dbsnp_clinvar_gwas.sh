#!/bin/bash
# 20_snpsift_annotate_dbsnp_clinvar_gwas.sh
#
# Purpose:
#   For each sample VCF (*filtered_snpeff.vcf), add external annotations in sequence:
#     1) dbSNP  -> *_dbsnp.vcf
#     2) ClinVar -> *_clinvar.vcf
#     3) GWAS Catalog -> *_gwas.vcf
#
# Output files are written in the SAME folder as the input VCF.
#
# Usage pattern (DSQ):
#   bash 20_snpsift_annotate_dbsnp_clinvar_gwas.sh > joblist_20.txt
#   module load dSQ
#   dsq --job-file joblist_20.txt -c 2 --mem-per-cpu=20G -t 10:00:00 --max-jobs <SET_ME> --output DSQ_results_20
#   sbatch dsq-joblist_20.slurm

for folder in ~/scratch60/*; do
  # Expect exactly one snpEff-annotated VCF per folder 
  file=($folder/*filtered_snpeff.vcf)

  # If glob doesn't match, skip folder quietly
  if [[ ! -e "${file[0]}" ]]; then
    continue
  fi

  input_vcf="${file[0]}"
  output_dbsnp="${input_vcf/.vcf/_dbsnp.vcf}"
  output_clinvar="${output_dbsnp/.vcf/_clinvar.vcf}"
  output_gwas="${output_clinvar/.vcf/_gwas.vcf}"

  # Emit one command line per sample for DSQ joblist
  echo "module load snpEff; \
~/tools/snpEff/scripts/snpSift annotate -v ~/project/AD_WES_yaro/reference_genome/dbSNP/00-All_chr.vcf.gz \"$input_vcf\" > \"$output_dbsnp\"; \
~/tools/snpEff/scripts/snpSift annotate -v ~/project/AD_WES_yaro/reference_genome/clinvar/clinvar.vcf \"$output_dbsnp\" > \"$output_clinvar\"; \
~/tools/snpEff/scripts/snpSift annotate -v ~/project/AD_WES_yaro/reference_genome/gwas/hg38_gwasCatalog.txt \"$output_clinvar\" > \"$output_gwas\""
done
