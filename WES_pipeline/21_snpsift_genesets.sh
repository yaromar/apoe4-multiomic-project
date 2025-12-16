#!/bin/bash
# 21_snpsift_genesets.sh
#
# Purpose:
#   Annotate the fully annotated VCF (*_snpeff_dbsnp_clinvar_gwas.vcf)
#   with gene set membership using MSigDB GMT.
#
# Notes from your runs:
#   - 6GB, 1 CPU was enough for ~1 hour for most files,
#     but more RAM and +1 CPU made it faster.
#   - we’ll use:
#       -c 2 --mem-per-cpu=20G -t 10:00:00
#
# Usage pattern (DSQ):
#   bash 21_snpsift_genesets.sh > joblist_21.txt
#   module load dSQ
#   dsq --job-file joblist_21.txt -c 2 --mem-per-cpu=20G -t 10:00:00 --max-jobs <SET_ME> --output DSQ_results_21
#   sbatch dsq-joblist_21.slurm

GMT=~/project/AD_WES_yaro/reference_genome/MSigDB/msigdb.v7.4.symbols.gmt

for folder in ~/scratch60/*; do
  file=($folder/*filtered_snpeff_dbsnp_clinvar_gwas.vcf)

  if [[ ! -e "${file[0]}" ]]; then
    continue
  fi

  input_vcf="${file[0]}"
  output_vcf="${input_vcf/.vcf/_genesets.vcf}"

  echo "module load snpEff; \
~/tools/snpEff/scripts/snpSift geneSets -v \"$GMT\" \"$input_vcf\" > \"$output_vcf\""
done
