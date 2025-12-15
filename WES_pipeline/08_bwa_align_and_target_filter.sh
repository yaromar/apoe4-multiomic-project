#!/bin/bash
# 08_bwa_align_and_target_filter.sh
#
# Purpose
# -------
# For each sample folder in ~/scratch60, align paired-end trimmed FASTQs to the
# GRCh38+decoy reference using bwa mem, add read group (RG) tags, and immediately
# filter alignments to on-target regions (IDT capture + spike-ins) using samtools view -L.
#
# DSQ usage
# ---------
# This script PRINTS one command per R1 file. Redirect to a joblist, then run via dSQ.
#   bash 08_bwa_align_and_target_filter.sh > joblist.txt
#   module load dSQ
#   dsq --job-file joblist.txt -c 8 --mem-per-cpu=7GB -t 10:00:00 --max-jobs 22 --output DSQ_results
#
# Slurm parameters (already validated)
# -----------------------------------
# The DSQ-generated .slurm should include something like:
#   #SBATCH --array 0-908%22
#   #SBATCH --output DSQ_results
#   #SBATCH --job-name dsq-joblist
#   #SBATCH -c 8 --mem-per-cpu=7GB -t 10:00:00
#
# Notes
# -----
# - Assumes TrimGalore output names contain "val_1"/"val_2".
# - RG fields:
#     PU derived from FASTQ header fields 3 and 4: "<flowcell>.<lane>"
#     ID derived from PU (first 5 chars of flowcell + lane)
#     SM derived from folder name (before _S###_) and cut to first two underscore-separated fields
#     LB derived from folder name (before _S###_) and last two underscore-separated fields
# - Output BAM is written as ${file/.fq.gz/.bam}. If your trimmed reads end in .fq.gz this works;
#   if they end in .fastq.gz, adjust the substitution below.
#
# Inputs
# ------
# - Reference:  ~/project/AD_WES_yaro/reference_genome/decoy.fa
# - Targets BED: ~/project/AD_WES_yaro/ycga_idt_target_regions/hg38_idt_and_spikein_regions.bed
#
# Outputs
# -------
# - One BAM per R1 file, filtered to targets.
# - Removes the paired trimmed FASTQs after successful BAM creation.

for folder in ~/scratch60/*; do
  # Only iterate over paired-end read 1 files
  for file in $(ls "$folder"/*R1* 2>/dev/null); do

    # Derive PU from the first FASTQ header line:
    # @<instrument>:<runid>:<flowcell>:<lane>:...
    # We use fields 3 and 4 => "<flowcell> <lane>" then make "<flowcell>.<lane>"
    fieldsARR=($(zcat "$file" | head -n1 | awk '{print $1}' | awk -F: '{print $3,$4}'))

    PU=${fieldsARR[0]}.${fieldsARR[1]}
    ID=${PU:0:5}.${fieldsARR[1]}

    # SM/LB derived from folder name
    # e.g. /.../scratch60/Sample_10260309_ST_BA22_I7_I5_S58_L004/...
    SM=$(echo "$file" | awk -F '/' '{print $6}' | awk -F '_S[0-9]+_' '{print $1}' | cut -f1,2 -d_)
    PL=ILLUMINA
    LB=$(echo "$file" | awk -F '/' '{print $6}' | awk -F '_S[0-9]+_' '{print $1}' | rev | cut -f1,2 -d_ | rev)

    RG="@RG\tID:$ID\tPU:$PU\tPL:$PL\tSM:$SM\tLB:$LB"

    # Expand R1 -> both mates; expand val_1 -> val_* (captures val_1 and val_2)
    temp=${file/R1/*}
    temp=${temp/val_1/val_*}

    # Align + tag RG + keep all alignments (-T 0) then filter to targets using -L BED
    # Write BAM and delete the input FASTQs after success.
    #
    echo "module load BWA; module load SAMtools; \
bwa mem -t 8 -R \"$RG\" -T 0 ~/project/AD_WES_yaro/reference_genome/decoy.fa $temp | \
samtools view -Shb -@ 5 -L ~/project/AD_WES_yaro/ycga_idt_target_regions/hg38_idt_and_spikein_regions.bed \
-o ${file/\.fq\.gz/.bam} - ; \
rm $temp"

  done
done
