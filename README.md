# APOEε4 Multi-Omic Project

**Lead author:** Yaroslav Markov  
**Primary affiliation (at time of work):** Yale University  

This repository contains code, minimal documentation, and analysis notebooks for a multi-omic study of APOEε4-associated proteomic and DNA methylation changes in Alzheimer disease using post-mortem human brain tissue.

---

## Project Motivation

Alzheimer disease (AD) is the leading cause of dementia worldwide and a major contributor to morbidity, mortality, and healthcare burden. The APOEε4 allele is the strongest common genetic risk factor for late-onset AD, but **how** APOEε4 shapes molecular changes in the brain—especially which changes are early, potentially compensatory versus late and overtly pathological—remains unclear.

This project addresses this using a uniquely **APOEε4-enriched subset of ROSMAP**, enabling us to dissect potential APOEε4-related risk and resilience mechanisms at the molecular level.

---

## Study Overview

We analyzed prefrontal cortex (BA10) tissue from a subset of the Religious Orders Study and Rush Memory and Aging Project (ROSMAP):

- **Proteomics (DIA LC–MS):** n = 302  
- **DNA methylation (EPIC array):** n = 310  
- **Genetic data / PRS (imputed WGS/WGS-derived genotypes):** n = 254  

Key cohort feature:  
- **69%** (n = 214) were APOEε4 carriers  
- **19.6% of carriers** (n = 42) had **no AD pathology** by NIA–Reagan criteria

We model:

- APOEε4 × AD interactions at the protein and CpG level  
- Mediation of APOEε4 → AD by selected proteins and CpG sites  
- Associations with neuropathology (CERAD, Braak) and polygenic risk scores (PRS)  
- Predictive performance of proteomic / epigenetic / PRS features stratified by APOEε4 status  

---

## Core Findings (Short Summary)

### 1. Biphasic model of APOEε4-driven pathology

In the absence of AD pathology, APOEε4 carriers show **lower levels of synaptic and metabolic proteins** (e.g., VAMP1, SYN3, CASKIN1, GLUD1, PI4KA), consistent with **early synaptic and metabolic vulnerability**.

In APOEε4 carriers with AD, we instead see **upregulation of inflammatory and proteostatic proteins**, including:

- GNAO1, AHNAK, FGG, HEBP1, APEX1, RAB4A, SLC12A5, LRP1, BAG6  
- Hypermethylation of cg06329447 in *ELAVL4*

Together, these support a **biphasic model**:  
early synaptic/metabolic disruptions → later inflammatory / proteostasis and epigenetic changes in APOEε4-associated AD.

### 2. APOEε4 × AD interaction and mediation

- Several proteins and CpG sites show **significant APOEε4 × AD interactions**, indicating genotype-specific molecular responses to pathology.  
- Mediation analyses identify:
  - Proteins such as **GRIPAP1, GSTK1, VAMP1, CASKIN1, DPP3, SYN3, FGG** mediating a substantial portion (~9–33%) of APOEε4’s effect on AD.
  - Hypermethylation at **cg06329447 in *ELAVL4*** mediating ~12% of the APOEε4 effect.

### 3. PRS: APOE-specific vs genome-wide risk

We computed two polygenic risk scores (PRS):

- **Non-APOE PRS:** genome-wide, excluding the APOE region  
- **APOE-inclusive PRS:** genome-wide, including APOE  

Key observations:

- The **non-APOE PRS** is strongly associated with AD neuropathology but **not** with the core APOEε4-related proteins or CpGs we identify.  
- The **APOE-inclusive PRS** is significantly associated with eight AD-related proteins in APOEε4 carriers, supporting that these markers are tied specifically to **APOEε4-related genetic risk**, not generalized AD polygenic risk outside APOE.

### 4. Predictive modeling stratified by APOEε4

In subjects with complete data for proteomics, DNAme, and PRS:

- In **non-carriers**, PRS performs best for AD classification (AUC ≈ 0.73).  
- In **carriers**, **proteomic and DNA methylation markers** outperform PRS (AUC up to ≈ 0.74).


## Repository Structure

```text
.
├── README.md                # This file
├── docs/                    # High-level project and methods documentation
│   ├── README.md            # Overview, methods, and FAQs
│   ├── data_overview.md     # Data types, cohorts, and key variables
│   ├── methods_summary.md   # Condensed methods (with refs to full paper)
│   └── authors_and_contributions.md
├── notebooks/
│   ├── 01_pheno_exploration_and_heatmaps.R
│   ├── 02_epigenetic_biomarker_associations.R
│   ├── 03_ewas_AD_status.R
│   ├── 04_wgcna_DNAm_modules.R
│   ├── 10_proteomics_import_and_cleaning.R
│   ├── 11_proteomics_normalization_and_integration.R
│   ├── 12_proteomics_ruv_normalization_and_differential.R
│   ├── 20_methylation_interactions_and_mediation.R
│   ├── 21_protein_interactions_and_mediation.R
│   ├── 30_prs_associations_and_modeling.R
│   └── 40_multiomic_integration_and_final_analyses.R
├── WES_pipeline/            
│   ├── 00_download_fastq.sh
│   ├── 01_fix_sample_names_and_paths.sh
│   ├── 02_merge_lanes_fastq.sh
│   ├── 03_quip_to_gz.slurm
│   ├── 04_fastqc_raw.slurm
│   ├── 05_multiqc_raw_QC.sh
│   ├── 06_cleanup_fastqc_reports.sh
│   ├── 07_trimgalore.slurm
│   ├── 08_compute_read_length_stats.slurm
│   ├── 09_build_decoy_reference.sh
│   ├── 10_bwa_align_and_target_filter.slurm
│   ├── 11_sort_bam_with_picard.slurm
│   ├── 12_merge_lanes_bam_with_picard.slurm
│   ├── 13_mark_duplicates.slurm
│   ├── 14_collect_alignment_metrics.slurm
│   ├── 15_base_recalibration.slurm
│   ├── 16_apply_bqsr.slurm
│   ├── 17_get_pileups_and_contamination.slurm
│   ├── 18_mutect2_calling.slurm
│   ├── 19_learn_read_orientation_and_sort_vcf.slurm
│   ├── 20_filter_mutect_calls.slurm
│   ├── 21_filter_variants_keep_pass.sh
│   ├── 22_snpeff_annotation.slurm
│   ├── 23_snpsift_annot_and_genesets.slurm
│   ├── 24_bgzip_and_tabix_vcfs.sh
│   ├── 25_merge_sample_vcfs.sh
│   └── readme_WES_pipeline.md            
├── config/
│   ├── sample_manifest.csv  # Sample IDs, APOE genotype, AD status, etc.
│   ├── proteomics_featureset.csv
│   ├── DNAm_featureset.csv
│   └── prs_metadata.csv
└── figures/
    ├── main_figures/        # Final figures used in the manuscript
    └── exploratory/         # Intermediate / exploratory plots
