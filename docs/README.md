# APOE ε4–focused multi-omic Alzheimer’s disease project (ROSMAP PFC)

This repository documents a multi-omic study of **APOE ε4–specific molecular mechanisms in Alzheimer’s disease (AD)** using post-mortem **prefrontal cortex (PFC)** brain tissue from **ROSMAP**. The goal is to provide a clear, reusable, notebook-style record of motivation, study design, analyses, and results **without distributing sensitive raw data**.


---

## Study motivation

APOE ε4 is the strongest common genetic risk factor for late-onset AD and is implicated in a broad set of biological processes (lipid/mitochondrial function, oxidative stress, inflammation, synaptic and cerebrovascular biology). 

The study uses a uniquely structured cohort enriched for APOE ε4 carriers, including a substantial subset of ε4 carriers with **no AD pathology**, enabling analyses of both **risk** and **resilience** mechanisms.

(For full context and results, see the associated manuscript/preprint.)

---

## Cohort and phenotypes

### Cohort source
- **ROSMAP** (Religious Orders Study + Memory and Aging Project)
- **Tissue:** post-mortem **prefrontal cortex (PFC)**

### AD case/control definition
AD status is based on a **neuropathology-only** classification using a dichotomized **NIA–Reagan** definition (CERAD plaque burden + Braak staging), evaluated blinded to clinical data.


---

## Molecular data types

This repo contains analysis notebooks for:

- **Proteomics (brain tissue):**
  - LC–MS with Data Independent Acquisition (DIA)
  - Protein-level QC and filtering prior to downstream analyses

- **DNA methylation:**
  - Illumina EPIC array
  - EWAS, interaction models, and co-methylation network analyses

- **Polygenic risk scores (PRS):**
  - AD PRS derived from genotype data using published GWAS weights

- **Neuropathology and cognitive phenotypes:**
  - CERAD, Braak, NIA-Reagan 


---

## High-level analysis plan

The analysis is organized around four core questions:

1. **Which molecular features differ in AD within APOE ε4 carriers vs non-carriers?**  
   (Genotype-stratified differential analyses + APOE×AD interaction models)

2. **Which features may mediate APOE ε4’s effect on AD pathology?**  
   (Mediation models on prioritized proteins and CpG sites)

3. **How do candidate features relate to neuropathology severity?**  
   (Associations with CERAD, Braak, NIA-Reagan)

4. **Do different modalities predict AD differently by APOE genotype?**  
   (Predictive models stratified by APOE ε4 status across proteomics, methylation, and PRS)

---

## Summary of key findings (paper-aligned)

### 1) Biphasic APOE ε4 molecular model
Results support a **dual-phase** framework:
- **Early phase (preclinical vulnerability in ε4 carriers):** reduced proteins linked to synaptic integrity and metabolic homeostasis.
- **Later phase (ε4 carriers with AD):** increased inflammatory / proteostatic / stress-response proteins and distinct epigenetic signals.

### 2) Candidate protein mediators of APOE ε4 → AD
Mediation analyses prioritize a small set of proteins as statistically supported intermediaries (e.g., GRIPAP1, GSTK1, VAMP1, CASKIN1, DPP3, SYN3, FGG), suggesting specific mechanisms that may link APOE genotype to pathological outcomes.

### 3) Epigenetic mediator and co-methylation structure
A CpG site associated with ELAVL4 (cg06329447) shows evidence of mediating a portion of APOE ε4’s effect and relates to neuropathology burden. Co-methylation network analyses (WGCNA) are used to capture broader epigenetic dysregulation beyond single-site signals.

### 4) Genotype-stratified predictive performance differs by modality
Predictive modeling suggests **PRS performs better in non-carriers**, while **non-genetic markers (especially DNA methylation summaries) are comparatively more informative in ε4 carriers**, motivating genotype-aware, multi-modal prediction strategies.

---

## How to navigate this repository

### Start here
- `README.md` (repo entry point)
- `docs/README.md` (this file: study overview + FAQs)
- `docs/data_overview.md` (what data exist, key variables, cohort structure)
- `docs/methods_summary.md` (methods snapshot)

### Analysis scripts
Notebooks in `notebooks/` are numbered to reflect the typical analysis flow:
- phenotype exploration → DNAm/proteomics preprocessing → interaction + mediation → PRS → multi-omic integration

### WES pipeline
The `WES_pipeline/` directory contains a Slurm + dSQ oriented pipeline with fixed parameters and joblist generation scripts.

---

## FAQs

### Does this repo include raw data?
No. This public release focuses on documentation and analysis structure. Raw omics files (e.g., IDATs, raw MS files, individual-level genotypes/VCFs) are not distributed here.

### Can I reproduce every figure from the paper using only this repo?
Not from scratch without controlled-access data. The repo is designed to make the *analysis logic* transparent and reusable.

### How can collaborators run the notebooks internally?
Recommended approach:
1. Replace file paths in the config section(s) of each notebook.
2. Provide a harmonized metadata table matching the variables described in `docs/data_overview.md`.

### What is meant by “APOE ε4 resilience” here?
Operationally: **APOE ε4 carriers with no AD pathology** under the neuropathology-based AD definition (NIA–Reagan dichotomization). Analyses use this group to probe molecular features associated with protection/resilience despite genetic risk.

### Where do I find the exact statistical models?
- High-level: `docs/methods_summary.md`
- Notebook-specific details: see the “Models” or “Methods” sections inside each notebook (interaction terms, covariates, multiple testing correction, and mediation setup).
