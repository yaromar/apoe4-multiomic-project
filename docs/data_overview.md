# Data overview

This document summarizes **what data exist**, **how the cohort is structured**, and the **key variables** used throughout the analyses and notebooks in `notebooks/`.

> **Important note:** Raw data are not distributed in this GitHub repository. Internal lab data locations are listed in the root `README.md` (Higgins-Chen lab).

---

## Cohort and tissue

### Study population
Data were derived from participants in **ROSMAP** (Religious Orders Study + Rush Memory and Aging Project). Participants enrolled without known dementia, underwent longitudinal clinical evaluation, and agreed to brain donation at death. Phenotype and neuropathology data were provided by the Rush Alzheimer’s Disease Center.

### Tissue source
- **Brain region:** prefrontal cortex, **Brodmann Area 10 (BA10)**  
- **Material:** frozen post-mortem tissue

---

## Core cohort labels

### APOE genotype
Participants are categorized by APOE genotype and typically represented as:
- **APOE ε4 carrier (APOEε4+)**: at least one ε4 allele (e.g., ε3/ε4 or ε4/ε4)
- **APOE ε4 non-carrier (APOEε4−)**: no ε4 allele (e.g., ε3/ε3)

> The cohort is intentionally **enriched for APOEε4 carriers**:
- **69.0%** of the total sample are ε4 carriers (**n = 214**)
- **19.6% of ε4 carriers** (**n = 42**) have **no AD pathology** (NIA–Reagan)

### AD pathology status (primary outcome)
**AD pathology** is defined using a dichotomized **NIA–Reagan** criterion, based solely on post-mortem neuropathology (blinded to clinical data):
- **AD**: intermediate or high likelihood
- **No AD**: low likelihood or no AD pathology

This criterion integrates:
- **CERAD** neuritic plaque score
- **Braak** neurofibrillary tangle stage

### Four genotype–phenotype subgroups
For interpretation and subgroup comparisons, analyses often use four biologically meaningful groups:
1. **APOEε4− / AD−**
2. **APOEε4− / AD+**
3. **APOEε4+ / AD−**
4. **APOEε4+ / AD+**

These subgroups support the project’s framing around **risk** and **resilience** among ε4 carriers.

---

## Neuropathology variables (key outcomes)

### NIA–Reagan (dichotomized)
- **Primary binary AD label** used across most analyses

### Braak stage (tau)
Often used in two forms:
- **ordinal/continuous** staging (0–VI)
- **dichotomized** for logistic regression:
  - **0–II** vs **≥III**
  - (0–II: tangles largely confined to entorhinal region; ≥III: limbic/neocortical spread)

### CERAD score (amyloid plaques)
Often used in two forms:
- **ordinal score** (none/sparse/moderate/frequent)
- **dichotomized** for logistic regression:
  - **none/sparse** vs **moderate/frequent**
  - (moderate/frequent corresponds to probable/definite AD plaque burden)

---

## Covariates used across modalities

These variables are included routinely to reduce confounding:
- **age at death**
- **sex**
- **post-mortem interval (PMI)**
- **estimated neuronal proportion** (derived from DNA methylation using CETS)

Additional covariates may appear in specific notebooks (e.g., batch indicators, technical QC metrics), but the list above is the “core set.”

---

## Molecular modalities and sample sizes

### Summary table (by modality)

| Modality | Platform / type | N (available) | Notes |
|---|---|---:|---|
| Proteomics | DIA LC–MS | 302 | Quantified with Scaffold DIA; protein grouping + filtering applied |
| DNA methylation | Illumina EPIC array | 310 | Used for EWAS, interaction models, neuron proportion estimation, WGCNA |
| Genetics / PRS | WGS-derived genotypes → PRS | 254 | PRS computed using published GWAS weights; APOE-inclusive and non-APOE PRS |
| Neuropathology | CERAD, Braak, NIA–Reagan | 310 | Used as outcomes and severity correlates |

> Overlap across modalities is substantial but not complete (e.g., PRS available in a subset). Some downstream analyses restrict to participants with complete multi-omic data.

---

## Proteomics data (DIA LC–MS)

### What exists
- Raw DIA mass spec files 

### Quantification and filtering (high level)
- DIA quantification identified **4,901** proteins (including isoforms/PTM variants).
- After grouping and filtering proteins with **>20% missingness**, **1,291** proteins were retained for downstream analyses.


---

## DNA methylation data (EPIC array)

### What exists
- EPIC methylation matrices (beta values; internal)
- Derived estimates:
  - neuronal proportion (CETS)
  - CpG module eigengenes (WGCNA)

---

## Genetics and PRS

### What exists
- Individual-level PRS score tables 
- Two PRS definitions:
  - **Non-APOE PRS**: genome-wide PRS excluding the APOE region
  - **APOE-inclusive PRS**: genome-wide PRS including the APOE region

---


## Cohort structure and overlap

- Each subject belongs to exactly one APOE group (carrier/non-carrier) and one NIA–Reagan group (AD/No AD).
- Molecular data are distributed across these groups with broadly comparable sampling proportions.
- Some analyses restrict to:
  - modality-specific availability (e.g., PRS subset)
  - complete-case multi-omic overlap (e.g., proteomics + DNAme + PRS)
