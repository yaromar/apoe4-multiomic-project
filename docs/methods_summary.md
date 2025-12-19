# Methods summary

This document is a **methods snapshot** aligned with the manuscript and the analysis notebooks in `notebooks/`.

---

## Study design and cohort

### Source cohort and tissue
- **Cohort:** Religious Orders Study + Rush Memory and Aging Project (**ROSMAP**).
- **Tissue:** frozen **prefrontal cortex (Brodmann Area 10 / BA10)**.

### Neuropathology-based AD definition
- AD pathology status was defined using **dichotomized NIA–Reagan criteria**, integrating:
  - **Braak stage** (tau tangles)
  - **CERAD** (neuritic plaques)
- Participants were categorized as **AD** if they met **high or intermediate likelihood**, and **No AD** if **low likelihood or no AD pathology**.
- For some logistic regression analyses, **Braak** and **CERAD** were also analyzed as dichotomized outcomes (e.g., Braak 0–II vs ≥III; CERAD sparse/none vs moderate/frequent).

### Key covariates
Unless otherwise noted, statistical models adjusted for:
- age at death
- sex
- estimated neuronal proportion
- post-mortem interval (PMI)

---

## Molecular data and preprocessing

### Proteomics (DIA LC–MS/MS)

#### Sample preparation and LC–MS/MS acquisition (summary)
- Frozen tissue lysis in RIPA buffer with protease/phosphatase inhibitors; methanol:chloroform:water precipitation.
- Reduction (DTT) and alkylation (iodoacetamide), trypsin digestion, C18 desalting.
- DIA LC–MS/MS on **Waters ACQUITY UPLC M-Class** coupled to **Thermo Q-Exactive HFX**.
- DIA settings included 10 Th isolation windows; Orbitrap MS1 and MS2 acquired at typical high-resolution settings.

#### Identification and quantification
- DIA spectra searched/quantified using:
  - **Scaffold DIA v2.1.0**
  - mzML conversion via **ProteoWizard v3.0.11748**
  - peptide filtering via **Percolator v3.01** (peptide-level FDR threshold 0.01)
  - quantification via **EncyclopeDIA v0.9.2**
- Protein groups were constructed under parsimony principles; proteins filtered to **≥2 peptides/protein** with **protein FDR 1%**.

#### Missingness, normalization, and QC
- Protein IDs consolidated into protein groups; missing values imputed using a minimum-value strategy.
- Abundances log2-transformed.
- Technical variation adjusted using **RUV4** (R package `ruv`) with **negative control peptides**.
- Proteins with **>20% missingness** removed; retained proteins were **scaled** (mean 0, SD 1).
- Outliers were identified using PCA and removed prior to downstream analyses.

---

### DNA methylation (Illumina EPIC array)
- DNA methylation profiled using the **Illumina EPIC array**.
- Filtering steps:
  - removed CpGs on sex chromosomes (per EPIC manifest)
  - removed CpGs in the **lowest 10% variance**
- Neuronal proportion was estimated from DNA methylation using **CETS** (R).

---

## Baseline AD molecular associations

### Proteomics differential abundance
- Differential protein abundance (AD vs No AD) was tested using **limma** (`lmFit`, `eBayes`) on normalized protein values.
- Multiple testing control: **FDR** (primary threshold typically 0.05; secondary exploratory threshold 0.10 used for enrichment expansion when noted).
- Functional enrichment for protein sets was performed using **g:Profiler** with FDR correction.

### DNA methylation differential methylation
- Differential methylation (AD vs No AD) tested using **limma** with the same core covariates.
- Multiple testing control: **FDR** (primary 0.05; secondary 0.10 when used for broader context).

---

## APOE ε4 interaction framework

### APOE ε4 definition and stratified groups
- APOE genotype was dichotomized as **ε4 carrier vs non-carrier** (≥1 ε4 allele).
- For interpretation and visualization, participants were grouped into four genotype–phenotype strata:
  - ε4− / AD−
  - ε4− / AD+
  - ε4+ / AD−
  - ε4+ / AD+

### Interaction models
- APOE-dependent molecular effects were tested using interaction models of the form:

  **AD ~ APOE(ε4 carrier) × molecular_feature + age + sex + neuronal_proportion + PMI**

- Implemented using standard R regression (logistic regression for AD status).

### Residual-based subgroup comparisons
- For subgroup comparisons of molecular levels, covariate-adjusted values were derived as residuals from:

  **molecular_feature ~ age + sex + neuronal_proportion + PMI**

- Residual distributions were compared across groups using Wilcoxon tests where applicable.

---

## Mediation analysis (APOE ε4 → molecular feature → AD)

- Mediation analyses were performed using the R **`mediation`** package.
- For prioritized proteins / CpGs, mediation quantified how much of the **total APOE ε4 association with AD** could be statistically explained by each molecular mediator (e.g., ACME with confidence intervals and p-values).

---

## Neuropathology association analyses (CERAD and Braak)

- Associations between prioritized molecular features and neuropathology measures were tested using regression models (and/or dichotomized logistic models where specified), typically adjusting for the same covariates:
  - age, sex, neuronal proportion, PMI

---

## Protein network analysis (STRING + clustering)

- Protein–protein interaction (PPI) networks were constructed using **STRING** interactions for selected proteins.
- **ELAVL4** was included as a gene-level node to contextualize the CpG cg06329447 mapping.
- Network clustering used the **Markov Cluster Algorithm (MCL)**.
- Resulting clusters/modules were annotated using protein function and enrichment terms (when applicable).

---

## Polygenic risk scores (PRS)

### PRS definition
Two PRS were computed:
- **Non-APOE PRS:** genome-wide PRS excluding the APOE region (defined as GRCh37 chr19:40–50 Mb)
- **APOE-inclusive PRS:** genome-wide PRS including the APOE region

### Weights and methods
- PRS weights derived from large AD GWAS summary statistics (e.g., Wightman et al.).
- PRS computed using two Bayesian methods without manual parameter tuning:
  - **PRS-CS-auto**
  - **SDPR**

### PRS association analyses
- PRS associations tested with:
  - AD pathology (NIA–Reagan), Braak, CERAD
  - prioritized proteins and CpGs
- Multiple testing correction applied when evaluating PRS–molecule associations.

---

## Predictive modeling (genotype-stratified)

### Feature construction
- **Proteomics:** PCA performed on protein sets (e.g., network modules); the first PC per module used as a compact feature representation. Unclustered proteins were grouped and summarized similarly.
- **DNA methylation:** evaluated both single-site predictors (e.g., cg06329447) and module-level predictors (WGCNA eigengenes).
- **PRS:** both PRS types evaluated as predictors.

### Models and evaluation
- Logistic regression models of the form:

  **AD ~ predictors + age + sex**

- Performance estimated using **caret** with **5-fold cross-validation repeated 5 times**, using out-of-fold predictions to compute ROC curves and AUC.
- Models were evaluated **separately within APOE ε4 carriers and non-carriers** to quantify genotype-dependent predictive signal.

---

## WGCNA on DNA methylation

- WGCNA performed on preselected CpGs nominally associated with AD (**p < 0.05**) to focus on disease-relevant signal and reduce dimensionality.
- Network construction (blockwiseModules) used signed network
- Module eigengenes (PC1 per module) were:
  - associated with AD / CERAD / Braak (linear models with covariates)
  - tested for APOE ε4 moderation using:

    **AD ~ APOE(ε4 carrier) × module_eigengene + covariates**

---

## Whole-exome sequencing (WES) pipeline (TBD; implemented in `WES_pipeline/`)

WES processing is implemented as a **Slurm + dSQ** pipeline in `WES_pipeline/` using fixed, pre-validated parameters. The high-level steps include:

1. FASTQ download from Yale tape (`00_download_fastq.sh`)
2. Sample name/path normalization, lane merging (FASTQ/BAM)
3. Raw QC (FastQC/MultiQC), trimming (Trim Galore)
4. Alignment to hg38 decoy reference (BWA)
5. Sorting/merging, duplicate marking (Picard)
6. Base quality recalibration (GATK)
7. Contamination estimation (GetPileupSummaries/CalculateContamination)
8. Variant calling and filtering (Mutect2 + FilterMutectCalls; PASS-focused outputs)
9. Annotation (SnpEff/SnpSift; optional gene set annotations)
10. Compression/indexing and merging across samples (bgzip/tabix; multi-sample merges)

---

## Where to find implementation details
Notebook-specific model formulas and covariate sets are documented inline in `notebooks/`.
