# =============================================================================
# 02_epigenetic_biomarker_associations.R
# Yaro
# Last update: 12/6/23
#
# PURPOSE
# -------
# Test associations between a set of epigenetic / mutational / SEM biomarkers
# (across brain regions) and AD-related phenotypes/pathology, while adjusting
# for covariates.
#
# Output:
#   1) A long-format results table (results_df) with one row per:
#        biomarker (Var1) x phenotype (Var2)
#   2) A coefficient heatmap of biomarker ~ phenotype relationships.
#
# IMPORTANT IMPLEMENTATION NOTES
# ------------------------------
# - For visualization convenience, non-significant associations are assigned
#   Coefficient = 0 and t_Value = 0 (P_Value still recorded).
# =============================================================================

# -----------------------------------------------------------------------------
# Libraries
# -----------------------------------------------------------------------------
library(lme4)           
library(dplyr)
library(ordinal)        
library(ggplot2)        
library(ComplexHeatmap) # heatmap
library(reshape2)       # acast
library(circlize)       # colorRamp2

# -----------------------------------------------------------------------------
# Inputs assumed to exist in the environment
# -----------------------------------------------------------------------------
# big_pheno_scaled must already exist (from 01_pheno_exploration_and_heatmaps.R)

# -----------------------------------------------------------------------------
# Biomarkers (predictors)
# -----------------------------------------------------------------------------
clockColumns <- c(
  "mutational_burden_st", "mutational_burden_pfc", "mutational_burden_cbm",
  "sum_hypo_pfc", "sum_hypo_st", "sum_hypo_cbm",
  "sum_hyper_pfc", "sum_hyper_cbm", "sum_hyper_st",
  "meth_hypo_pfc", "meth_hypo_cbm", "meth_hypo_st",
  "mid_hyper_pfc", "mid_hyper_cbm", "mid_hyper_st",
  "mid_hypo_pfc", "mid_hypo_cbm", "mid_hypo_st",
  "unmeth_hyper_pfc", "unmeth_hyper_cbm", "unmeth_hyper_st",
  "PCPhenoAge_pfc", "PCPhenoAge_cbm", "PCPhenoAge_st",
  "Horvath1_pfc", "Horvath1_cbm", "Horvath1_st",
  "dunedinPACE_pfc", "dunedinPACE_cbm", "dunedinPACE_st",
  "DNAmClockCortical_cbm", "DNAmClockCortical_pfc", "DNAmClockCortical_st"
)

# -----------------------------------------------------------------------------
# Covariates
# -----------------------------------------------------------------------------
covariates <- c("Age", "Female", "pmi", "neurons_st")

# -----------------------------------------------------------------------------
# Phenotypes / outcomes (choose ONE MODE explicitly)
# -----------------------------------------------------------------------------
MODE <- "binary"   # options: "numeric", "ordinal", "binary"

cognitive_cols_numeric <- c(
  "gpath",
  "plaq_d", "plaq_n", "nft", "count", "amyloid",
  "cognep_random_slope", "cogng_random_slope", "cognpo_random_slope",
  "cognps_random_slope", "cognse_random_slope", "cognwo_random_slope",
  "tangles", "neurons_pfc", "neurons_cbm", "neurons_st",
  "ca1hipdp", "entodp", "infpardp", "midfrontdp", "midtempdp",
  "ca1hipnp", "entonp", "infparnp", "midfrontnp", "midtempnp",
  "ca1hipnft", "entonft", "infparnft", "midfrontnft", "midtempnft",
  "amyloid_ag", "amyloid_calc", "amyloid_cg", "amyloid_ec", "amyloid_hip",
  "amyloid_it", "amyloid_mf", "amyloid_mt", "amyloid_sf",
  "tangles_ag", "tangles_calc", "tangles_cg", "tangles_ec",
  "tangles_hip", "tangles_it", "tangles_mf", "tangles_mt",
  "age_first_ad_dx", "age_first_dem_dx"
)

cognitive_cols_ordinal <- c(
  "tdp_stage4", "tomm40_hap", "AD_Clinical", "cogdx", "braaksc"
)

cognitive_cols_binary <- c(
  "ci_num2_gct", "ci_num2_mct",
  "AD_Reagan", "ceradsc", "tdp_stage4_recode", "apoe_genotype"
)

cognitive_cols <- switch(
  MODE,
  numeric = cognitive_cols_numeric,
  ordinal = cognitive_cols_ordinal,
  binary  = cognitive_cols_binary
)

# -----------------------------------------------------------------------------
# Basic checks (fail early if something is missing)
# -----------------------------------------------------------------------------
stopifnot(exists("big_pheno_scaled"))
missing_predictors <- setdiff(clockColumns, colnames(big_pheno_scaled))
missing_outcomes    <- setdiff(cognitive_cols, colnames(big_pheno_scaled))
missing_covs        <- setdiff(covariates, colnames(big_pheno_scaled))

if (length(missing_predictors) > 0) stop("Missing biomarker columns: ", paste(missing_predictors, collapse = ", "))
if (length(missing_outcomes) > 0)   stop("Missing outcome columns: ", paste(missing_outcomes, collapse = ", "))
if (length(missing_covs) > 0)       stop("Missing covariate columns: ", paste(missing_covs, collapse = ", "))

# -----------------------------------------------------------------------------
# Results container
# -----------------------------------------------------------------------------
results_df <- data.frame(
  Var1 = character(),
  Var2 = character(),
  Coefficient = numeric(),
  P_Value = numeric(),
  t_Value = numeric(),
  stringsAsFactors = FALSE
)

# -----------------------------------------------------------------------------
# Modeling loop
# -----------------------------------------------------------------------------
# For each biomarker (predictor), fit one model per outcome.
# Model form:
#   outcome ~ biomarker + covariates
#
# Model type depends on MODE:
#   numeric:  lm()
#   binary:   glm(family = binomial)   (default here matches MODE)
#   ordinal:  clm()                    (ordinal package)
# -----------------------------------------------------------------------------

for (col in clockColumns) {
  for (compare_col in cognitive_cols) {
    
    formula_str <- paste(compare_col, "~", col, "+", paste(covariates, collapse = " + "))
    print(formula_str)
    
    fit <- switch(
      MODE,
      numeric = lm(as.formula(formula_str), data = big_pheno_scaled),
      binary  = glm(as.formula(formula_str), data = big_pheno_scaled, family = "binomial"),
      ordinal = clm(as.formula(formula_str), data = big_pheno_scaled)
    )
    
    sm <- summary(fit)
    
    coef_tab <- sm$coefficients
    
    # We expect a row named exactly like the biomarker column (col).
    # If not found, stop immediately so we don't silently index the wrong thing.
    if (!(col %in% rownames(coef_tab))) {
      stop(
        "Could not find biomarker term in coefficient table.\n",
        "Outcome: ", compare_col, "\n",
        "Biomarker: ", col, "\n",
        "Available terms: ", paste(rownames(coef_tab), collapse = ", ")
      )
    }
    
    # P-value column name differs by model
    p_col <- intersect(colnames(coef_tab), c("Pr(>|t|)", "Pr(>|z|)"))
    if (length(p_col) != 1) {
      stop("Could not uniquely identify p-value column in coefficients table for model: ", MODE)
    }
    p_col <- p_col[[1]]
    
    # Extract estimate / test-stat / pvalue for the biomarker term
    est <- coef_tab[col, 1]
    stat <- coef_tab[col, 3]
    pval <- coef_tab[col, p_col]
    
    # Store results:
    # - significant: keep estimate/stat
    # - non-significant: store zeros for coefficient/stat (but preserve p-value)
    if (!is.na(pval) && pval < 0.05) {
      results_df <- rbind(
        results_df,
        data.frame(Var1 = col, Var2 = compare_col, Coefficient = est, t_Value = stat, P_Value = pval)
      )
    } else {
      results_df <- rbind(
        results_df,
        data.frame(Var1 = col, Var2 = compare_col, Coefficient = 0, t_Value = 0, P_Value = pval)
      )
    }
  }
}

# -----------------------------------------------------------------------------
# Heatmap construction
# -----------------------------------------------------------------------------
heatmap_matrix <- acast(results_df, Var1 ~ Var2, value.var = "Coefficient")

# Color scale centered at 0, symmetric by max absolute coefficient
mx <- max(abs(heatmap_matrix), na.rm = TRUE)
color_fun <- colorRamp2(c(-mx, 0, mx), c("blue", "white", "red"))

Heatmap(
  heatmap_matrix,
  name = "estimated coefficient",
  col = color_fun,
  show_column_names = TRUE,
  show_row_names = TRUE,
  cluster_rows = TRUE,
  cluster_columns = TRUE
)
