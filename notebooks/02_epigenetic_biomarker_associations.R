# Yaro
# Last update: 12/6/23
#
# Script Purpose:
#   Test associations between a set of biomarkers ("clockColumns") and a set of
#   pathology / cognitive variables ("cognitive_cols") using regression models.
#   Save effect size + p-value for each biomarker–outcome pair and visualize
#   the resulting association matrix as a heatmap.
#
# Output:
#   - `results_df`: long table with coefficient and p-value for each pair.
#   - A ComplexHeatmap heatmap of coefficients (optionally masked by significance).

# ---- Libraries ----
library(dplyr)          # data manipulation
library(ggplot2)        # optional visualization (not strictly used here)
library(ComplexHeatmap) # heatmap
library(reshape2)       # acast
library(circlize)       # colorRamp2
library(ordinal)        # clm (if using ordinal outcomes)

# ---- User settings ----

# Biomarker variables to test (predictors)
clockColumns <- c(
  "mutational_burden_st", "mutational_burden_pfc", "mutational_burden_cbm",
  "sum_hypo_pfc", "sum_hypo_st", "sum_hypo_cbm", "sum_hyper_pfc", "sum_hyper_cbm", "sum_hyper_st",
  "meth_hypo_pfc", "meth_hypo_cbm", "meth_hypo_st", "mid_hyper_pfc", "mid_hyper_cbm", "mid_hyper_st",
  "mid_hypo_pfc", "mid_hypo_cbm", "mid_hypo_st", "unmeth_hyper_pfc", "unmeth_hyper_cbm", "unmeth_hyper_st",
  "PCPhenoAge_pfc", "PCPhenoAge_cbm", "PCPhenoAge_st",
  "Horvath1_pfc", "Horvath1_cbm", "Horvath1_st", "dunedinPACE_pfc", "dunedinPACE_cbm", "dunedinPACE_st",
  "DNAmClockCortical_cbm", "DNAmClockCortical_pfc", "DNAmClockCortical_st"
)

# Covariates included in each model.
# Note: this is a fixed set; if you want region-matched neuron proportion, adapt as needed.
covariates <- c("Age", "Female", "pmi") # eg "neurons_st"

# Choose the outcome type for THIS run.
# Valid values: "numeric", "binary", "ordinal"
outcome_type <- "binary"

# Define candidate outcomes for each type.
# You should set `cognitive_cols` ONCE based on outcome_type.
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

cognitive_cols_ordinal <- c("tdp_stage4", "tomm40_hap", "AD_Clinical", "cogdx", "braaksc")

cognitive_cols_binary <- c("ci_num2_gct", "ci_num2_mct", "AD_Reagan", "ceradsc", "tdp_stage4_recode", "apoe_genotype")

cognitive_cols <- switch(
  outcome_type,
  numeric = cognitive_cols_numeric,
  binary  = cognitive_cols_binary,
  ordinal = cognitive_cols_ordinal
)

# Significance threshold for optional masking of coefficients
alpha <- 0.05

# If TRUE, non-significant coefficients are replaced with 0 in the heatmap.
# If FALSE, the true coefficients are shown regardless of significance.
mask_nonsignificant <- TRUE


# ---- Initialize results storage ----
# We store one row per (biomarker, outcome) pair.
results_df <- data.frame(
  Var1 = character(),    # biomarker
  Var2 = character(),    # outcome
  Coefficient = numeric(),
  P_Value = numeric(),
  Stat = numeric(),      # t-value (lm) or z-value (glm/clm)
  stringsAsFactors = FALSE
)

# ---- Helper function: fit model and extract coefficient/pvalue for biomarker ----
fit_and_extract <- function(dat, outcome, biomarker, covariates, outcome_type) {

  # Build the model formula
  formula_str <- paste(outcome, "~", biomarker, "+", paste(covariates, collapse = " + "))
  fml <- as.formula(formula_str)

  # Fit model depending on outcome type
  if (outcome_type == "numeric") {
    fit <- lm(fml, data = dat)
    coefs <- summary(fit)$coefficients
  } else if (outcome_type == "binary") {
    fit <- glm(fml, data = dat, family = "binomial")
    coefs <- summary(fit)$coefficients
  } else if (outcome_type == "ordinal") {
    # Ensure outcome is ordered
    dat[[outcome]] <- as.ordered(dat[[outcome]])
    fit <- clm(fml, data = dat)
    coefs <- summary(fit)$coefficients
  } else {
    stop("Unknown outcome_type: ", outcome_type)
  }

  # If the biomarker term is missing (e.g., dropped due to singularity), return NA
  if (!(biomarker %in% rownames(coefs))) {
    return(list(beta = NA_real_, stat = NA_real_, p = NA_real_, formula = formula_str))
  }

  beta <- coefs[biomarker, 1]
  stat <- coefs[biomarker, 3]
  pval <- coefs[biomarker, 4]

  list(beta = beta, stat = stat, p = pval, formula = formula_str)
}


# ---- Main loop over biomarker–outcome pairs ----
for (bio in clockColumns) {
  for (out in cognitive_cols) {

    # Skip if columns are missing (defensive; helpful when reusing script)
    needed <- c(bio, out, covariates)
    if (!all(needed %in% colnames(big_pheno_scaled))) {
      results_df <- rbind(
        results_df,
        data.frame(Var1 = bio, Var2 = out, Coefficient = NA_real_, P_Value = NA_real_, Stat = NA_real_)
      )
      next
    }

    # Fit model + extract statistics
    res <- fit_and_extract(
      dat = big_pheno_scaled,
      outcome = out,
      biomarker = bio,
      covariates = covariates,
      outcome_type = outcome_type
    )

    # Optionally mask non-significant coefficients
    beta_to_store <- res$beta
    if (mask_nonsignificant && !is.na(res$p) && res$p >= alpha) {
      beta_to_store <- 0
    }

    results_df <- rbind(
      results_df,
      data.frame(
        Var1 = bio,
        Var2 = out,
        Coefficient = beta_to_store,
        P_Value = res$p,
        Stat = res$stat,
        stringsAsFactors = FALSE
      )
    )
  }
}

# ---- Convert results to a matrix for heatmap ----
# Rows = biomarkers (Var1), Columns = outcomes (Var2)
heatmap_matrix <- acast(results_df, Var1 ~ Var2, value.var = "Coefficient")

# ---- Define heatmap colors ----
# Centered at 0. Use symmetric limits so blue/red are comparable.
mx <- max(abs(heatmap_matrix), na.rm = TRUE)
if (!is.finite(mx) || mx == 0) mx <- 1

color_fun <- colorRamp2(c(-mx, 0, mx), c("blue", "white", "red"))

# ---- Plot heatmap ----
Heatmap(
  heatmap_matrix,
  name = "estimated coefficient",
  col = color_fun,
  show_column_names = TRUE,
  show_row_names = TRUE,
  cluster_rows = TRUE,
  cluster_columns = TRUE
)
