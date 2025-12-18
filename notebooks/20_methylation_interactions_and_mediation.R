# =============================================================================
# 20_methylation_interactions_and_mediation.R
#
# Expected objects in environment:
#   - DNAm_all: matrix/data.frame [samples x CpGs] 
#   - big_pheno_all: data.frame with at least:
#       SampleID (optional but recommended), niareagansc_recode, apoe_genotype,
#       Female, Age, pmi, neurons_pfc, braaksc_num, ceradsc_num
#   - (optional) WGCNA: MEs_DNAme, moduleColors_DNAme, phenotypes (trait names)
#
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(parallel)
  library(mediation)
  library(MASS)
})

# -----------------------------------------------------------------------------
# 0) Harmonize phenotype variables
# -----------------------------------------------------------------------------

# AD outcome: must be 0/1 for binomial glm
big_pheno_all$AD_status <- as.integer(as.character(big_pheno_all$niareagansc_recode))
if (any(is.na(big_pheno_all$AD_status))) {
  # try factor levels fallback
  big_pheno_all$AD_status <- as.integer(big_pheno_all$niareagansc_recode) - 1
}

# APOE4 exposure: 0/1
big_pheno_all$APOE4 <- as.integer(as.character(big_pheno_all$apoe_genotype))
if (any(is.na(big_pheno_all$APOE4))) {
  big_pheno_all$APOE4 <- as.integer(big_pheno_all$apoe_genotype) - 1
}

# Covariates: handle missing
big_pheno_all$pmi[is.na(big_pheno_all$pmi)] <- mean(big_pheno_all$pmi, na.rm = TRUE)

# Basic alignment check
stopifnot(nrow(DNAm_all) == nrow(big_pheno_all))

# -----------------------------------------------------------------------------
# 1) CpG-wise APOE4 × CpG interaction (logistic AD)
# -----------------------------------------------------------------------------

# Robust coefficient extractor (doesn't assume exact names)
extract_term <- function(coefs, pattern) {
  rn <- rownames(coefs)
  idx <- grep(pattern, rn)
  if (length(idx) == 0) return(c(Estimate=NA, `Std. Error`=NA, `Pr(>|z|)`=NA))
  # If multiple matches, take first
  coefs[idx[1], c("Estimate", "Std. Error", "Pr(>|z|)")]
}

perform_interaction_one_cpg <- function(cpg) {
  mol <- DNAm_all[, cpg]

  df <- data.frame(
    AD_status = big_pheno_all$AD_status,
    APOE4 = big_pheno_all$APOE4,
    Molecule = as.numeric(mol),
    Female = big_pheno_all$Female,
    Age = big_pheno_all$Age,
    pmi = big_pheno_all$pmi,
    Neurons = big_pheno_all$neurons_pfc
  )

  # Drop incomplete rows 
  df <- df[complete.cases(df), ]

  # AD logistic
  fit <- glm(
    AD_status ~ APOE4 * Molecule + Female + Age + pmi + Neurons,
    data = df,
    family = binomial()
  )
  coefs <- summary(fit)$coefficients

  apoe  <- extract_term(coefs, "^APOE4$")
  molt  <- extract_term(coefs, "^Molecule$")
  inter <- extract_term(coefs, "APOE4:Molecule|Molecule:APOE4")

  data.frame(
    CpG = cpg,
    Beta_APOE = apoe["Estimate"], SE_APOE = apoe["Std. Error"], P_APOE = apoe["Pr(>|z|)"],
    Beta_Molecule = molt["Estimate"], SE_Molecule = molt["Std. Error"], P_Molecule = molt["Pr(>|z|)"],
    Beta_Interaction = inter["Estimate"], SE_Interaction = inter["Std. Error"], P_Interaction = inter["Pr(>|z|)"],
    N = nrow(df),
    stringsAsFactors = FALSE
  )
}

# Parallel run
run_interactions_parallel <- function(cpg_list, ncores = max(1, detectCores() - 1)) {
  cl <- makeCluster(ncores)
  clusterExport(cl, c("DNAm_all", "big_pheno_all", "perform_interaction_one_cpg", "extract_term"), envir = environment())
  clusterEvalQ(cl, {
    library(stats)
  })
  res <- parLapply(cl, cpg_list, function(cpg) {
    tryCatch(perform_interaction_one_cpg(cpg), error = function(e) NULL)
  })
  stopCluster(cl)
  bind_rows(res)
}

methylation_interaction_results <- run_interactions_parallel(colnames(DNAm_all))

methylation_interaction_results <- methylation_interaction_results %>%
  mutate(Adj_P_Interaction = p.adjust(P_Interaction, method = "BH"))


# Example: significant interactions
significant_interactions <- methylation_interaction_results %>%
  filter(!is.na(Adj_P_Interaction) & Adj_P_Interaction <= 0.10)

# -----------------------------------------------------------------------------
# 2) Mediation: APOE4 -> CpG -> AD (binary)
# -----------------------------------------------------------------------------

perform_mediation_one_cpg <- function(cpg, sims = 1000) {
  df <- data.frame(
    APOE4 = big_pheno_all$APOE4,
    Methylation = as.numeric(DNAm_all[, cpg]),
    AD = big_pheno_all$AD_status,
    Female = big_pheno_all$Female,
    Age = big_pheno_all$Age,
    pmi = big_pheno_all$pmi,
    neurons_pfc = big_pheno_all$neurons_pfc
  )
  df <- df[complete.cases(df), ]

  med_model <- lm(Methylation ~ APOE4 + Female + Age + pmi + neurons_pfc, data = df)
  out_model <- glm(AD ~ APOE4 + Methylation + Female + Age + pmi + neurons_pfc,
                   family = binomial(), data = df)

  mediate(med_model, out_model,
          treat = "APOE4", mediator = "Methylation",
          boot = TRUE, sims = sims)
}

run_mediation_parallel <- function(cpg_list, sims = 1000, ncores = max(1, detectCores() - 2)) {
  cl <- makeCluster(ncores)
  clusterExport(cl, c("DNAm_all", "big_pheno_all", "perform_mediation_one_cpg", "sims"), envir = environment())
  clusterEvalQ(cl, {
    library(mediation)
    library(stats)
  })
  res <- parLapply(cl, cpg_list, function(cpg) {
    tryCatch(perform_mediation_one_cpg(cpg, sims = sims), error = function(e) NULL)
  })
  stopCluster(cl)
  setNames(res, cpg_list)
}

# Example: run on a subset — e.g., significant CpGs list:
# mediation_results <- run_mediation_parallel(cpgs_to_test, sims = 1000)

extract_mediation_table <- function(med_list, p_cut = 0.10) {
  out <- lapply(names(med_list), function(cpg) {
    mr <- med_list[[cpg]]
    if (is.null(mr)) return(NULL)
    s <- summary(mr)

    # ACME = average causal mediation effect
    data.frame(
      CpG = cpg,
      ACME = mr$d.avg,
      ACME_p = mr$d.avg.p,
      ADE = mr$z.avg,
      ADE_p = mr$z.avg.p,
      TotalEffect = mr$tau.coef,
      TotalEffect_p = mr$tau.p,
      PropMediated = mr$n.avg,
      PropMediated_p = mr$n.avg.p,
      stringsAsFactors = FALSE
    )
  })
  bind_rows(out) %>%
    filter(!is.na(ACME_p) & ACME_p <= p_cut)
}

# -----------------------------------------------------------------------------
# 3) Moderated mediation (APOE4 × CpG in outcome model)
# -----------------------------------------------------------------------------

perform_moderated_mediation_one_cpg <- function(cpg, sims = 1000) {
  df <- data.frame(
    APOE4 = big_pheno_all$APOE4,
    Methylation = as.numeric(scale(DNAm_all[, cpg])),
    AD = big_pheno_all$AD_status,
    Female = big_pheno_all$Female,
    Age = big_pheno_all$Age,
    pmi = big_pheno_all$pmi,
    neurons_pfc = big_pheno_all$neurons_pfc
  )
  df <- df[complete.cases(df), ]

  med_model <- lm(Methylation ~ APOE4 + Female + Age + pmi + neurons_pfc, data = df)
  out_model <- glm(AD ~ APOE4 * Methylation + Female + Age + pmi + neurons_pfc,
                   family = binomial(), data = df)

  mediate(med_model, out_model,
          treat = "APOE4", mediator = "Methylation",
          boot = TRUE, sims = sims)
}

# Example single CpG:
# mm_res <- perform_moderated_mediation_one_cpg("cg06329447", sims = 1000)
# print(summary(mm_res))

# -----------------------------------------------------------------------------
# 4) Association of a CpG with neuropathology (Braak/CERAD)
# -----------------------------------------------------------------------------
# If Braak/CERAD are ordinal, use  polr. If binary, use glm binomial.
test_neuropath_polr <- function(cpg, outcome_col) {
  df <- data.frame(
    outcome = big_pheno_all[[outcome_col]],
    CpG = as.numeric(scale(DNAm_all[, cpg])),
    Female = big_pheno_all$Female,
    Age = big_pheno_all$Age,
    pmi = big_pheno_all$pmi,
    neurons_pfc = big_pheno_all$neurons_pfc
  )
  df <- df[complete.cases(df), ]
  df$outcome <- as.ordered(df$outcome)

  fit <- MASS::polr(outcome ~ CpG + Female + Age + pmi + neurons_pfc,
                    data = df, Hess = TRUE)

  # Wald p for CpG
  co <- coef(summary(fit))
  p <- 2 * pnorm(abs(co["CpG", "t value"]), lower.tail = FALSE)

  list(
    outcome = outcome_col,
    beta_CpG = co["CpG", "Value"],
    se_CpG = co["CpG", "Std. Error"],
    p_CpG = p,
    n = nrow(df)
  )
}

# Example:
# braak_res <- test_neuropath_polr("cg06329447", "braaksc_num")
# cerad_res <- test_neuropath_polr("cg06329447", "ceradsc_num")
# print(braak_res); print(cerad_res)

# -----------------------------------------------------------------------------
# 5) WGCNA methylation module interactions (optional)
# -----------------------------------------------------------------------------
# Requires: phenotypes (vector of trait column names), MEs_DNAme (data.frame-like),
#           and big_pheno_all has those traits and apoe_genotype
if (exists("phenotypes") && exists("MEs_DNAme")) {

  interaction_results_WGCNA_methylation <- list()

  for (trait in phenotypes) {
    for (module in colnames(MEs_DNAme)) {

      df <- cbind(big_pheno_all, MEs_DNAme)
      df <- df[complete.cases(df[, c(trait, module, "apoe_genotype", "Female", "Age")]), ]

      # Choose model family based on trait type
      y <- df[[trait]]
      is_binary <- length(unique(y)) == 2

      if (is_binary) {
        fit <- glm(as.formula(paste(trait, "~", module, "* apoe_genotype + Female + Age")),
                   data = df, family = binomial())
      } else {
        fit <- lm(as.formula(paste(trait, "~", module, "* apoe_genotype + Female + Age")),
                  data = df)
      }

      sm <- summary(fit)
      coef_name <- paste0(module, ":apoe_genotypeTRUE")

      if (coef_name %in% rownames(sm$coefficients)) {
        interaction_results_WGCNA_methylation[[paste(module, trait, sep = "_")]] <-
          c(sm$coefficients[coef_name, ], module = module, trait = trait)
      }
    }
  }

  interaction_results_WGCNA_methylation_df <- as.data.frame(do.call(rbind, interaction_results_WGCNA_methylation))
  interaction_results_WGCNA_methylation_df[, 1:4] <- sapply(interaction_results_WGCNA_methylation_df[, 1:4], as.numeric)

  methylation_significant_module_interactions <-
    interaction_results_WGCNA_methylation_df[interaction_results_WGCNA_methylation_df$`Pr(>|t|)` <= 0.1 |
                                              interaction_results_WGCNA_methylation_df$`Pr(>|z|)` <= 0.1, ]
}

# -----------------------------------------------------------------------------
# 6) Module mediation (optional)
# -----------------------------------------------------------------------------
# APOE4 -> module eigengene -> AD
perform_mediation_module <- function(module_vec, sims = 1000) {
  df <- data.frame(
    APOE4 = big_pheno_all$APOE4,
    Mediator = as.numeric(module_vec),
    AD = big_pheno_all$AD_status,
    Female = big_pheno_all$Female,
    Age = big_pheno_all$Age,
    pmi = big_pheno_all$pmi,
    neurons_pfc = big_pheno_all$neurons_pfc
  )
  df <- df[complete.cases(df), ]

  med_model <- lm(Mediator ~ APOE4 + Female + Age + pmi + neurons_pfc, data = df)
  out_model <- glm(AD ~ APOE4 + Mediator + Female + Age + pmi + neurons_pfc,
                   family = binomial(), data = df)

  mediate(med_model, out_model, treat = "APOE4", mediator = "Mediator",
          boot = TRUE, sims = sims)
}

# -----------------------------------------------------------------------------
# 7) Enrichment with missMethyl (optional)
# -----------------------------------------------------------------------------
# Example usage:
#   significant_cpgs <- significant_interactions$CpG
#
# Requires:
#   library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
#   library(missMethyl)
#
# gst <- gometh(sig.cpg = significant_cpgs, all.cpg = colnames(DNAm_all), array.type = "EPIC")
# topGSA(gst[gst$FDR <= 0.05, ])
