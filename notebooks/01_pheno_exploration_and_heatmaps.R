# 01_pheno_exploration_and_heatmaps.R
# Purpose:
#   Build a merged phenotype + molecular table across brain regions (PFC/CBM/ST),
#   then compute pairwise association strengths and display an interactive clustered heatmap.
#
# Notes:
#   - This is EDA only (not for inference).
#   - Associations are standardized to [-1, 1] where possible.
#   - Categorical / continuous associations use a correlation-ratio style metric (eta),
#     mapped to [-1, 1] by applying sign = sign(Spearman on numeric-coded groups) when feasible.
#

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(janitor)
  library(vcd)        # assocstats for Cramer's V (optional)
  library(plotly)
  library(ggplot2)
  library(reshape2)
})

# ----------------------------
# 0) Helper functions
# ----------------------------

# Safe Cramer's V for factor-factor
cramers_v_safe <- function(x, y) {
  x <- droplevels(as.factor(x))
  y <- droplevels(as.factor(y))
  tbl <- table(x, y)
  # Need at least 2x2
  if (nrow(tbl) < 2 || ncol(tbl) < 2) return(NA_real_)
  out <- tryCatch({
    as.numeric(assocstats(tbl)$cramer)
  }, error = function(e) NA_real_)
  out
}

# Correlation ratio (eta): numeric outcome explained by categorical predictor
# Returns in [0, 1]
eta_ratio <- function(y_numeric, x_factor) {
  y <- y_numeric
  x <- droplevels(as.factor(x_factor))
  ok <- complete.cases(y, x)
  y <- y[ok]
  x <- x[ok]
  if (length(y) < 3) return(NA_real_)
  if (nlevels(x) < 2) return(NA_real_)

  grand_mean <- mean(y)
  ss_between <- sum(tapply(y, x, function(v) length(v) * (mean(v) - grand_mean)^2))
  ss_total <- sum((y - grand_mean)^2)
  if (ss_total == 0) return(NA_real_)
  sqrt(ss_between / ss_total)
}

# Signed eta: eta in [0,1] plus a sign heuristic to map to [-1,1]
# For multi-level factors, sign is not inherently defined; we approximate using Spearman
# between y and ordered group means (works best if factor levels have ordinal meaning).
signed_eta <- function(y_numeric, x_factor) {
  eta <- eta_ratio(y_numeric, x_factor)
  if (is.na(eta)) return(NA_real_)

  # sign heuristic: order levels by mean(y) and correlate with those ordered means
  x <- droplevels(as.factor(x_factor))
  ok <- complete.cases(y_numeric, x)
  y <- y_numeric[ok]
  x <- x[ok]

  means <- tapply(y, x, mean)
  # Map each sample to its group's mean and correlate with y
  yhat <- means[as.character(x)]
  s <- suppressWarnings(cor(y, yhat, method = "spearman"))
  if (is.na(s)) s <- 1
  sign(s) * eta
}

# Spearman correlation for numeric-numeric
spearman_safe <- function(x, y) {
  out <- suppressWarnings(cor(x, y, method = "spearman", use = "pairwise.complete.obs"))
  as.numeric(out)
}

# Utility: scale numeric columns only
scale_numeric_cols <- function(df) {
  num_cols <- names(df)[vapply(df, is.numeric, logical(1))]
  df %>% mutate(across(all_of(num_cols), ~ as.numeric(scale(.x))))
}

# ----------------------------
# 1) Standardize/clean inputs
# ----------------------------

# Ensure consistent key naming
# NOTE: if projid is not character everywhere, merging can silently misbehave
as_key <- function(x) as.character(x)

# ---- PFC core table: DNAm + metadata + clocks ----
pfc_core <- DNAmAge_pfc %>%
  mutate(projid = as_key(projid)) %>%
  select(
    projid,
    Age, Female, pmi, study, gpath,
    plaq_d, plaq_n, nft, count, amyloid,
    starts_with("cogn"), tangles, neurons,
    ca1hipdp, entodp, infpardp, midfrontdp, midtempdp,
    ca1hipnp, entonp, infparnp, midfrontnp, midtempnp,
    ca1hipnft, entonft, infparnft, midfrontnft, midtempnft,
    starts_with("amyloid_"), starts_with("tangles_"),
    ci_num2_gct, ci_num2_mct,
    AD_Reagan, ceradsc, tdp_stage4,
    apoe_genotype, apoe_long, tomm40_hap,
    age_first_ad_dx, age_first_dem_dx, cogdx,
    braaksc, caa_4gp, AD_Clinical,
    PCPhenoAge, dunedinPoAm_git, DNAmClockCortical, Horvath1
  ) %>%
  rename(
    neurons_pfc           = neurons,
    PCPhenoAge_pfc        = PCPhenoAge,
    dunedinPACE_pfc       = dunedinPoAm_git,
    DNAmClockCortical_pfc = DNAmClockCortical,
    Horvath1_pfc          = Horvath1
  )

# ---- Somatic burden ----
som_pfc <- somatic_pfc %>%
  transmute(projid = as_key(INDV), mutational_burden_pfc = N_SITES)

som_cbm <- somatic_cbm %>%
  transmute(projid = as_key(INDV), mutational_burden_cbm = N_SITES)

som_st <- somatic_st %>%
  transmute(projid = as_key(INDV), mutational_burden_st = N_SITES)

# ---- SEM helper: rename SEM columns with region suffix ----
rename_sem <- function(df, region) {
  df %>%
    mutate(projid = as_key(projid)) %>%
    rename_with(~ paste0(.x, "_", region), -projid)
}

pfc_sem <- list(SEM_pfc_1, SEM_pfc_2, SEM_pfc_3) %>%
  lapply(rename_sem, region = "pfc") %>%
  reduce(full_join, by = "projid")

cbm_sem <- list(SEM_cbm_1, SEM_cbm_2, SEM_cbm_3) %>%
  lapply(rename_sem, region = "cbm") %>%
  reduce(full_join, by = "projid")

st_sem <- list(SEM_st_1, SEM_st_2, SEM_st_3) %>%
  lapply(rename_sem, region = "st") %>%
  reduce(full_join, by = "projid")

# Add SEM summaries (hypo/hyper)
add_sem_summaries <- function(df, region) {
  hypo_cols  <- paste0(c("unmeth_hypo", "mid_hypo", "meth_hypo"), "_", region)
  hyper_cols <- paste0(c("unmeth_hyper", "mid_hyper", "meth_hyper"), "_", region)

  df %>%
    mutate(
      !!paste0("sum_hypo_", region)  := rowSums(across(all_of(hypo_cols)),  na.rm = TRUE),
      !!paste0("sum_hyper_", region) := rowSums(across(all_of(hyper_cols)), na.rm = TRUE)
    )
}

pfc_sem <- add_sem_summaries(pfc_sem, "pfc")
cbm_sem <- add_sem_summaries(cbm_sem, "cbm")
st_sem  <- add_sem_summaries(st_sem,  "st")

# ---- CBM DNAm clocks ----
cbm_core <- DNAmAge_cbm %>%
  mutate(projid = as_key(projid)) %>%
  select(projid, neurons, PCPhenoAge, dunedinPoAm_git, DNAmClockCortical, Horvath1) %>%
  rename(
    neurons_cbm           = neurons,
    PCPhenoAge_cbm        = PCPhenoAge,
    dunedinPACE_cbm       = dunedinPoAm_git,
    DNAmClockCortical_cbm = DNAmClockCortical,
    Horvath1_cbm          = Horvath1
  )

# ---- ST DNAm clocks ----
st_core <- DNAmAge_st %>%
  mutate(projid = as_key(projid)) %>%
  select(projid, neurons, PCPhenoAge, dunedinPoAm_git, DNAmClockCortical, Horvath1) %>%
  rename(
    neurons_st           = neurons,
    PCPhenoAge_st        = PCPhenoAge,
    dunedinPACE_st       = dunedinPoAm_git,
    DNAmClockCortical_st = DNAmClockCortical,
    Horvath1_st          = Horvath1
  )

# ---- ROSMAP basic ----
basic <- ROSMAP.BasicPheno %>%
  rename(projid = 1) %>%
  mutate(projid = as_key(projid)) %>%
  select(projid, niareagansc)

# ----------------------------
# 2) Build big_pheno
# ----------------------------
big_pheno <- pfc_core %>%
  full_join(pfc_sem, by = "projid") %>%
  full_join(som_pfc, by = "projid") %>%
  full_join(cbm_sem, by = "projid") %>%
  full_join(som_cbm, by = "projid") %>%
  full_join(cbm_core, by = "projid") %>%
  full_join(st_sem, by = "projid") %>%
  full_join(som_st, by = "projid") %>%
  full_join(st_core, by = "projid") %>%
  left_join(basic, by = "projid") %>%
  clean_names()

# ----------------------------
# 3) Derived phenotype variables
# ----------------------------
big_pheno <- big_pheno %>%
  mutate(
    niareagansc = as.factor(niareagansc),

    tdp_stage4_recode = case_when(
      tdp_stage4 %in% c(0, 1) ~ FALSE,
      tdp_stage4 %in% c(2, 3) ~ TRUE,
      TRUE ~ NA
    ) %>% as.factor(),

    # NIA-Reagan dichotomy:
    niareagansc_recode = case_when(
      niareagansc %in% c(1, 2) ~ 1,
      niareagansc %in% c(3, 4) ~ 0,
      TRUE ~ NA_real_
    ) %>% as.factor(),

    length_ad = age - age_first_ad_dx,

    # AD indicator: strict clinical+pathology
    ad_indicator_binary = case_when(
      cogdx == 3 & niareagansc %in% c(1, 2) ~ 1,
      cogdx == 1 ~ 0,
      TRUE ~ NA_real_
    ) %>% as.factor(),

    # AD indicator: any cognitive impairment + pathology
    ad_indicator_binary_anycog = case_when(
      cogdx %in% c(2, 3) & niareagansc %in% c(1, 2) ~ 1,
      cogdx == 1 ~ 0,
      TRUE ~ NA_real_
    ) %>% as.factor()
  )

# ----------------------------
# 4) Scale numeric columns for comparability (optional)
# ----------------------------
big_pheno_scaled <- scale_numeric_cols(big_pheno)

# ----------------------------
# 5) Pairwise association matrix
# ----------------------------
vars <- names(big_pheno_scaled)

assoc_matrix <- matrix(NA_real_, nrow = length(vars), ncol = length(vars),
                       dimnames = list(vars, vars))

for (i in seq_along(vars)) {
  for (j in seq_along(vars)) {
    if (i == j) next

    xi <- big_pheno_scaled[[vars[i]]]
    xj <- big_pheno_scaled[[vars[j]]]

    if (is.numeric(xi) && is.numeric(xj)) {
      assoc_matrix[i, j] <- spearman_safe(xi, xj)

    } else if (is.factor(xi) && is.factor(xj)) {
      assoc_matrix[i, j] <- cramers_v_safe(xi, xj)

    } else if (is.numeric(xi) && is.factor(xj)) {
      assoc_matrix[i, j] <- signed_eta(xi, xj)

    } else if (is.factor(xi) && is.numeric(xj)) {
      assoc_matrix[i, j] <- signed_eta(xj, xi)

    } else {
      assoc_matrix[i, j] <- NA_real_
    }
  }
}

# For clustering, replace NA with 0 (EDA choice; avoids dist() errors)
assoc_for_cluster <- assoc_matrix
assoc_for_cluster[is.na(assoc_for_cluster)] <- 0

hc_rows <- hclust(dist(assoc_for_cluster))
hc_cols <- hclust(dist(t(assoc_for_cluster)))

ordered_matrix <- assoc_matrix[hc_rows$order, hc_cols$order]

# ----------------------------
# 6) Plot interactive heatmap
# ----------------------------
assoc_long <- reshape2::melt(ordered_matrix)

p <- ggplot(assoc_long, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue", high = "red", mid = "white",
    midpoint = 0, limits = c(-1, 1), na.value = "grey90",
    name = "Association"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank()
  )

interactive_heatmap <- ggplotly(p, tooltip = c("Var1", "Var2", "value"))
interactive_heatmap
