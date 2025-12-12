# =============================================================================
# 01_pheno_exploration_and_heatmaps.R
# Yaro
# last update: 12/6/23
#
# PURPOSE
# -------
# This script merges several phenotype/omics-derived summary tables across 3 brain
# regions (PFC, CBM, ST), including:
#   - Stochastic epigenetic mutations (SEM) features (SEM_* tables)
#   - DNA methylation aging clocks (DNAmAge_* tables)
#   - Somatic mutational burden per region (somatic_* tables)
#   - Core ROSMAP phenotypes (ROSMAP.BasicPheno)
#
# Then it computes a pairwise "association strength" matrix between ALL columns:
#   - numeric vs numeric: correlation via cor() 
#   - factor vs factor: Cramer's V of chisq.test()
#   - mixed: correlation of as.numeric() codes 
#
# Finally, it hierarchically clusters the association matrix and visualizes it as an
# interactive plotly heatmap (ggplotly over ggplot2 tiles).


library(confintr)  # cramersv
library(plotly)    # ggplotly
library(ggplot2)   # heatmap base plot
library(reshape2)  # melt

# -----------------------------------------------------------------------------
# Small helpers 
# -----------------------------------------------------------------------------

# Lightweight wrapper so repetitive merge() calls are easier to scan.
mmerge <- function(x, y, ...) merge(x, y, ...)

# Index-based renaming
rename_by_index <- function(df, mapping) {
  # mapping is a named integer vector, e.g. c(meth_hypo_pfc=4, neurons_pfc=30)
  for (nm in names(mapping)) colnames(df)[mapping[[nm]]] <- nm
  df
}

# Association function
assoc_value <- function(x, y) {
  if (is.numeric(x) && is.numeric(y)) {
    cor(x, y, use = "complete.obs")
  } else if (is.factor(x) && is.factor(y)) {
    cramersv(chisq.test(x, y))
  } else {
    #  mixed types -> cor(as.numeric(x), as.numeric(y))
    cor(as.numeric(x), as.numeric(y), use = "complete.obs")
  }
}

# -----------------------------------------------------------------------------
# Column selection for DNAmAge_pfc merge
# -----------------------------------------------------------------------------

dnamage_pfc_cols <- c(
  "PCPhenoAge", "dunedinPoAm_git", "DNAmClockCortical", "Horvath1",
  "projid", "Age", "Female", "pmi", "study", "gpath",
  "plaq_d", "plaq_n", "nft", "count", "amyloid",
  "cognep_random_slope", "cogng_random_slope", "cognpo_random_slope",
  "cognps_random_slope", "cognse_random_slope", "cognwo_random_slope",
  "tangles", "neurons",
  "ca1hipdp", "entodp", "infpardp", "midfrontdp", "midtempdp",
  "ca1hipnp", "entonp", "infparnp", "midfrontnp", "midtempnp",
  "ca1hipnft", "entonft", "infparnft", "midfrontnft", "midtempnft",
  "amyloid_ag", "amyloid_calc", "amyloid_cg", "amyloid_ec", "amyloid_hip",
  "amyloid_it", "amyloid_mf", "amyloid_mt", "amyloid_sf",
  "tangles_ag", "tangles_calc", "tangles_cg", "tangles_ec",
  "tangles_hip", "tangles_it", "tangles_mf", "tangles_mt",
  "ci_num2_gct", "ci_num2_mct",
  "AD_Reagan", "ceradsc", "tdp_stage4", "apoe_genotype", "apoe_long",
  "tomm40_hap", "age_first_ad_dx", "age_first_dem_dx", "cogdx",
  "braaksc", "caa_4gp", "AD_Clinical"
)

# =============================================================================
# CLEANING / MERGING PFC
# =============================================================================

# Merge SEM_pfc_1 with selected DNAmAge_pfc columns by projid
big_pheno <- mmerge(SEM_pfc_1, DNAmAge_pfc[dnamage_pfc_cols], by = "projid")

# Add SEM_pfc_2 and SEM_pfc_3 (in that order)
big_pheno <- mmerge(SEM_pfc_2, big_pheno, by = "projid")
big_pheno <- mmerge(SEM_pfc_3, big_pheno, by = "projid")

# Region-specific SEM-derived summary features
big_pheno$sum_hypo_pfc  <- big_pheno$unmeth_hypo  + big_pheno$mid_hypo  + big_pheno$meth_hypo
big_pheno$sum_hyper_pfc <- big_pheno$unmeth_hyper + big_pheno$mid_hyper + big_pheno$meth_hyper

# Merge somatic burden 
temp <- somatic_pfc[c("INDV", "N_SITES")]
colnames(temp) <- c("projid", "mutational_burden_pfc")
big_pheno <- mmerge(temp, big_pheno, by = "projid", all.y = TRUE)

# Index-based renames 
big_pheno <- rename_by_index(big_pheno, c(
  meth_hypo_pfc         = 4,
  mid_hyper_pfc         = 5,
  mid_hypo_pfc          = 6,
  unmeth_hyper_pfc      = 7,
  PCPhenoAge_pfc        = 9,
  dunedinPACE_pfc       = 10,
  DNAmClockCortical_pfc = 11,
  Horvath1_pfc          = 12,
  neurons_pfc           = 30
))

# Optional non-invasive guardrail: confirm those indices still exist
stopifnot(ncol(big_pheno) >= 30)

# =============================================================================
# CLEANING / MERGING CBM
# =============================================================================

big_pheno <- mmerge(SEM_cbm_1, big_pheno, by = "projid")
big_pheno <- mmerge(SEM_cbm_2, big_pheno, by = "projid")
big_pheno <- mmerge(SEM_cbm_3, big_pheno, by = "projid")

big_pheno$sum_hypo_cbm  <- big_pheno$unmeth_hypo.x + big_pheno$mid_hypo  + big_pheno$meth_hypo
big_pheno$sum_hyper_cbm <- big_pheno$unmeth_hyper  + big_pheno$mid_hyper + big_pheno$meth_hyper.x

temp <- somatic_cbm[c("INDV", "N_SITES")]
colnames(temp) <- c("projid", "mutational_burden_cbm")
big_pheno <- mmerge(temp, big_pheno, by = "projid", all.y = TRUE)

big_pheno <- mmerge(
  big_pheno,
  DNAmAge_cbm[c("neurons", "projid", "PCPhenoAge", "dunedinPoAm_git", "DNAmClockCortical", "Horvath1")],
  by = "projid"
)

big_pheno <- rename_by_index(big_pheno, c(
  meth_hypo_cbm         = 4,
  mid_hyper_cbm         = 5,
  mid_hypo_cbm          = 6,
  unmeth_hyper_cbm      = 7,
  neurons_cbm           = 88,
  PCPhenoAge_cbm        = 89,
  dunedinPACE_cbm       = 90,
  DNAmClockCortical_cbm = 91,
  Horvath1_cbm          = 92
))

stopifnot(ncol(big_pheno) >= 92)

# =============================================================================
# MERGE SEM ST + SOMATIC ST (then merge DNAmAge_st)
# =============================================================================

big_pheno <- mmerge(SEM_st_1, big_pheno, by = "projid")
big_pheno <- mmerge(SEM_st_2, big_pheno, by = "projid")
big_pheno <- mmerge(SEM_st_3, big_pheno, by = "projid")

big_pheno$sum_hypo_st  <- big_pheno$unmeth_hypo  + big_pheno$mid_hypo  + big_pheno$meth_hypo
big_pheno$sum_hyper_st <- big_pheno$unmeth_hyper + big_pheno$mid_hyper + big_pheno$meth_hyper

temp <- somatic_st[c("INDV", "N_SITES")]
colnames(temp) <- c("projid", "mutational_burden_st")
big_pheno <- mmerge(temp, big_pheno, by = "projid", all.y = TRUE)
rm(temp)

# Merge DNAmAge_st LAST for ST region
big_pheno <- mmerge(
  big_pheno,
  DNAmAge_st[c("neurons", "projid", "PCPhenoAge", "dunedinPoAm_git", "DNAmClockCortical", "Horvath1")],
  by = "projid"
)

big_pheno <- rename_by_index(big_pheno, c(
  meth_hypo_st         = 4,
  mid_hyper_st         = 5,
  mid_hypo_st          = 6,
  unmeth_hyper_st      = 7,
  neurons_st           = 102,
  PCPhenoAge_st        = 103,
  dunedinPACE_st       = 104,
  DNAmClockCortical_st = 105,
  Horvath1_st          = 106
))

stopifnot(ncol(big_pheno) >= 106)

# =============================================================================
# CLEANING THE MERGED DATA (phenotype recodes + scaling)
# =============================================================================

colnames(ROSMAP.BasicPheno)[1] <- "projid"
big_pheno <- mmerge(big_pheno, ROSMAP.BasicPheno[c("projid", "niareagansc")], by = "projid")
big_pheno$niareagansc <- as.factor(big_pheno$niareagansc)

# Remove redundant columns created by merge suffixing 
big_pheno <- big_pheno[, !(names(big_pheno) %in% c(
  "meth_hyper", "unmeth_hypo",
  "meth_hyper.x", "unmeth_hypo.x",
  "meth_hyper.y", "unmeth_hypo.y"
))]

# Recode TDP stage into binary factor 
big_pheno$tdp_stage4_recode <- ifelse(big_pheno$tdp_stage4 %in% c(0, 1), FALSE,
                                      ifelse(big_pheno$tdp_stage4 %in% c(2, 3), TRUE, NA))
big_pheno$tdp_stage4_recode <- as.factor(big_pheno$tdp_stage4_recode)

# Recode niareagansc into binary factor 
big_pheno$niareagansc_recode <- ifelse(big_pheno$niareagansc %in% c(1, 2), 1,
                                       ifelse(big_pheno$niareagansc %in% c(3, 4), 0, NA))
big_pheno$niareagansc_recode <- as.factor(big_pheno$niareagansc_recode)

# Disease duration estimate 
big_pheno$lengthAD <- big_pheno$Age - big_pheno$age_first_ad_dx

# AD indicator: strict AD (clinical AD + neuropathologic AD)
big_pheno$AD_indicator_binary <- ifelse(big_pheno$cogdx == 3 & big_pheno$niareagansc %in% c(1, 2), 1,
                                        ifelse(big_pheno$cogdx == 1, 0, NA))
big_pheno$AD_indicator_binary <- as.factor(big_pheno$AD_indicator_binary)

# AD indicator: "any cognitive impairment" + neuropathologic AD
big_pheno$AD_indicator_binary_anycog <- ifelse(big_pheno$cogdx %in% c(2, 3) & big_pheno$niareagansc %in% c(1, 2), 1,
                                               ifelse(big_pheno$cogdx == 1, 0, NA))
big_pheno$AD_indicator_binary_anycog <- as.factor(big_pheno$AD_indicator_binary_anycog)

# Scale numeric features only 
numeric_columns <- colnames(big_pheno)[sapply(big_pheno, is.numeric)]
big_pheno_scaled <- big_pheno
big_pheno_scaled[, numeric_columns] <- scale(big_pheno[, numeric_columns])

# =============================================================================
# ASSOCIATION MATRIX (pairwise over all columns)
# =============================================================================

assoc_matrix <- matrix(
  NA,
  ncol = ncol(big_pheno_scaled),
  nrow = ncol(big_pheno_scaled),
  dimnames = list(names(big_pheno_scaled), names(big_pheno_scaled))
)

for (i in 1:ncol(big_pheno_scaled)) {
  for (j in 1:ncol(big_pheno_scaled)) {
    if (i != j) {
      assoc_matrix[i, j] <- assoc_value(big_pheno_scaled[[i]], big_pheno_scaled[[j]])
      print(paste(colnames(big_pheno_scaled)[i], colnames(big_pheno_scaled)[j], assoc_matrix[i, j]))
    }
  }
}

# =============================================================================
# CLUSTERING + INTERACTIVE HEATMAP
# =============================================================================

hc_rows <- hclust(dist(t(assoc_matrix)))
hc_cols <- hclust(dist(assoc_matrix))

ordered_matrix <- assoc_matrix[hc_rows$order, hc_cols$order]
assoc_long <- reshape2::melt(ordered_matrix)

p <- ggplot(assoc_long, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue", high = "red", mid = "white",
    midpoint = 0, limit = c(-1, 1), space = "Lab",
    name = "Association"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank()
  )

interactive_heatmap <- ggplotly(p, tooltip = c("Var1", "Var2", "value"))

# Launch interactive heatmap
interactive_heatmap
