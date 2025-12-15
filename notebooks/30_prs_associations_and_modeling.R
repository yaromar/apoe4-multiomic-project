# =============================================================================
# 30_prs_associations_and_modeling.R
#
# Purpose:
#   1) Load PRS and merge into phenotype objects
#   2) PRS associations with CpG + proteins
#   3) PRS associations with AD/neuropath outcomes
#   4) Predictive modeling comparison (PRS vs DNAme vs Proteomics) with repeated CV
#      and ROC/AUC stratified by APOE4
#
# Expected objects already in environment:
#   - big_pheno_all: phenotype df keyed by projid (for DNAm + AD outcomes)
#   - pheno_PFC: phenotype df keyed by projid (for proteomics subjects)
#   - DNAm_all: matrix/data.frame with rownames = projid
#   - proteins_PFC_scaled: df/matrix with rownames = projid
#   - datExpr: DNAm expression matrix for WGCNA eigengenes, rownames = projid
#   - moduleColors_DNAme: DNAm module colors for WGCNA
#
#
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(broom)
  library(pROC)
  library(tibble)
  library(purrr)
  library(caret)
  library(ggplot2)
  library(gridExtra)
  library(WGCNA)
})

# -----------------------------------------------------------------------------
# 0) Helpers
# -----------------------------------------------------------------------------

read_prs_file <- function(path) {
  prs <- read.table(path, header = TRUE, stringsAsFactors = FALSE)
  stopifnot(all(c("IID") %in% colnames(prs)))
  prs %>%
    mutate(projid = substring(IID, 4)) %>%           # remove first 3 chars
    mutate(
      PRScs_z = as.numeric(scale(PRScs_score)),
      SDPR_z  = as.numeric(scale(SDPR_score))
    ) %>%
    dplyr::select(projid, PRScs_z, SDPR_z)
}

ensure_binary_outcomes <- function(df) {
  # AD binary and label factor for caret
  out <- df

  # AD indicator (0/1)
  if (!("AD_bin" %in% names(out))) {
    out$AD_bin <- as.integer(as.character(out$niareagansc_recode))
    if (any(is.na(out$AD_bin))) out$AD_bin <- as.integer(out$niareagansc_recode) - 1
  }

  # Factor with stable levels: Control first, AD second (positive class = "AD")
  out$AD_status <- factor(ifelse(out$AD_bin == 1, "AD", "Control"),
                          levels = c("Control", "AD"))

  # APOE4 carrier (0/1)
  out$APOE4 <- as.integer(as.character(out$apoe_genotype))
  if (any(is.na(out$APOE4))) out$APOE4 <- as.integer(out$apoe_genotype) - 1

  out
}

extract_gene_from_prot <- function(colname) {
  if (!grepl("--", colname)) return(colname)
  strsplit(colname, "--", fixed = TRUE)[[1]][2]
}

# -----------------------------------------------------------------------------
# 1) Load PRS + merge into phenotypes
# -----------------------------------------------------------------------------
prs_path <- "~/Desktop/research/ROSMAP/ROSMAP_PRS_Leqi/AD_include_APOE_include_MHC_EUR_baseline_AMP_AD.txt"
prs_data <- read_prs_file(prs_path)

# Clear old columns if present (safe)
big_pheno_all <- big_pheno_all %>% dplyr::select(-any_of(c("PRScs_z", "SDPR_z")))
pheno_PFC     <- pheno_PFC %>% dplyr::select(-any_of(c("PRScs_z", "SDPR_z")))

big_pheno_all <- big_pheno_all %>% left_join(prs_data, by = "projid")
pheno_PFC     <- pheno_PFC %>% left_join(prs_data, by = "projid")

# Add consistent outcomes/codings
big_pheno_all <- ensure_binary_outcomes(big_pheno_all)
pheno_PFC     <- ensure_binary_outcomes(pheno_PFC)

# -----------------------------------------------------------------------------
# 2) PRS association with molecular features
# -----------------------------------------------------------------------------

# 2A) CpG association (cg06329447)
stopifnot("cg06329447" %in% colnames(DNAm_all))
cpg_df <- data.frame(
  projid = rownames(DNAm_all),
  CpG = as.numeric(DNAm_all[, "cg06329447"])
) %>%
  left_join(big_pheno_all, by = "projid")

model_cpg <- lm(CpG ~ PRScs_z + Female + Age + pmi + neurons_pfc, data = cpg_df)
print(summary(model_cpg))

# 2B) Protein association for a targeted set of genes
target_genes <- c(
  "GRM2","GSTK1","GLUD1","SYN3","LRRC47","SFXN1","CAP2","HPRT1",
  "PI4KA","DMXL2","ARFGEF3","KIF5C","OXR1","DPP3","VAMP1","FARSB",
  "CASKIN1","SARS2","GNAI2","HSPA8","VPS26B","VGF","SLC25A22",
  "NRXN3","ME3","GNAI3","NCOA7","GNAO1","AHNAK","FGG","HEBP1",
  "APEX1","RAB4A","SLC12A5","LRP1","BAG6","ARPC4","WDR1","EEF2",
  "ALDH7A1","ALDOC","EIF4B","SSBP1","NADK2","KBTBD11","IDH1","CRKL",
  "RYR2","GRIPAP1","BSN","APEH"
)

# Ensure proteins_PFC_scaled rownames are projid (recommended). If not, map from CleanedSampleID.
if (!all(rownames(proteins_PFC_scaled) %in% pheno_PFC$projid)) {
  if ("CleanedSampleID" %in% colnames(pheno_PFC)) {
    message("Mapping proteins_PFC_scaled rows via CleanedSampleID -> projid")
    map_df <- pheno_PFC %>% dplyr::select(projid, CleanedSampleID)
    idx <- match(rownames(proteins_PFC_scaled), map_df$CleanedSampleID)
    rownames(proteins_PFC_scaled) <- map_df$projid[idx]
  }
}

protein_genes <- sapply(colnames(proteins_PFC_scaled), extract_gene_from_prot)
selected_proteins <- colnames(proteins_PFC_scaled)[protein_genes %in% target_genes]

protein_prs_assoc <- bind_rows(lapply(selected_proteins, function(prot) {
  df <- data.frame(
    projid = rownames(proteins_PFC_scaled),
    Protein = as.numeric(proteins_PFC_scaled[, prot])
  ) %>%
    left_join(pheno_PFC %>%
                dplyr::select(projid, SDPR_z, PRScs_z, Female, Age, pmi, neurons_pfc),
              by = "projid") %>%
    filter(complete.cases(.))

  fit <- lm(Protein ~ SDPR_z + Female + Age + pmi + neurons_pfc, data = df)
  broom::tidy(fit) %>%
    mutate(ProteinName = prot)
}))

protein_prs_assoc_sdpr <- protein_prs_assoc %>%
  filter(term == "SDPR_z") %>%
  mutate(P_adj = p.adjust(p.value, method = "BH")) %>%
  arrange(P_adj)


# -----------------------------------------------------------------------------
# 3) PRS association with AD / neuropath outcomes
# -----------------------------------------------------------------------------
# Example using CERAD/Braak binaries 

# Ensure your ceradsc_num/braaksc_num exist or create them
if (!("braaksc_num" %in% names(big_pheno_all)) && "braaksc" %in% names(big_pheno_all)) {
  big_pheno_all$braaksc_num <- ifelse(big_pheno_all$braaksc %in% c(3,4), 1,
                                     ifelse(big_pheno_all$braaksc %in% c(1,2), 0, NA))
}
if (!("ceradsc_num" %in% names(big_pheno_all)) && "ceradsc" %in% names(big_pheno_all)) {
  big_pheno_all$ceradsc_num <- ifelse(big_pheno_all$ceradsc == 2, 1,
                                     ifelse(big_pheno_all$ceradsc == 1, 0, NA))
}

fit_ad <- glm(AD_bin ~ PRScs_z + Female + Age + pmi + neurons_pfc,
              data = big_pheno_all, family = binomial())
print(summary(fit_ad))

if ("braaksc_num" %in% names(big_pheno_all)) {
  fit_braak <- glm(braaksc_num ~ PRScs_z + Female + Age + pmi + neurons_pfc,
                   data = big_pheno_all, family = binomial())
  print(summary(fit_braak))
}

if ("ceradsc_num" %in% names(big_pheno_all)) {
  fit_cerad <- glm(ceradsc_num ~ PRScs_z + Female + Age + pmi + neurons_pfc,
                   data = big_pheno_all, family = binomial())
  print(summary(fit_cerad))
}

# -----------------------------------------------------------------------------
# 4) Prediction modeling comparison (Repeated CV) stratified by APOE4
# -----------------------------------------------------------------------------

# 4A) Build DNAme module eigengenes (subset + align to projid)
stopifnot(exists("datExpr"), exists("moduleColors_DNAme"))

# Keep only subjects present across pheno_PFC + datExpr + DNAm_all
common_projid <- Reduce(intersect, list(
  pheno_PFC$projid,
  rownames(DNAm_all),
  rownames(datExpr)
))

# CpG data aligned
cpg_common <- data.frame(
  projid = common_projid,
  CpG = as.numeric(DNAm_all[common_projid, "cg06329447"])
)

# DNAme MEs aligned
MEs0 <- moduleEigengenes(datExpr[common_projid, , drop = FALSE],
                         colors = moduleColors_DNAme)$eigengenes
MEs  <- orderMEs(MEs0)

# Keep the specific modules of interest
keep_MEs <- intersect(colnames(MEs), c("MElightgreen","MEsalmon","MEmidnightblue","MEpink"))
MEs_keep <- MEs[, keep_MEs, drop = FALSE]
MEs_df <- as.data.frame(MEs_keep) %>% rownames_to_column("projid")

# 4B) Proteomic PCs aligned to projid
prot_common <- proteins_PFC_scaled[common_projid, , drop = FALSE]

groups_list <- list(
  Cell_activation_AND_Chemical_synaptic_transmission = c("LRP1","VGF","FGG","BAG6","HEBP1","SLC12A5","BSN"),
  Negative_regulation_of_adenylate_cyclase = c("GNAO1","GRM2","GNAI2","CASKIN1"),
  MET_receptor_recycling_AND_Lipid_metabolism = c("PI4KA","CRKL","RAB4A","GRIPAP1"),
  TCA_cycle = c("GLUD1","DLST","ME3"),
  Protein_translation = c("FARSB","SARS2"),
  Vesicle_mediated_transport = c("SYN3","VAMP1"),
  Protein_folding = c("HSPA8","DPP3")
)

prot_genes <- sapply(colnames(prot_common), extract_gene_from_prot)

protein_PC1_df <- data.frame(projid = common_projid)
for (g in names(groups_list)) {
  genes <- groups_list[[g]]
  cols <- which(prot_genes %in% genes)
  if (length(cols) < 2) {
    message("Skipping proteomics group ", g, " (<2 proteins found).")
    next
  }
  X <- prot_common[, cols, drop = FALSE]
  pca <- prcomp(X, center = TRUE, scale. = TRUE)
  protein_PC1_df[[g]] <- pca$x[,1]
}

# 4C) PRS aligned
prs_common <- prs_data %>% filter(projid %in% common_projid)

# 4D) Merge into pred_data
pred_data <- pheno_PFC %>%
  filter(projid %in% common_projid) %>%
  dplyr::select(projid, Female, Age, pmi, neurons_pfc, AD_bin, AD_status, APOE4) %>%
  left_join(cpg_common, by = "projid") %>%
  left_join(MEs_df, by = "projid") %>%
  left_join(protein_PC1_df, by = "projid") %>%
  left_join(prs_common, by = "projid")

pred_data$AD_status <- factor(pred_data$AD_status, levels = c("Control", "AD"))

# -----------------------------------------------------------------------------
# 4E) Define model formulas
# -----------------------------------------------------------------------------
covariates <- c("Female","Age","pmi","neurons_pfc")

proteomics_vars <- setdiff(names(protein_PC1_df), "projid")
me_vars <- keep_MEs

formula_list <- list(
  PRS = as.formula(paste("AD_status ~ PRScs_z +", paste(covariates, collapse=" + "))),
  DNAme_CpG = as.formula(paste("AD_status ~ CpG +", paste(covariates, collapse=" + "))),
  DNAme_MEs = if (length(me_vars) > 0)
    as.formula(paste("AD_status ~", paste(me_vars, collapse=" + "), "+", paste(covariates, collapse=" + ")))
  else NULL,
  Proteomics = if (length(proteomics_vars) > 0)
    as.formula(paste("AD_status ~", paste(proteomics_vars, collapse=" + "), "+", paste(covariates, collapse=" + ")))
  else NULL
)

formula_list <- formula_list[!sapply(formula_list, is.null)]

# -----------------------------------------------------------------------------
# 4F) Fair comparison: complete-cases across ALL predictors used in ANY model
# -----------------------------------------------------------------------------
all_predictors <- unique(c(
  "AD_status","APOE4",
  covariates,
  "PRScs_z","CpG",
  me_vars,
  proteomics_vars
))
all_predictors <- intersect(all_predictors, names(pred_data))

pred_data_cc <- pred_data %>%
  dplyr::select(all_of(all_predictors)) %>%
  filter(complete.cases(.)) %>%
  bind_cols(pred_data %>% dplyr::select(projid), .)

message("N complete-case for all models: ", nrow(pred_data_cc))

# -----------------------------------------------------------------------------
# 4G) Train repeated CV models + get out-of-fold probs
# -----------------------------------------------------------------------------
ctrl <- trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final"
)

set.seed(1)
models_cv <- lapply(formula_list, function(f) {
  train(f, data = pred_data_cc, method = "glm", family = binomial(),
        metric = "ROC", trControl = ctrl)
})

# Gather out-of-fold predictions
preds <- bind_rows(lapply(names(models_cv), function(mn) {
  m <- models_cv[[mn]]
  df <- m$pred

  # Identify probability column for "AD" robustly
  prob_col <- "AD"
  if (!(prob_col %in% names(df))) {
    # fallback: take last factor level probability column if caret named it differently
    levs <- levels(pred_data_cc$AD_status)
    prob_col <- tail(levs, 1)
  }

  df %>%
    mutate(
      Model = mn,
      APOE4 = pred_data_cc$APOE4[rowIndex],
      ProbAD = .data[[prob_col]],
      TrueAD = obs
    ) %>%
    dplyr::select(rowIndex, Model, APOE4, TrueAD, ProbAD)
}))

# ROC curves
roc_dfs <- preds %>%
  group_by(Model, APOE4) %>%
  nest() %>%
  mutate(
    roc_obj = map(data, ~ roc(.x$TrueAD, .x$ProbAD, quiet = TRUE, levels = c("Control","AD"))),
    roc_df = map(roc_obj, ~ tibble(
      FPR = 1 - .$specificities,
      TPR = .$sensitivities
    ))
  ) %>%
  dplyr::select(Model, APOE4, roc_df) %>%
  unnest(roc_df)

auc_tbl <- preds %>%
  group_by(Model, APOE4) %>%
  summarize(
    AUC = if (n_distinct(TrueAD) > 1) as.numeric(auc(roc(TrueAD, ProbAD, levels=c("Control","AD")))) else NA_real_,
    .groups = "drop"
  ) %>%
  group_by(APOE4) %>%
  mutate(y_pos = seq(0.95, 0.75, length.out = n())) %>%
  ungroup()

# Plot ROC
p <- ggplot(roc_dfs, aes(FPR, TPR, color = Model)) +
  geom_line(linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  facet_wrap(~ APOE4, ncol = 2,
             labeller = as_labeller(c(`0`="Non-carriers", `1`="Carriers"))) +
  geom_text(
    data = auc_tbl,
    aes(x = 0.2, y = y_pos, label = paste0(Model, ": AUC=", round(AUC, 2)), color = Model),
    inherit.aes = FALSE, size = 4, hjust = 0
  ) +
  labs(x = "False Positive Rate", y = "True Positive Rate",
       title = "ROC Curves (out-of-fold) by APOE4 Status") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

print(p)
