# =============================================================================
# 12_proteomics_ruv_normalization_and_differential.R
#
# Inputs expected in environment:
#   - peptide_intensities_df_grouped_norm  (sample x gene/protein intensities + SampleID/Batch/Region)
#   - temp_rt_calibration_intensities_df   (sample x RT peptide intensities; must include Sample + control_peptide_columns)
#   - control_peptide_columns              (character vector of RTC intensity column names)
#   - pheno_all                            (phenotype table incl SampleID, Region, niareagansc/niareagansc_recode, Female, Age, apoe_long, neurons_pfc, pmi, apoe_genotype)
#
# Outputs:
#   - proteins_PFC, pheno_PFC  (RUV-adjusted protein matrix + phenotype, aligned)
# =============================================================================

suppressPackageStartupMessages({
  library(ruv)
  library(limma)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
})


# -----------------------------------------------------------------------------
# 1) Helper: robustly parse protein column names to "ACCESSION--GENE"
# -----------------------------------------------------------------------------
make_feature_name <- function(x) {
  # Expecting UniProt-ish strings with |ACC| and GN=
  acc  <- str_extract(x, "(?<=\\|)[A-Z0-9]+(?=\\|)")
  gene <- str_extract(x, "GN=([^ ]+)")
  gene <- ifelse(is.na(gene), NA_character_, sub("^GN=", "", gene))

  # If accession missing, fall back to original name (shortened)
  acc <- ifelse(is.na(acc), x, acc)

  ifelse(is.na(gene), acc, paste0(acc, "--", gene))
}

# -----------------------------------------------------------------------------
# 2) Build protein matrix (samples x features), keep SampleID rownames
# -----------------------------------------------------------------------------
normalized_columns <- setdiff(
  names(peptide_intensities_df_grouped_norm),
  c("SampleID", "Batch", "Region")
)

Y_proteins <- peptide_intensities_df_grouped_norm %>%
  select(SampleID, all_of(normalized_columns)) %>%
  column_to_rownames("SampleID")

# Coerce to numeric matrix
Y_proteins <- as.matrix(sapply(as.data.frame(Y_proteins), as.numeric))
rownames(Y_proteins) <- rownames(peptide_intensities_df_grouped_norm)

sample_ids <- rownames(Y_proteins)

# Replace zeros with NA (consistent with earlier)
Y_proteins[Y_proteins == 0] <- NA

# Rename features (columns)
colnames(Y_proteins) <- vapply(colnames(Y_proteins), make_feature_name, character(1))

# Collapse duplicated feature names by summing 
if (any(duplicated(colnames(Y_proteins)))) {
  Y_proteins_df <- as.data.frame(Y_proteins) %>%
    rownames_to_column("SampleID") %>%
    pivot_longer(-SampleID, names_to = "feature", values_to = "value") %>%
    group_by(SampleID, feature) %>%
    summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = feature, values_from = value)

  Y_proteins <- Y_proteins_df %>%
    column_to_rownames("SampleID") %>%
    as.matrix()
}

# -----------------------------------------------------------------------------
# 3) Build control matrix from RTC table and merge samples
# -----------------------------------------------------------------------------
# Controls should be RTC intensity columns only; ensure they exist:
control_peptide_columns <- intersect(control_peptide_columns, colnames(temp_rt_calibration_intensities_df))
stopifnot("Sample" %in% colnames(temp_rt_calibration_intensities_df))
stopifnot(length(control_peptide_columns) > 0)

Y_controls <- temp_rt_calibration_intensities_df %>%
  select(Sample, all_of(control_peptide_columns)) %>%
  mutate(across(all_of(control_peptide_columns), ~ as.numeric(as.character(.)))) %>%
  column_to_rownames("Sample") %>%
  as.matrix()

# Intersect samples present in BOTH
common_samples <- intersect(rownames(Y_proteins), rownames(Y_controls))

Y_proteins <- Y_proteins[common_samples, , drop = FALSE]
Y_controls <- Y_controls[common_samples, , drop = FALSE]

# -----------------------------------------------------------------------------
# 4) Restrict to PFC samples with non-missing trait info
# -----------------------------------------------------------------------------
present_samples <- pheno_all$SampleID[
  !is.na(pheno_all$niareagansc) &
    pheno_all$Region == "PFC" &
    pheno_all$SampleID %in% common_samples
]

Y_proteins <- Y_proteins[present_samples, , drop = FALSE]
Y_controls <- Y_controls[present_samples, , drop = FALSE]

# -----------------------------------------------------------------------------
# 5) Missingness filter on PROTEINS only (controls excluded)
# -----------------------------------------------------------------------------
missing_percent <- colSums(is.na(Y_proteins)) / nrow(Y_proteins)

# adjust the threshold
threshold_missing <- 0.99
keep_features <- names(missing_percent)[missing_percent < threshold_missing]
Y_proteins <- Y_proteins[, keep_features, drop = FALSE]
missing_percent <- missing_percent[keep_features]

message(sprintf("Retained %d / %d protein features after missingness filter.",
                length(keep_features), length(missing_percent)))

# -----------------------------------------------------------------------------
# 6) Impute (choose the method)
# -----------------------------------------------------------------------------
impute_min <- function(x) {
  x[is.na(x)] <- min(x, na.rm = TRUE)
  x
}
# Alternative:
# impute_median <- function(x) { x[is.na(x)] <- median(x, na.rm = TRUE); x }

Y_proteins_imp  <- apply(Y_proteins,  2, impute_min)
Y_controls_imp  <- apply(Y_controls,  2, impute_min)

Y_proteins_imp <- as.matrix(apply(as.data.frame(Y_proteins_imp), 2, as.numeric))
Y_controls_imp <- as.matrix(apply(as.data.frame(Y_controls_imp), 2, as.numeric))

rownames(Y_proteins_imp) <- rownames(Y_proteins)
rownames(Y_controls_imp) <- rownames(Y_controls)

# -----------------------------------------------------------------------------
# 7) Build combined matrix for RUV (samples x features_total)
# -----------------------------------------------------------------------------
Y_combined <- cbind(Y_proteins_imp, Y_controls_imp)

# Define control features for RUV: exactly the RTC columns
ctl <- colnames(Y_combined) %in% colnames(Y_controls_imp)

# -----------------------------------------------------------------------------
# 8) Align phenotype subset to samples
# -----------------------------------------------------------------------------
pheno_subset <- pheno_all %>%
  filter(SampleID %in% rownames(Y_combined)) %>%
  arrange(match(SampleID, rownames(Y_combined))) %>%
  column_to_rownames("SampleID")

stopifnot(identical(rownames(pheno_subset), rownames(Y_combined)))

# filter APOE genotypes 
target_genotypes <- c("33", "34", "44")
pheno_subset <- pheno_subset %>% filter(apoe_long %in% target_genotypes)

# Drop outlier samples if desired 
outliers <- c("146P", "90P_r")
pheno_subset <- pheno_subset[!(rownames(pheno_subset) %in% outliers), , drop = FALSE]
Y_combined   <- Y_combined[rownames(pheno_subset), , drop = FALSE]

# -----------------------------------------------------------------------------
# 9) RUV4
# -----------------------------------------------------------------------------
# Design for RUV: factor of interest + covariates; adjust the covariates
X_ruv <- model.matrix(~ niareagansc + Female + Age + apoe_long, data = pheno_subset)

Y_ruv <- t(log2(Y_combined + 1))

# Choose k 
k <- 13 # eg max is the number of negative controls; adjust
ruv4_results <- RUV4(Y = Y_ruv, X = X_ruv, ctl = ctl, k = k, Z = NULL)

# Adjusted Y (features x samples)
W <- ruv4_results$W
alpha <- ruv4_results$alpha
unwanted_variation <- W %*% alpha
Y_adjusted_ruv <- Y_ruv - unwanted_variation

# Convert back to samples x features for PCA, etc.
Y_adjusted <- t(Y_adjusted_ruv)
rownames(Y_adjusted) <- rownames(Y_combined)
colnames(Y_adjusted) <- colnames(Y_combined)

# Keep only protein features as "proteins_PFC" output at end
protein_feature_names <- colnames(Y_proteins_imp)
proteins_PFC_adjusted <- Y_adjusted[, protein_feature_names, drop = FALSE]

# -----------------------------------------------------------------------------
# 10) PCA quick QC plots (optional)
# -----------------------------------------------------------------------------
pca_result <- prcomp(proteins_PFC_adjusted, scale. = TRUE)
pca_data <- data.frame(
  Sample = rownames(proteins_PFC_adjusted),
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2]
)

p1 <- ggplot(pca_data, aes(PC1, PC2, color = pheno_subset$apoe_long)) +
  geom_point(alpha = 0.8, size = 3) + theme_minimal() + labs(color = "APOE")

p2 <- ggplot(pca_data, aes(PC1, PC2, color = pheno_subset$niareagansc)) +
  geom_point(alpha = 0.8, size = 3) + theme_minimal() + labs(color = "NIA AD")

p3 <- ggplot(pca_data, aes(PC1, PC2, color = pheno_subset$Female)) +
  geom_point(alpha = 0.8, size = 3) + theme_minimal() + labs(color = "Sex")

p4 <- ggplot(pca_data, aes(PC1, PC2, color = pheno_subset$Age)) +
  geom_point(alpha = 0.8, size = 3) + theme_minimal() + labs(color = "Age")

print(p1); print(p2); print(p3); print(p4)

# -----------------------------------------------------------------------------
# 11) Differential testing with limma
# -----------------------------------------------------------------------------
run_limma <- function(Y_samples_by_features, pheno_df, design_formula, coef_name) {
  # Y_samples_by_features: samples x features
  Y_fxs <- t(Y_samples_by_features)  # features x samples
  design <- model.matrix(design_formula, data = pheno_df)

  fit <- lmFit(Y_fxs, design)
  fit <- eBayes(fit)

  tt <- topTable(fit, coef = coef_name, number = Inf, adjust.method = "BH")
  tt$Protein <- rownames(tt)
  tt
}

# ---- A) NIA Reagan AD effect (uses niareagansc_recode + covariates)
pheno_subset$pmi[is.na(pheno_subset$pmi)] <- mean(pheno_subset$pmi, na.rm = TRUE)

protein_results_reagan <- run_limma(
  Y_samples_by_features = proteins_PFC_adjusted[, protein_feature_names, drop = FALSE],
  pheno_df = pheno_subset,
  design_formula = ~ niareagansc_recode + Female + Age + neurons_pfc + pmi,
  coef_name = "niareagansc_recode1"
)

# Build volcano table + export
volcano_reagan <- protein_results_reagan %>%
  transmute(
    Protein,
    logFC = logFC,
    P.Value = P.Value,
    adj.P.Val = adj.P.Val,
    minusLog10P = -log10(P.Value)
  )
volcano_reagan$missing <- missing_percent[volcano_reagan$Protein]
volcano_reagan$gene <- sub(".*--(.*)", "\\1", volcano_reagan$Protein)


# ---- B) APOE effect (apoe_genotype TRUE/FALSE + covariates)
protein_results_apoe <- run_limma(
  Y_samples_by_features = proteins_PFC_adjusted[, protein_feature_names, drop = FALSE],
  pheno_df = pheno_subset,
  design_formula = ~ apoe_genotype + Female + Age + neurons_pfc + pmi,
  coef_name = "apoe_genotypeTRUE"
)

volcano_apoe <- protein_results_apoe %>%
  transmute(
    Protein,
    logFC = logFC,
    P.Value = P.Value,
    adj.P.Val = adj.P.Val,
    minusLog10P = -log10(P.Value)
  )
volcano_apoe$missing <- missing_percent[volcano_apoe$Protein]
volcano_apoe$gene <- sub(".*--(.*)", "\\1", volcano_apoe$Protein)


# -----------------------------------------------------------------------------
# 12) Volcano plot helper (optional)
# -----------------------------------------------------------------------------
plot_volcano <- function(volcano_df, title, logFC_threshold = 0.58, fdr_line = 0.10) {
  # Compute y-threshold corresponding to FDR line (approx via smallest p among those below FDR)
  idx <- which(volcano_df$adj.P.Val < fdr_line)
  y_fdr <- if (length(idx) == 0) NA_real_ else -log10(min(volcano_df$P.Value[idx], na.rm = TRUE))

  volcano_df$Significant <- with(volcano_df,
                                 ifelse(adj.P.Val <= 0.05 & logFC >= logFC_threshold, "Upregulated",
                                        ifelse(adj.P.Val <= 0.05 & logFC <= -logFC_threshold, "Downregulated",
                                               "Not Significant")))

  p <- ggplot(volcano_df, aes(x = logFC, y = minusLog10P)) +
    geom_point(alpha = 0.4) +
    geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed") +
    theme_minimal() +
    labs(title = title, x = "log2 Fold Change", y = "-log10(P)")

  if (!is.na(y_fdr)) {
    p <- p + geom_hline(yintercept = y_fdr, linetype = "dashed")
  }

  p + geom_text_repel(
    data = subset(volcano_df, adj.P.Val <= 0.05 & abs(logFC) >= logFC_threshold),
    aes(label = Protein),
    size = 2,
    max.overlaps = 30
  )
}

print(plot_volcano(volcano_reagan, "PFC: Reagan AD vs control (RUV4-adjusted)"))
print(plot_volcano(volcano_apoe, "PFC: APOE genotype effect (RUV4-adjusted)"))

# -----------------------------------------------------------------------------
# 13) Final objects for downstream scripts
# -----------------------------------------------------------------------------
proteins_PFC <- proteins_PFC_adjusted
pheno_PFC <- pheno_subset
