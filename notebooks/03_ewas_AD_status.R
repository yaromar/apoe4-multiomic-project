# Yaro
# Last update: 12/16/23
#
# Script Purpose:
#   Identify differentially methylated positions (DMPs) associated with AD pathology
#   using limma on EPIC array beta values (probe-wise differential methylation).
#   Includes:
#     1) Preprocessing (sample alignment, filtering sex chr, QC filtering, variance filter)
#     2) EWAS via limma (AD vs control, adjusted covariates)
#     3) Manhattan + volcano plots
#     4) Optional: APOE+ only EWAS
#     5) Optional: CpG-level APOE×Methylation interaction scan via logistic regression
#     6) Optional: GO/KEGG enrichment and DMR analysis (DMRcate)
#

# ---- Libraries ----
library(limma)
library(minfi)
library(WGCNA)

library(dplyr)
library(stringr)

library(qqman)
library(ggplot2)
library(ggrepel)

library(missMethyl)
library(DMRcate)
library(GenomicRanges)

set.seed(0)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

# ---- Annotation ----
annEPIC <- getAnnotation("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")

# ---- Inputs ----
# Use all available phenotype data to maximize power
big_pheno_all <- big_pheno_scaled

# ---- Align phenotype to PFC methylation data ----
# Subset to samples that exist in both pheno and methylation
DNAm_all <- PFCSamples[rownames(PFCSamples) %in% big_pheno_all$projid, , drop = FALSE]
DNAm_all <- DNAm_all[match(big_pheno_all$projid, rownames(DNAm_all)), , drop = FALSE]

# Keep only phenotypes that actually have methylation data after matching
big_pheno_all <- big_pheno_all[big_pheno_all$projid %in% rownames(DNAm_all), , drop = FALSE]
DNAm_all <- DNAm_all[match(big_pheno_all$projid, rownames(DNAm_all)), , drop = FALSE]

stopifnot(identical(rownames(DNAm_all), as.character(big_pheno_all$projid)))

# ---- Remove sex chromosome CpGs ----
# annEPIC$Name contains CpG ids; annEPIC rownames are also usually CpG ids, but safer to use $Name.
sex_cpgs <- annEPIC$Name[annEPIC$chr %in% c("chrX", "chrY")]
DNAm_all <- DNAm_all[, !(colnames(DNAm_all) %in% sex_cpgs), drop = FALSE]

# ---- Remove bad CpGs/samples (WGCNA helper) ----
gsg <- goodSamplesGenes(DNAm_all, verbose = 3)
if (!gsg$allOK) {
  if (sum(!gsg$goodGenes) > 0) {
    printFlush(paste("Removing genes:", paste(colnames(DNAm_all)[!gsg$goodGenes], collapse = ", ")))
  }
  if (sum(!gsg$goodSamples) > 0) {
    printFlush(paste("Removing samples:", paste(rownames(DNAm_all)[!gsg$goodSamples], collapse = ", ")))
  }
  DNAm_all <- DNAm_all[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
  big_pheno_all <- big_pheno_all[big_pheno_all$projid %in% rownames(DNAm_all), , drop = FALSE]
  DNAm_all <- DNAm_all[match(big_pheno_all$projid, rownames(DNAm_all)), , drop = FALSE]
}
stopifnot(identical(rownames(DNAm_all), as.character(big_pheno_all$projid)))

# ---- OPTIONAL: Outlier filtering separately within AD and control, then recombine ----
do_outlier_filter <- TRUE

if (do_outlier_filter) {

  DNAm_allAD <- DNAm_all[big_pheno_all$niareagansc_recode == 1, , drop = FALSE]
  DNAm_allControl <- DNAm_all[big_pheno_all$niareagansc_recode == 0, , drop = FALSE]

  # --- AD cluster filtering ---
  sampleTreeAD <- hclust(dist(DNAm_allAD), method = "average")
  clustAD <- cutreeStatic(sampleTreeAD, cutHeight = 85, minSize = 1)
  keepAD <- (clustAD == 1)
  DNAm_allAD <- DNAm_allAD[keepAD, , drop = FALSE]

  # --- Control cluster filtering ---
  sampleTreeCtrl <- hclust(dist(DNAm_allControl), method = "average")
  clustCtrl <- cutreeStatic(sampleTreeCtrl, cutHeight = 55, minSize = 1)
  keepCtrl <- (clustCtrl %in% c(1, 2))
  DNAm_allControl <- DNAm_allControl[keepCtrl, , drop = FALSE]

  # Recombine
  DNAm_all <- rbind(DNAm_allAD, DNAm_allControl)

  # Re-align phenotype
  big_pheno_all <- big_pheno_all[big_pheno_all$projid %in% rownames(DNAm_all), , drop = FALSE]
  DNAm_all <- DNAm_all[match(big_pheno_all$projid, rownames(DNAm_all)), , drop = FALSE]
  stopifnot(identical(rownames(DNAm_all), as.character(big_pheno_all$projid)))
}

# ---- Variance filter: remove lowest 10% SD CpGs ----
sd_cpgs <- apply(DNAm_all, 2, sd, na.rm = TRUE)
variance_threshold <- quantile(sd_cpgs, 0.1, na.rm = TRUE)
cpgs_to_keep <- which(sd_cpgs > variance_threshold)
DNAm_all <- DNAm_all[, cpgs_to_keep, drop = FALSE]

# ---- Restrict to APOE long genotypes 33/34/44 ----
target_genotypes <- c("33", "34", "44")
big_pheno_all <- big_pheno_all %>% filter(apoe_long %in% target_genotypes)
DNAm_all <- DNAm_all[match(big_pheno_all$projid, rownames(DNAm_all)), , drop = FALSE]
stopifnot(identical(rownames(DNAm_all), as.character(big_pheno_all$projid)))

# ---- Handle missing covariates ----
big_pheno_all$pmi[is.na(big_pheno_all$pmi)] <- mean(big_pheno_all$pmi, na.rm = TRUE)

# ---- EWAS with limma ----
# Model: methylation ~ AD + age + sex + neurons + PMI
# niareagansc_recode is assumed to be coded 0/1; we rely on model.matrix to generate coef name.
designAll <- model.matrix(~ niareagansc_recode + Age + Female + neurons_pfc + pmi, data = big_pheno_all)

# limma wants features in rows, samples in columns => transpose
fit <- lmFit(t(DNAm_all), designAll)
fit <- eBayes(fit)

# Attach annotation columns to topTable output (using CpG IDs in rownames(fit))
anno_subset <- annEPIC[match(rownames(fit), annEPIC$Name), c(1:4, 12:19, 24:ncol(annEPIC))]

results_niareagansc_recode_all <- topTable(
  fit,
  num = Inf,
  coef = "niareagansc_recode1",
  genelist = anno_subset
)

# Make CpG ID explicit as rownames
rownames(results_niareagansc_recode_all) <- results_niareagansc_recode_all$Name

# Quick check for “near-significant” set
results_niareagansc_recode_all[results_niareagansc_recode_all$adj.P.Val <= 0.11, ]

# ---- Gene label cleanup for plotting (only for “highlighted” CpGs) ----
sig_for_labels <- results_niareagansc_recode_all %>% filter(adj.P.Val <= 0.11)

sig_for_labels <- sig_for_labels %>%
  mutate(Combined_Names = paste(GencodeBasicV12_NAME, GencodeCompV12_NAME, sep = ";")) %>%
  mutate(Unique_Genes = sapply(
    str_split(Combined_Names, ";"),
    function(x) {
      cleaned_names <- str_replace_all(x, "\\bRP\\d+-.+", "")
      unique_names <- unique(cleaned_names)
      final_names <- paste(unique_names, collapse = ";")
      str_remove_all(final_names, "^;+")
    }
  ))

results_niareagansc_recode_all <- left_join(
  results_niareagansc_recode_all,
  sig_for_labels[, c("Name", "Unique_Genes")],
  by = "Name"
)

# ---- Manhattan plot prep ----
results_niareagansc_recode_all$chr_numeric <- sub("chr", "", results_niareagansc_recode_all$chr)
results_niareagansc_recode_all <- results_niareagansc_recode_all[results_niareagansc_recode_all$chr != "chrY", ]
results_niareagansc_recode_all$chr_numeric <- suppressWarnings(as.numeric(results_niareagansc_recode_all$chr_numeric))

# ---- Manhattan plot ----
pdf("Manhattan_Plot_PFC.pdf", width = 10, height = 6)
manhattan(
  results_niareagansc_recode_all,
  chr = "chr_numeric",
  bp = "pos",
  p = "P.Value",
  snp = "Unique_Genes",
  suggestiveline = FALSE,
  genomewideline = -log10(3.61e-6),
  annotatePval = 3.61e-6,
  annotateTop = FALSE,
  main = "Differential Methylation Manhattan Plot (PFC)"
)
abline(h = -log10(4.86e-7), col = "red", lty = 3)
dev.off()

# ---- Volcano plot (single version) ----
volcano_plot <- ggplot(results_niareagansc_recode_all, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(alpha = 0.4, color = "gray40") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(3.61e-6), linetype = "dashed", color = "red", linewidth = 0.5) +
  geom_hline(yintercept = -log10(4.86e-7), linetype = "dotted", color = "red", linewidth = 0.5) +
  geom_text_repel(
    aes(label = ifelse(P.Value <= 3.61e-6, Unique_Genes, "")),
    segment.color = "grey50",
    size = 2,
    max.overlaps = 30
  ) +
  labs(
    x = expression(Delta~beta~"(AD - Control)"),
    y = expression(-log[10](italic(P))),
    title = "Volcano Plot of Differential Methylation (PFC)"
  ) +
  theme_minimal(base_size = 12)

ggsave("Volcano_Plot_PFC.pdf", plot = volcano_plot, width = 10, height = 6, units = "in", dpi = 300)

# ---- Export CpGs (example) ----
write.csv(
  data.frame(CpG = rownames(results_niareagansc_recode_all[results_niareagansc_recode_all$P.Value < 0.05, ])),
  "niareagansc_recode_cpgs_PFC.csv",
  row.names = FALSE
)

# ---- OPTIONAL: APOE+ only EWAS ----
do_apoe_positive_only <- TRUE

if (do_apoe_positive_only) {

  big_pheno_APOE <- big_pheno_scaled[big_pheno_scaled$apoe_genotype == TRUE, , drop = FALSE]

  DNAm_APOE <- PFCSamples[rownames(PFCSamples) %in% big_pheno_APOE$projid, , drop = FALSE]
  DNAm_APOE <- DNAm_APOE[match(big_pheno_APOE$projid, rownames(DNAm_APOE)), , drop = FALSE]

  # Keep only samples retained in main filtered set (post-QC + genotype filtering)
  big_pheno_APOE <- big_pheno_APOE[big_pheno_APOE$projid %in% big_pheno_all$projid, , drop = FALSE]
  DNAm_APOE <- DNAm_APOE[match(big_pheno_APOE$projid, rownames(DNAm_APOE)), , drop = FALSE]

  # Reuse sex chr removal and variance filter for this subset
  DNAm_APOE <- DNAm_APOE[, !(colnames(DNAm_APOE) %in% sex_cpgs), drop = FALSE]
  sd_cpgs_apoe <- apply(DNAm_APOE, 2, sd, na.rm = TRUE)
  var_thresh_apoe <- quantile(sd_cpgs_apoe, 0.1, na.rm = TRUE)
  DNAm_APOE <- DNAm_APOE[, which(sd_cpgs_apoe > var_thresh_apoe), drop = FALSE]

  big_pheno_APOE$pmi[is.na(big_pheno_APOE$pmi)] <- mean(big_pheno_APOE$pmi, na.rm = TRUE)

  designAPOE <- model.matrix(~ niareagansc_recode + Age + Female + neurons_pfc + pmi, data = big_pheno_APOE)

  fit_apoe <- lmFit(t(DNAm_APOE), designAPOE)
  fit_apoe <- eBayes(fit_apoe)

  anno_subset_apoe <- annEPIC[match(rownames(fit_apoe), annEPIC$Name), c(1:4, 12:19, 24:ncol(annEPIC))]
  results_niareagansc_recode_APOE <- topTable(
    fit_apoe,
    num = Inf,
    coef = "niareagansc_recode1",
    genelist = anno_subset_apoe
  )
  rownames(results_niareagansc_recode_APOE) <- results_niareagansc_recode_APOE$Name

  # Example: subset for feature selection
  results_niareagansc_recode_APOE <- results_niareagansc_recode_APOE[results_niareagansc_recode_APOE$adj.P.Val < 0.06, ]
  cpg_niareagansc_recode_APOE <- rownames(results_niareagansc_recode_APOE)
}

# ---- OPTIONAL: APOE × methylation interaction scan across CpGs ----
# This is independent of limma EWAS and can be slow.
# It scans: AD_status ~ APOE_status * CpG + covariates
do_interaction_scan <- TRUE

if (do_interaction_scan) {

  # Make sure these are factors as expected by glm
  big_pheno_all$AD_status <- as.factor(big_pheno_all$niareagansc_recode)
  big_pheno_all$APOE_status <- as.factor(big_pheno_all$apoe_genotype)

  big_pheno_all$pmi[is.na(big_pheno_all$pmi)] <- mean(big_pheno_all$pmi, na.rm = TRUE)

  perform_analysis <- function(cpg) {

    molecule_values <- DNAm_all[, cpg]

    model_data <- data.frame(
      AD_status = big_pheno_all$AD_status,
      APOE_status = big_pheno_all$APOE_status,
      Molecule = molecule_values,
      Female = big_pheno_all$Female,
      Age = big_pheno_all$Age,
      pmi = big_pheno_all$pmi,
      Neurons = big_pheno_all$neurons_pfc
    )

    model <- glm(
      AD_status ~ APOE_status * Molecule + Female + Age + pmi + Neurons,
      data = model_data,
      family = binomial(link = "logit")
    )

    coefs <- summary(model)$coefficients

    # Defensive: if factor encoding differs, these rows may not exist
    get_row <- function(rn) {
      if (rn %in% rownames(coefs)) coefs[rn, ] else c(NA, NA, NA, NA)
    }

    apoe_row <- get_row("APOE_statusTRUE")
    mol_row <- get_row("Molecule")
    int_row <- get_row("APOE_statusTRUE:Molecule")

    data.frame(
      CpG = cpg,
      Beta_APOE = apoe_row[1],
      SE_APOE = apoe_row[2],
      P_APOE = apoe_row[4],
      Beta_Molecule = mol_row[1],
      SE_Molecule = mol_row[2],
      P_Molecule = mol_row[4],
      Beta_Interaction = int_row[1],
      SE_Interaction = int_row[2],
      P_Interaction = int_row[4],
      stringsAsFactors = FALSE
    )
  }

  # Parallelize with mclapply (works on mac/linux)
  numCores <- max(1, parallel::detectCores() - 1)

  methylation_interaction_results <- do.call(
    rbind,
    parallel::mclapply(colnames(DNAm_all), perform_analysis, mc.cores = numCores)
  )

  methylation_interaction_results$Adj_P_Interaction <- p.adjust(
    methylation_interaction_results$P_Interaction,
    method = "BH"
  )
}

# ---- OPTIONAL: Merge EWAS and interaction results ----
do_merge_interaction <- TRUE

if (do_merge_interaction && exists("methylation_interaction_results")) {

  combined_results <- results_niareagansc_recode_all %>%
    mutate(Name = Name) %>%
    full_join(
      methylation_interaction_results %>%
        transmute(Name = CpG, APOE_Interaction_Significance = P_Interaction),
      by = "Name"
    )

  combined_results <- combined_results %>%
    mutate(Combined_Names = paste(GencodeBasicV12_NAME, GencodeCompV12_NAME, sep = ";")) %>%
    mutate(Unique_Genes = sapply(
      str_split(Combined_Names, ";"),
      function(x) {
        cleaned_names <- str_replace_all(x, "\\bRP\\d+-.+", "")
        unique_names <- unique(cleaned_names)
        final_names <- paste(unique_names, collapse = ";")
        str_remove_all(final_names, "^;+")
      }
    ))

  write.csv(
    combined_results,
    "niareagansc_recode_cpgs_PFC_significant.csv",
    row.names = FALSE
  )
}

# ---- OPTIONAL: Enrichment (CpG-level) ----
do_cpg_enrichment <- TRUE

if (do_cpg_enrichment) {

  sig_cpgs <- results_niareagansc_recode_all$Name[results_niareagansc_recode_all$adj.P.Val < 0.06]

  if (length(sig_cpgs) > 0) {

    gst_go <- gometh(sig.cpg = sig_cpgs, all.cpg = colnames(DNAm_all), array.type = "EPIC")
    topGSA(gst_go[gst_go$P.DE < 0.05, ])

    gst_kegg <- gometh(sig.cpg = sig_cpgs, all.cpg = colnames(DNAm_all), array.type = "EPIC", collection = "KEGG")
    topGSA(gst_kegg[gst_kegg$P.DE < 0.05, ])
  }
}

# ---- OPTIONAL: DMR analysis via DMRcate ----
do_dmr <- TRUE

if (do_dmr) {

  myAnnotationAll <- cpg.annotate(
    object = t(DNAm_all),
    datatype = "array",
    what = "Beta",
    analysis.type = "differential",
    design = designAll,
    contrasts = FALSE,
    coef = "niareagansc_recode1",
    arraytype = "EPIC"
  )

  DMRsAll <- dmrcate(myAnnotationAll)
  results.rangesAll <- extractRanges(DMRsAll)

  # Region enrichment (GO/KEGG)
  gst_reg_go <- goregion(regions = results.rangesAll, all.cpg = colnames(DNAm_all), array.type = "EPIC")
  topGSA(gst_reg_go[gst_reg_go$P.DE < 0.05 & gst_reg_go$ONTOLOGY == "BP", ], num = Inf)

  gst_reg_kegg <- goregion(regions = results.rangesAll, all.cpg = colnames(DNAm_all), array.type = "EPIC", collection = "KEGG")
  topGSA(gst_reg_kegg[gst_reg_kegg$P.DE < 0.05, ], num = Inf)
}
