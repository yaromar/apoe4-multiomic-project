# =============================================================================
# 03_ewas_AD_status.R
# Yaro
# Last update: 12/16/23
#
# Reference tutorial (conceptual): 
# https://nbis-workshop-epigenomics.readthedocs.io/en/latest/content/tutorials/methylationArray/Array_Tutorial.html
#
# PURPOSE
# -------
# Identify differentially methylated probes (DMPs) associated with AD pathology
# in ROSMAP PFC methylation data using limma (probe-wise linear modeling).
#
# Additionally:
# - remove sex chromosomes
# - remove outlier samples within AD and Control groups using WGCNA clustering
# - filter low-variance CpGs (bottom 10% SD)
# - generate Manhattan + Volcano plots
# - optional APOE+ subgroup EWAS
# - optional APOE interaction sensitivity scan (logistic regression across CpGs)
# - enrichment analyses (gometh / KEGG) and DMRcate region analysis
# - Venn diagram across regions (if you have region-specific CpG files)
#
# IMPORTANT INPUTS ASSUMED TO EXIST IN THE ENVIRONMENT
# ---------------------------------------------------
# big_pheno_scaled : phenotype table containing at least:
#   projid, niareagansc_recode, Age, Female, neurons_pfc, pmi,
#   apoe_genotype (logical or factor), apoe_long (string genotype code)
#
# PFCSamples : matrix/data.frame of methylation beta values with:
#   rows = sample IDs matching projid
#   cols = CpG IDs matching EPIC annotation (annEPIC$Name)
#
# =============================================================================

# -----------------------------------------------------------------------------
# Libraries
# -----------------------------------------------------------------------------
library(limma)
library(parallel)
library(BiocParallel)
library(minfi)
library(missMethyl)
library(DMRcate)
library(GenomicRanges)

library(WGCNA)       # goodSamplesGenes + cutreeStatic + clustering plots
library(dplyr)
library(stringr)

library(qqman)       # Manhattan plot
library(ggplot2)
library(ggrepel)

library(VennDiagram)
library(grid)        # grid.draw for VennDiagram output


library(MASS)       # rlm
library(ggpubr)     # stat_compare_means
library(ggeffects)  # ggpredict

set.seed(0)

# -----------------------------------------------------------------------------
# Working directory + annotation
# -----------------------------------------------------------------------------
setwd("~/Desktop/research/ROSMAP/DNAme/feature_selection/")

annEPIC <- getAnnotation("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")

# -----------------------------------------------------------------------------
# Parallel settings (BiocParallel for some Bioconductor functions)
# -----------------------------------------------------------------------------
numCores <- detectCores()
register(BiocParallel::MulticoreParam(numCores - 1))

# -----------------------------------------------------------------------------
# STEP 1: Assemble full analysis set (use all available data)
# -----------------------------------------------------------------------------
big_pheno_all <- big_pheno_scaled

DNAm_all <- PFCSamples[(rownames(PFCSamples) %in% big_pheno_all$projid), ]
DNAm_all <- DNAm_all[match(big_pheno_all$projid, rownames(DNAm_all)), ]

# -----------------------------------------------------------------------------
# STEP 2: Remove sex chromosomes (chrX and/or chrY)
# -----------------------------------------------------------------------------
DNAm_all <- DNAm_all[, !(colnames(DNAm_all) %in% rownames(annEPIC[annEPIC$chr == "chrX", ]))]
DNAm_all <- DNAm_all[, !(colnames(DNAm_all) %in% rownames(annEPIC[annEPIC$chr == "chrY", ]))]

# -----------------------------------------------------------------------------
# STEP 3: WGCNA goodSamplesGenes QC
# -----------------------------------------------------------------------------
allowWGCNAThreads()
options(stringsAsFactors = FALSE)

gsg <- goodSamplesGenes(DNAm_all, verbose = 3)
gsg$allOK

if (!gsg$allOK) {
  if (sum(!gsg$goodGenes) > 0) {
    printFlush(paste("Removing genes:", paste(colnames(DNAm_all)[!gsg$goodGenes], collapse = ", ")))
  }
  if (sum(!gsg$goodSamples) > 0) {
    printFlush(paste("Removing samples:", paste(rownames(DNAm_all)[!gsg$goodSamples], collapse = ", ")))
  }
  DNAm_all <- DNAm_all[gsg$goodSamples, gsg$goodGenes]
}

# -----------------------------------------------------------------------------
# STEP 4: Split by AD/Control (based on niareagansc_recode)
# -----------------------------------------------------------------------------
DNAm_allAD      <- DNAm_all[big_pheno_all$niareagansc_recode == 1, ]
DNAm_allControl <- DNAm_all[big_pheno_all$niareagansc_recode == 0, ]

# -----------------------------------------------------------------------------
# STEP 5: Outlier removal via clustering (hard cut heights preserved)
# -----------------------------------------------------------------------------
cluster_filter_samples <- function(mat, cutHeight, keep_rule, plot_title, abline_h = NULL) {
  sampleTree <- hclust(dist(mat), method = "average")
  
  # Plot 
  sizeGrWindow(12, 9)
  par(cex = 0.6)
  par(mar = c(0, 4, 2, 0))
  plot(sampleTree, main = plot_title, sub = "", xlab = "",
       cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
  
  # Diagnostics 
  print(summary(sampleTree$height))
  print(sd(sampleTree$height))
  print(quantile(sampleTree$height, 0.95))
  
  if (!is.null(abline_h)) abline(h = abline_h, col = "red")
  
  clust <- cutreeStatic(sampleTree, cutHeight = cutHeight, minSize = 1)
  print(table(clust))
  
  keepSamples <- keep_rule(clust)
  mat[keepSamples, , drop = FALSE]
}

# AD group: play with the parameters
DNAm_allAD <- cluster_filter_samples(
  mat = DNAm_allAD,
  cutHeight = 85,
  keep_rule = function(cl) (cl == 1),
  plot_title = "Sample clustering to detect outliers (AD)",
  abline_h = 55
)

# Control group: play with the parameters
DNAm_allControl <- cluster_filter_samples(
  mat = DNAm_allControl,
  cutHeight = 55,
  keep_rule = function(cl) (cl %in% c(1, 2)),
  plot_title = "Sample clustering to detect outliers (Control)",
  abline_h = 50
)

# Recombine after outlier removal
DNAm_all <- rbind(DNAm_allAD, DNAm_allControl)

# Re-sync phenotype rows to methylation rows (keep same ordering logic)
big_pheno_all <- big_pheno_all[big_pheno_all$projid %in% rownames(DNAm_all), ]
DNAm_all <- DNAm_all[as.character(big_pheno_all$projid), ]
rownames(DNAm_all) == big_pheno_all$projid

# -----------------------------------------------------------------------------
# STEP 6: Variance filter (bottom 10% SD removed) (kept)
# -----------------------------------------------------------------------------
sd_cpgs <- apply(DNAm_all, 2, sd)
variance_threshold <- quantile(sd_cpgs, 0.1)
cpgs_to_keep <- which(sd_cpgs > variance_threshold)
DNAm_all <- DNAm_all[, cpgs_to_keep]

# spot-check a CpG
DNAm_all[, "cg06329447"]

# -----------------------------------------------------------------------------
# STEP 7: Restrict to specific APOE genotypes via apoe_long 
# -----------------------------------------------------------------------------
target_genotypes <- c("33", "34", "44")  # Excluding "24"

print(unique(big_pheno_all$apoe_long))

big_pheno_all <- big_pheno_all %>%
  filter(apoe_long %in% target_genotypes)

DNAm_all <- DNAm_all[rownames(DNAm_all) %in% big_pheno_all$projid, ]
DNAm_all <- DNAm_all[big_pheno_all$projid, ]

dim(big_pheno_all)
dim(DNAm_all)
table(big_pheno_all$apoe_long)

# -----------------------------------------------------------------------------
# STEP 8: Impute missing pmi 
# -----------------------------------------------------------------------------
big_pheno_all[is.na(big_pheno_all[, "pmi"]), "pmi"] <- mean(big_pheno_all[, "pmi"], na.rm = TRUE)

# -----------------------------------------------------------------------------
# STEP 9: limma EWAS (AD vs Control via niareagansc_recode + covariates)
# -----------------------------------------------------------------------------
designAll <- model.matrix(~ niareagansc_recode + Age + Female + neurons_pfc + pmi, data = big_pheno_all)

fit <- lmFit(t(DNAm_all), designAll)
fit <- eBayes(fit)

# Attach annotation rows matching probes present in fit
anno_match <- annEPIC[match(rownames(fit), annEPIC$Name), c(1:4, 12:19, 24:ncol(annEPIC))]

results_niareagansc_recode_all <- topTable(
  fit,
  num = Inf,
  coef = "niareagansc_recode1",
  genelist = anno_match
)

rownames(results_niareagansc_recode_all) <- results_niareagansc_recode_all$Name

# Quick check for “enough signal” for enrichment (kept threshold)
results_niareagansc_recode_all[results_niareagansc_recode_all$adj.P.Val <= 0.11, ]

# -----------------------------------------------------------------------------
# STEP 10: Build "Unique_Genes" labels for significant entries 
# -----------------------------------------------------------------------------
significant_entries <- results_niareagansc_recode_all %>%
  filter(adj.P.Val <= 0.11) %>%
  mutate(Combined_Names = paste(GencodeBasicV12_NAME, GencodeCompV12_NAME, sep = ";")) %>%
  mutate(
    Unique_Genes = sapply(str_split(Combined_Names, ";"), function(x) {
      cleaned_names <- str_replace_all(x, "\\bRP\\d+-.+", "")
      unique_names <- unique(cleaned_names)
      final_names <- paste(unique_names, collapse = ";")
      str_remove_all(final_names, "^;+")
    })
  )

results_niareagansc_recode_all <- left_join(
  results_niareagansc_recode_all,
  significant_entries[, c("Name", "Unique_Genes")],
  by = "Name"
)

# -----------------------------------------------------------------------------
# STEP 11: Manhattan plot 
# -----------------------------------------------------------------------------
results_niareagansc_recode_all$chr_numeric <- sub("chr", "", results_niareagansc_recode_all$chr)
results_niareagansc_recode_all <- results_niareagansc_recode_all[results_niareagansc_recode_all$chr != "chrY", ]
results_niareagansc_recode_all$chr_numeric <- as.numeric(results_niareagansc_recode_all$chr_numeric)

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
  main = "A. Differential Methylation Manhattan Plot"
)

# Second line for FDR 0.05 (kept numeric)
abline(h = -log10(4.86e-7), col = "red", lty = 3)

dev.off()

# -----------------------------------------------------------------------------
# STEP 12: Volcano plot (keep one version; your script had two near-duplicates)
# -----------------------------------------------------------------------------
volcano_plot <- ggplot(results_niareagansc_recode_all, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(alpha = 0.4, color = "gray40") +
  labs(
    x = expression(Delta~"β (AD - Control)"),
    y = expression(-log[10](italic(P))),
    title = "B. Volcano Plot of Differential Methylation"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(3.61e-6), linetype = "dashed", color = "red", linewidth = 0.5) +
  geom_hline(yintercept = -log10(4.86e-7), linetype = "dotted", color = "red", linewidth = 0.5) +
  geom_text_repel(
    aes(
      label = ifelse(P.Value <= 3.61e-6, Unique_Genes, ""),
      color = ifelse(logFC >= 0, "red", "blue")
    ),
    segment.color = "grey50",
    size = 2,
    show.legend = FALSE,
    force = 2,
    force_pull = 0.5,
    max.overlaps = 30
  ) +
  scale_color_identity() +
  theme_minimal(base_size = 12)

volcano_plot
ggsave("Volcano_Plot_PFC.pdf", plot = volcano_plot, width = 10, height = 6, units = "in", dpi = 300)


# -----------------------------------------------------------------------------
# STEP 13: Optional APOE+ subgroup EWAS (kept structure)
# -----------------------------------------------------------------------------
big_pheno_APOE <- big_pheno_scaled[big_pheno_scaled$apoe_genotype == TRUE, ]

DNAm_APOE <- PFCSamples[(rownames(PFCSamples) %in% big_pheno_APOE$projid), ]
DNAm_APOE <- DNAm_APOE[match(big_pheno_APOE$projid, rownames(DNAm_APOE)), ]

# Subset to samples that survived earlier filtering
big_pheno_APOE <- big_pheno_APOE[big_pheno_APOE$projid %in% big_pheno_all$projid, ]
DNAm_APOE <- DNAm_APOE[rownames(DNAm_APOE) %in% rownames(DNAm_all), ]

rownames(DNAm_APOE) == big_pheno_APOE$projid

sd_cpgs <- apply(DNAm_APOE, 2, sd)
variance_threshold <- quantile(sd_cpgs, 0.1)
cpgs_to_keep <- which(sd_cpgs > variance_threshold)
DNAm_APOE <- DNAm_APOE[, cpgs_to_keep]

big_pheno_APOE[is.na(big_pheno_APOE[, "pmi"]), "pmi"] <- mean(big_pheno_APOE[, "pmi"], na.rm = TRUE)

designAPOE <- model.matrix(~ niareagansc_recode + Age + Female + neurons_pfc + pmi, data = big_pheno_APOE)

fit <- lmFit(t(DNAm_APOE), designAPOE)
fit <- eBayes(fit)

anno_match_APOE <- annEPIC[match(rownames(fit), annEPIC$Name), c(1:4, 12:19, 24:ncol(annEPIC))]

results_niareagansc_recode_APOE <- topTable(
  fit,
  num = Inf,
  coef = "niareagansc_recode1",
  genelist = anno_match_APOE
)

results_niareagansc_recode_APOE[results_niareagansc_recode_APOE$adj.P.Val < 0.06, ]
results_niareagansc_recode_APOE <- results_niareagansc_recode_APOE[results_niareagansc_recode_APOE$adj.P.Val < 0.06, ]

cpg_niareagansc_recode_APOE <- rownames(results_niareagansc_recode_APOE)

# -----------------------------------------------------------------------------
# STEP 14: Optional APOE interaction scan (logistic regression per CpG)
# -----------------------------------------------------------------------------
# Prepare statuses (kept)
big_pheno_all$AD_status <- as.factor(big_pheno_all$niareagansc_recode)
big_pheno_all$APOE_status <- as.factor(big_pheno_all$apoe_genotype)
big_pheno_all$pmi[is.na(big_pheno_all$pmi)] <- mean(big_pheno_all$pmi, na.rm = TRUE)

numCores <- detectCores() - 1
cl <- makeCluster(numCores)

clusterExport(cl, list("DNAm_all", "big_pheno_all"))
clusterEvalQ(cl, { library(stats) })

perform_analysis <- function(cpg_id) {
  molecule_values <- DNAm_all[, cpg_id]

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

  # Keep original coefficient row names (APOE_statusTRUE etc.)
  list(
    CpG = cpg_id,
    Beta_APOE = coefs["APOE_statusTRUE", "Estimate"],
    SE_APOE = coefs["APOE_statusTRUE", "Std. Error"],
    P_APOE = coefs["APOE_statusTRUE", "Pr(>|z|)"],
    Beta_Molecule = coefs["Molecule", "Estimate"],
    SE_Molecule = coefs["Molecule", "Std. Error"],
    P_Molecule = coefs["Molecule", "Pr(>|z|)"],
    Beta_Interaction = coefs["APOE_statusTRUE:Molecule", "Estimate"],
    SE_Interaction = coefs["APOE_statusTRUE:Molecule", "Std. Error"],
    P_Interaction = coefs["APOE_statusTRUE:Molecule", "Pr(>|z|)"]
  )
}

results <- parLapply(cl, colnames(DNAm_all), perform_analysis)
stopCluster(cl)

methylation_interaction_results <- do.call(rbind.data.frame, results)
methylation_interaction_results$Adj_P_Interaction <- p.adjust(
  methylation_interaction_results$P_Interaction,
  method = "BH"
)

# Merge “main” EWAS output with interaction p-values (kept structure)
significant_main_results <- results_niareagansc_recode_all

significant_sensitivity_results <- methylation_interaction_results %>%
  mutate(APOE_Interaction_Significance = P_Interaction) %>%
  mutate(Name = CpG) %>%
  dplyr::select(Name, APOE_Interaction_Significance)

combined_results <- significant_main_results %>%
  full_join(significant_sensitivity_results, by = "Name")

combined_results <- combined_results %>%
  mutate(Combined_Names = paste(GencodeBasicV12_NAME, GencodeCompV12_NAME, sep = ";")) %>%
  mutate(
    Unique_Genes = sapply(str_split(Combined_Names, ";"), function(x) {
      cleaned_names <- str_replace_all(x, "\\bRP\\d+-.+", "")
      unique_names <- unique(cleaned_names)
      final_names <- paste(unique_names, collapse = ";")
      str_remove_all(final_names, "^;+")
    })
  )



# -----------------------------------------------------------------------------
# STEP 15: Enrichment (gometh / KEGG) (kept thresholds)
# -----------------------------------------------------------------------------
gst <- gometh(
  sig.cpg = results_niareagansc_recode_all[results_niareagansc_recode_all$adj.P.Val < 0.06, ]$Name,
  all.cpg = colnames(DNAm_all),
  array.type = "EPIC"
)
topGSA(gst[gst$P.DE < 0.05, ])

gst <- gometh(
  sig.cpg = results_niareagansc_recode_all[results_niareagansc_recode_all$adj.P.Val < 0.06, ]$Name,
  all.cpg = colnames(DNAm_all),
  array.type = "EPIC",
  collection = "KEGG"
)
topGSA(gst[gst$P.DE < 0.05, ])

# -----------------------------------------------------------------------------
# STEP 16: DMRcate region analysis + region enrichment (kept)
# -----------------------------------------------------------------------------
myAnnotationAPOE <- cpg.annotate(
  object = t(DNAm_all),
  datatype = "array",
  what = "Beta",
  analysis.type = "differential",
  design = designAll,
  contrasts = FALSE,
  coef = "niareagansc_recode1",
  arraytype = "EPIC"
)

DMRsAPOE <- dmrcate(myAnnotationAPOE)
results.rangesAPOE <- extractRanges(DMRsAPOE)
results.rangesAPOE

gst <- goregion(
  regions = results.rangesAPOE,
  all.cpg = colnames(DNAm_all),
  array.type = "EPIC"
)
topGSA(gst[gst$P.DE < 0.05 & gst$ONTOLOGY == "BP", ], num = Inf)

gst <- goregion(
  regions = results.rangesAPOE,
  all.cpg = colnames(DNAm_all),
  array.type = "EPIC",
  collection = "KEGG"
)
topGSA(gst[gst$P.DE < 0.05, ], num = Inf)

# -----------------------------------------------------------------------------
# STEP 17: Venn diagram across region CpG lists
# -----------------------------------------------------------------------------
# NOTE: Your original code read PFC twice and labeled it twice.
# If that duplication was intentional for presentation, revert to the original.
cpgs_pfc <- read.csv("~/Desktop/research/ROSMAP/DNAme/feature_selection/niareagansc_recode_cpgs_PFC.csv")
cpgs_pfc <- na.omit(cpgs_pfc$APOE_TRUE_only)

cpgs_st <- read.csv("~/Desktop/research/ROSMAP/DNAme/feature_selection/niareagansc_recode_cpgs_ST.csv")
cpgs_st <- na.omit(cpgs_st$APOE_TRUE_only)

venn.plot <- venn.diagram(
  x = list(
    PFC = cpgs_pfc,
    ST  = cpgs_st
  ),
  filename = NULL,
  output = TRUE,
  height = 800,
  width = 800,
  resolution = 300,
  compression = "lzw",
  lwd = 2,
  lty = "blank",
  fill = c("cornflowerblue", "yellow"),
  alpha = 0.50,
  cex = 1.5,
  fontfamily = "serif",
  fontface = "bold",
  cat.col = c("darkblue", "darkred"),
  cat.cex = 1.2,
  cat.fontfamily = "serif",
  cat.dist = c(0.06, 0.06),
  cat.pos = c(-27, 135),
  cat.default.pos = "outer"
)

grid.draw(venn.plot)
