# Yaro
# Script: 04_wgcna_DNAm_modules.R
#
# Purpose:
#   Build WGCNA modules from selected CpGs (feature-selected CpGs),
#   relate module eigengenes to AD-related phenotypes (adjusted regression),
#   optionally run moderation (module × APOE interaction),
#   and run enrichment (GO/KEGG) for CpGs in modules of interest.
#
#
# Outputs:
#   - WGCNA module assignment CSV (cpg_assignment_PFC.csv)
#   - Module–trait heatmaps (main effects and moderation effects)
#   - Optional: enrichment tables and plots

setwd("~/Desktop/research/ROSMAP/DNAme/WGCNA/")
set.seed(0)

# ---- Libraries ----
library(WGCNA)
library(dplyr)
library(reshape2)
library(minfi)
library(missMethyl)
library(ggplot2)

allowWGCNAThreads()
options(stringsAsFactors = FALSE)

# ---- Load CpG list (feature-selected CpGs) ----
cpgs_df <- read.csv("~/Desktop/research/ROSMAP/DNAme/feature_selection/niareagansc_recode_cpgs_PFC.csv")
cpgs <- na.omit(cpgs_df$APOE_TRUE_only)

# Safety checks: CpGs exist in methylation matrix
cpgs <- intersect(cpgs, colnames(DNAm_all))
stopifnot(length(cpgs) > 0)

# ---- Align phenotype to DNAm_all ----
# Ensure big_pheno_all is matched to DNAm_all row order.
big_pheno_all <- big_pheno_all[big_pheno_all$projid %in% rownames(DNAm_all), , drop = FALSE]
DNAm_all <- DNAm_all[match(big_pheno_all$projid, rownames(DNAm_all)), , drop = FALSE]
stopifnot(identical(rownames(DNAm_all), as.character(big_pheno_all$projid)))

# ---- Expression-like matrix for WGCNA ----
# WGCNA expects samples × features. Here: samples × CpGs.
datExpr <- scale(DNAm_all[, cpgs, drop = FALSE])

# Drop any samples/CpGs with too many missing values / invalid entries.
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
  big_pheno_all <- big_pheno_all[big_pheno_all$projid %in% rownames(datExpr), , drop = FALSE]
  datExpr <- datExpr[match(big_pheno_all$projid, rownames(datExpr)), , drop = FALSE]
}
stopifnot(identical(rownames(datExpr), as.character(big_pheno_all$projid)))

nSamples <- nrow(datExpr)
nGenes <- ncol(datExpr)

# ---- Soft-threshold power selection ----
powers <- c(seq(2, 10, by = 2))

sft <- pickSoftThreshold(
  datExpr,
  powerVector = powers,
  verbose = 5,
  networkType = "signed"
)

# Plot scale-free fit and mean connectivity to choose power.
sizeGrWindow(9, 5)
par(mfrow = c(1, 2))
cex1 <- 0.9

plot(
  sft$fitIndices[, 1],
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  xlab = "Soft Threshold (power)",
  ylab = "Scale Free Topology Model Fit, signed R^2",
  type = "n",
  main = "Scale independence"
)
text(
  sft$fitIndices[, 1],
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  labels = powers,
  cex = cex1,
  col = "red"
)
abline(h = 0.90, col = "red")

plot(
  sft$fitIndices[, 1],
  sft$fitIndices[, 5],
  xlab = "Soft Threshold (power)",
  ylab = "Mean Connectivity",
  type = "n",
  main = "Mean connectivity"
)
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")

# ---- Network construction ----
chosen_power <- 6

net_PFC_DNAme <- blockwiseModules(
  datExpr,
  power = chosen_power,
  networkType = "signed",
  minModuleSize = 30,
  maxBlockSize = max(1000, nGenes),
  numericLabels = TRUE,
  saveTOMs = FALSE,
  saveTOMFileBase = "PFC_AD_anycog_TOM",
  verbose = 3,
  nThreads = 8
)

# Module colors and eigengenes
moduleLabels_DNAme <- net_PFC_DNAme$colors
moduleColors_DNAme <- labels2colors(moduleLabels_DNAme)

MEs_DNAme0 <- moduleEigengenes(datExpr, colors = moduleColors_DNAme)$eigengenes
MEs_DNAme <- orderMEs(MEs_DNAme0)

geneTree_DNAme <- net_PFC_DNAme$dendrograms[[1]]

# Module membership table
table(moduleColors_DNAme)

# ---- Optional: override module assignments from saved CSV ----
use_saved_assignments <- FALSE

if (use_saved_assignments) {
  module_assignments <- read.csv("./cpg_assignment_PFC.csv", stringsAsFactors = FALSE)
  stopifnot(all(module_assignments$probes %in% colnames(datExpr)))
  moduleColors_DNAme <- module_assignments$moduleColor
  names(moduleColors_DNAme) <- module_assignments$probes
  moduleColors_DNAme <- moduleColors_DNAme[colnames(datExpr)]
  MEs_DNAme0 <- moduleEigengenes(datExpr, colors = moduleColors_DNAme)$eigengenes
  MEs_DNAme <- orderMEs(MEs_DNAme0)
}

# ---- Phenotype engineering (binary encodings you used later) ----
big_pheno_all$braaksc_num <- ifelse(big_pheno_all$braaksc %in% c(3,4), 1,
                                   ifelse(big_pheno_all$braaksc %in% c(1,2), 0, NA))
big_pheno_all$braaksc_num <- as.factor(big_pheno_all$braaksc_num)

big_pheno_all$apoe_genotype_num <- ifelse(big_pheno_all$apoe_genotype == TRUE, 1,
                                         ifelse(big_pheno_all$apoe_genotype == FALSE, 0, NA))
big_pheno_all$apoe_genotype_num <- as.factor(big_pheno_all$apoe_genotype_num)

big_pheno_all$ceradsc_num <- ifelse(big_pheno_all$ceradsc == 2, 1,
                                   ifelse(big_pheno_all$ceradsc == 1, 0, NA))
big_pheno_all$ceradsc_num <- as.factor(big_pheno_all$ceradsc_num)

big_pheno_all$ADvsNorm <- ifelse(big_pheno_all$cogdx == 3, 1,
                                ifelse(big_pheno_all$cogdx == 1, 0, NA))
big_pheno_all$ADvsNorm <- as.factor(big_pheno_all$ADvsNorm)

# Ensure phenotype data matches datExpr samples exactly
big_pheno_all <- big_pheno_all[match(rownames(datExpr), big_pheno_all$projid), , drop = FALSE]
stopifnot(identical(rownames(datExpr), as.character(big_pheno_all$projid)))

# Scale numeric columns (optional; you did this in your original)
numeric_columns <- colnames(big_pheno_all)[sapply(big_pheno_all, is.numeric)]
numeric_columns <- setdiff(numeric_columns, "projid")
big_pheno_all[, numeric_columns] <- scale(big_pheno_all[, numeric_columns])

# ---- Module–trait association via adjusted linear models ----
phenotypes <- c("niareagansc_recode", "ADvsNorm",
                "apoe_genotype_num", "ceradsc_num", "braaksc_num")

covars <- c("Female", "Age", "pmi", "neurons_pfc")

# Fit: scale(ME) ~ trait + covariates
results_main <- list()

for (trait in phenotypes) {
  for (module in names(MEs_DNAme)) {

    df <- cbind(MEs_DNAme, big_pheno_all)

    # Build model formula
    fml <- as.formula(
      paste0("scale(", module, ") ~ ", trait, " + ", paste(covars, collapse = " + "))
    )

    model <- lm(fml, data = df)
    sm <- summary(model)

    # Determine coefficient name for trait
    # If trait is factor, coefficient often ends with '1' when coded 0/1 as factor.
    # Keep your existing assumption but guard in case it differs.
    if (is.factor(df[[trait]])) {
      coef_name <- paste0(trait, "1")
    } else {
      coef_name <- trait
    }

    if (coef_name %in% rownames(sm$coefficients)) {
      results_main[[paste(module, trait, sep = "_")]] <- c(
        sm$coefficients[coef_name, ],
        module = module,
        trait = trait
      )
    }
  }
}

results_main_df <- as.data.frame(do.call(rbind, results_main))
results_main_df[, 1:4] <- sapply(results_main_df[, 1:4], as.numeric)

# Create estimate and p-value matrices for heatmap
est_main <- acast(results_main_df, module ~ trait, value.var = "Estimate")
p_main <- acast(results_main_df, module ~ trait, value.var = "Pr(>|t|)")

# Remove grey module
est_main <- est_main[rownames(est_main) != "MEgrey", , drop = FALSE]
p_main <- p_main[rownames(p_main) != "MEgrey", , drop = FALSE]

# Reorder columns according to your phenotype order
est_main <- est_main[, phenotypes, drop = FALSE]
p_main <- p_main[, phenotypes, drop = FALSE]

# Cluster modules based on their estimate patterns
row_hclust <- hclust(dist(est_main, method = "euclidean"), method = "complete")
row_order <- order.dendrogram(as.dendrogram(row_hclust))
est_main <- est_main[row_order, , drop = FALSE]
p_main <- p_main[row_order, , drop = FALSE]

# Text matrix: estimate + p-value
text_main <- paste0(signif(est_main, 2), "\n(", signif(p_main, 2), ")")
dim(text_main) <- dim(est_main)

# Plot
par(mar = c(11, 2, 4, 1))
max_val <- 1  # keep your fixed scaling for comparability
labeledHeatmap(
  Matrix = est_main,
  xLabels = c("Reagan", "Clinical dx", "APOE4", "CERAD", "Braak"),
  yLabels = rownames(est_main),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = text_main,
  setStdMargins = FALSE,
  cex.text = 1.4,
  cex.lab = 2,
  cex.main = 1.8,
  zlim = c(-max_val, max_val),
  main = "Module–trait relationships adjusted for Sex, Age, PMI, Neurons (PFC)"
)

# ---- Moderation analysis: trait ~ ME * apoe_genotype + covariates ----
# This tests whether the association between module eigengene and trait differs by APOE carrier status.

binary_traits_for_moderation <- phenotypes  # keep your list, but see suggestions below

results_mod <- list()

for (trait in binary_traits_for_moderation) {
  for (module in names(MEs_DNAme)) {

    df <- cbind(MEs_DNAme, big_pheno_all)

    fml <- as.formula(
      paste0(trait, " ~ scale(", module, ") * apoe_genotype + ",
             paste(covars, collapse = " + "))
    )

    model <- glm(fml, data = df, family = binomial(link = "logit"))
    sm <- summary(model)

    coef_name <- paste0("scale(", module, "):apoe_genotypeTRUE")

    if (coef_name %in% rownames(sm$coefficients)) {
      results_mod[[paste(module, trait, sep = "_")]] <- c(
        sm$coefficients[coef_name, ],
        module = module,
        trait = trait
      )
    }
  }
}

results_mod_df <- as.data.frame(do.call(rbind, results_mod))
results_mod_df[, 1:4] <- sapply(results_mod_df[, 1:4], as.numeric)

est_mod <- acast(results_mod_df, module ~ trait, value.var = "Estimate")
p_mod <- acast(results_mod_df, module ~ trait, value.var = "Pr(>|z|)")

est_mod <- est_mod[rownames(est_mod) != "MEgrey", , drop = FALSE]
p_mod <- p_mod[rownames(p_mod) != "MEgrey", , drop = FALSE]

est_mod <- est_mod[, phenotypes, drop = FALSE]
p_mod <- p_mod[, phenotypes, drop = FALSE]

row_hclust <- hclust(dist(est_mod, method = "euclidean"), method = "complete")
row_order <- order.dendrogram(as.dendrogram(row_hclust))
est_mod <- est_mod[row_order, , drop = FALSE]
p_mod <- p_mod[row_order, , drop = FALSE]

text_mod <- paste0(signif(est_mod, 2), "\n(", signif(p_mod, 2), ")")
dim(text_mod) <- dim(est_mod)

par(mar = c(11, 2, 4, 1))
max_val_mod <- max(abs(est_mod), na.rm = TRUE)

labeledHeatmap(
  Matrix = est_mod,
  xLabels = c("Reagan", "Clinical+Reagan", "APOE4", "CERAD", "Braak"),
  yLabels = rownames(est_mod),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = text_mod,
  setStdMargins = FALSE,
  cex.text = 1.4,
  cex.lab = 2,
  cex.main = 1.8,
  zlim = c(-max_val_mod, max_val_mod),
  main = "Moderation: (Module eigengene × APOE) effect on trait"
)

# ---- Enrichment for CpGs in modules of interest ----
# Example: choose module colors of interest.
modules_of_interest <- c("pink")

cpgs_to_test <- colnames(datExpr)[moduleColors_DNAme %in% modules_of_interest]

if (length(cpgs_to_test) > 0) {

  gst_go <- gometh(sig.cpg = cpgs_to_test, all.cpg = colnames(datExpr), array.type = "EPIC")

  # Filter and sort for export/plotting
  filtered_gst <- gst_go[gst_go$FDR <= 0.1, ]
  sorted_gst <- filtered_gst[order(filtered_gst$FDR), ]
  top_gst <- head(sorted_gst, 10)

  # Save enrichment results
  write.csv(top_gst, "./enrichment_moderator_positive_PFC.csv", row.names = FALSE)

  # Plot top terms (note: your original y-axis label was -log10(FDR) but you plotted FDR directly)
  ggplot(top_gst, aes(x = reorder(TERM, FDR), y = -log10(FDR))) +
    geom_bar(stat = "identity", color = "black") +
    coord_flip() +
    labs(x = "GO Term", y = "-log10(FDR)") +
    theme_minimal()
}

# ---- Save CpG assignment table (module membership + trait significance) ----
# This recreates the canonical WGCNA output structure used for hub selection and follow-up.

niareagansc_recode_num <- as.data.frame(as.numeric(big_pheno_all$niareagansc_recode))
names(niareagansc_recode_num) <- "niareagansc_recode"

modNames <- substring(names(MEs_DNAme), 3)

geneModuleMembership <- as.data.frame(cor(datExpr, MEs_DNAme, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) <- paste0("MM", modNames)
names(MMPvalue) <- paste0("p.MM", modNames)

geneTraitSignificance <- as.data.frame(cor(datExpr, niareagansc_recode_num, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) <- "GS.niareagansc_recode"
names(GSPvalue) <- "p.GS.niareagansc_recode"

annEPIC <- getAnnotation("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
probes <- colnames(datExpr)
probes2annEPIC <- match(probes, annEPIC$Name)

geneInfo0 <- data.frame(
  probes = probes,
  geneSymbol = annEPIC$UCSC_RefGene_Name[probes2annEPIC],
  RefGene_ID = annEPIC$UCSC_RefGene_Accession[probes2annEPIC],
  moduleColor = moduleColors_DNAme,
  geneTraitSignificance,
  GSPvalue,
  stringsAsFactors = FALSE
)

# Order modules by their correlation with the trait (optional)
modOrder <- order(-abs(cor(MEs_DNAme, niareagansc_recode_num, use = "p")))

for (mod in 1:ncol(geneModuleMembership)) {
  oldnames <- names(geneInfo0)
  geneInfo0 <- data.frame(
    geneInfo0,
    geneModuleMembership[, modOrder[mod]],
    MMPvalue[, modOrder[mod]],
    stringsAsFactors = FALSE
  )
  names(geneInfo0) <- c(
    oldnames,
    paste0("MM.", modNames[modOrder[mod]]),
    paste0("p.MM.", modNames[modOrder[mod]])
  )
}

geneOrder <- order(geneInfo0$moduleColor, -abs(geneInfo0$GS.niareagansc_recode))
geneInfo <- geneInfo0[geneOrder, ]

write.csv(geneInfo, file = "./cpg_assignment_PFC.csv", row.names = FALSE)
