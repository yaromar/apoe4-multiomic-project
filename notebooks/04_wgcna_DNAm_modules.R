# =============================================================================
# 04_wgcna_DNAm_modules.R
#
# PURPOSE
# -------
# Run WGCNA on selected CpGs,
# define co-methylation modules, compute module eigengenes, and test:
#   (A) module ~ trait associations (linear models; adjusted for covariates)
#   (B) moderation: trait ~ module * APOE4 (logistic models; interaction term)
#
# Then:
# - module GO/KEGG enrichment via missMethyl::gometh
# - compute CpG "module membership" and "gene significance"
# - export CpG-to-module assignment file
# - optional: load TOM and compute hub CpGs (top 10% TOM-sum) for enrichment
#
# =============================================================================

setwd("~/Desktop/research/ROSMAP/DNAme/WGCNA/")
set.seed(0)

# -----------------------------------------------------------------------------
# Libraries
# -----------------------------------------------------------------------------
library(WGCNA)
library(dplyr)
library(reshape2)
library(missMethyl)
library(minfi)
library(VennDiagram)
library(grid)

allowWGCNAThreads()
options(stringsAsFactors = FALSE)

# -----------------------------------------------------------------------------
# STEP 1: Choose CpGs used for WGCNA
# -----------------------------------------------------------------------------
# Feature-selected CpGs from prior EWAS step:
cpgs <- read.csv("~/Desktop/research/ROSMAP/DNAme/feature_selection/niareagansc_recode_cpgs_PFC.csv")
cpgs <- na.omit(cpgs$APOE_TRUE_only)

# Alternative: CpGs significant for APOE interaction scan
# cpgs <- methylation_interaction_results[methylation_interaction_results$P_Interaction < 0.05, ]$Protein

# -----------------------------------------------------------------------------
# STEP 2: Build expression-like matrix for WGCNA (CpGs as "genes")
# -----------------------------------------------------------------------------
# Keeping exactly that.
datExpr <- scale(DNAm_all[, cpgs])

# Ensure no missingness / WGCNA compatibility (optional sanity check)
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

nSamples <- nrow(datExpr)
nGenes   <- ncol(datExpr)

# -----------------------------------------------------------------------------
# STEP 3: Pick soft-threshold
# -----------------------------------------------------------------------------
powers <- c(seq(2, 4, by = 2))  # adjust the range

sft <- pickSoftThreshold(
  datExpr,
  powerVector = powers,
  verbose = 5,
  networkType = "signed",
  blockSize = 32000
)

# Plot scale independence and mean connectivity
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

# -----------------------------------------------------------------------------
# STEP 4: Build modules
# -----------------------------------------------------------------------------
net_PFC_DNAme <- blockwiseModules(
  datExpr,
  power = 6,                      # 6 for regular; 4 for moderated 
  networkType = "signed",
  minModuleSize = 30,
  maxBlockSize = 34000,
  numericLabels = TRUE,
  # saveTOMs = TRUE,
  saveTOMFileBase = "PFC_AD_anycog_TOM",
  verbose = 3,
  nThreads = 8
)

table(net_PFC_DNAme$colors)

# Convert numeric module labels to color names
moduleLabels_DNAme <- net_PFC_DNAme$colors
moduleColors_DNAme <- labels2colors(net_PFC_DNAme$colors)

# Dendrogram for reference
geneTree_DNAme <- net_PFC_DNAme$dendrograms[[1]]

# Compute module eigengenes from moduleColors_DNAme 
MEs_DNAme0 <- moduleEigengenes(datExpr, colors = moduleColors_DNAme)$eigengenes
MEs_DNAme  <- orderMEs(MEs_DNAme0)

# -----------------------------------------------------------------------------
# OPTIONAL: Override module assignments from an external file
# -----------------------------------------------------------------------------
use_external_module_assignments <- FALSE

if (use_external_module_assignments) {
  module_assignments <- read.csv("./cpg_assignment_PFC.csv", stringsAsFactors = FALSE)

  # Sanity checks
  stopifnot(all(module_assignments$probes %in% colnames(datExpr)))
  stopifnot(all(colnames(datExpr) %in% module_assignments$probes))

  moduleColors_DNAme <- module_assignments$moduleColor
  names(moduleColors_DNAme) <- module_assignments$probes
  moduleColors_DNAme <- moduleColors_DNAme[colnames(datExpr)]

  MEs_DNAme0 <- moduleEigengenes(datExpr, colors = moduleColors_DNAme)$eigengenes
  MEs_DNAme  <- orderMEs(MEs_DNAme0)
}

# -----------------------------------------------------------------------------
# STEP 5: Prepare phenotype data for modeling 
# -----------------------------------------------------------------------------
# Ensure big_pheno_all is aligned to datExpr rows 
big_pheno_all <- big_pheno_all[rownames(datExpr), ]

# Scale numeric phenotypes 
numeric_columns <- colnames(big_pheno_all)[sapply(big_pheno_all, is.numeric)]
if (length(numeric_columns) >= 2) {
  numeric_columns_to_scale <- numeric_columns[2:length(numeric_columns)]
  big_pheno_all[, numeric_columns_to_scale] <- scale(big_pheno_all[, numeric_columns_to_scale])
}

# Create derived binary traits 
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

# Traits analyzed in module-trait regression
phenotypes <- c(
  "niareagansc_recode",
  "ADvsNorm",
  "apoe_genotype_num",
  "ceradsc_num",
  "braaksc_num"
)

# -----------------------------------------------------------------------------
# STEP 6: Linear models: scale(ME) ~ trait + covariates
# -----------------------------------------------------------------------------
# Shows how each module eigengene associates with each trait (adjusted)
results_linear <- list()

for (trait in phenotypes) {
  for (module in names(MEs_DNAme)) {

    formula_str <- paste0(
      "scale(", module, ") ~ ", trait, " + Female + Age + pmi + neurons_pfc"
    )
    model <- lm(as.formula(formula_str), data = cbind(MEs_DNAme, big_pheno_all))
    sm <- summary(model)

    coef_name <- if (is.factor(big_pheno_all[[trait]])) paste0(trait, "1") else trait

    if (coef_name %in% rownames(sm$coefficients)) {
      results_linear[[paste(module, trait, sep = "_")]] <-
        c(sm$coefficients[coef_name, ], module = module, trait = trait)
    }
  }
}

results_df <- as.data.frame(do.call(rbind, results_linear))
results_df[, 1:4] <- sapply(results_df[, 1:4], as.numeric)

results_matrix <- acast(results_df, module ~ trait, value.var = "Estimate")
pvalue_matrix  <- acast(results_df, module ~ trait, value.var = "Pr(>|t|)")

# Drop MEgrey
results_matrix <- results_matrix[rownames(results_matrix) != "MEgrey", , drop = FALSE]
pvalue_matrix  <- pvalue_matrix[rownames(pvalue_matrix)  != "MEgrey", , drop = FALSE]

# Enforce column order
results_matrix <- results_matrix[, phenotypes, drop = FALSE]
pvalue_matrix  <- pvalue_matrix[, phenotypes, drop = FALSE]

# Cluster rows (Euclidean on estimated effects)
row_hclust <- hclust(dist(results_matrix, method = "euclidean"), method = "complete")
row_order  <- order.dendrogram(as.dendrogram(row_hclust))

results_matrix <- results_matrix[row_order, , drop = FALSE]
pvalue_matrix  <- pvalue_matrix[row_order,  , drop = FALSE]

textMatrix <- paste(signif(results_matrix, 2), "\n(", signif(pvalue_matrix, 1), ")", sep = "")
dim(textMatrix) <- dim(results_matrix)

par(mar = c(11, 2, 4, 1))
max_val <- 1  # you hard-set this; keep it

labeledHeatmap(
  Matrix = results_matrix,
  xLabels = c("Reagan", "Clinical dx", "APOE4", "CERAD", "Braak"),
  yLabels = rownames(results_matrix),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 1.4,
  cex.lab  = 2,
  cex.main = 1.8,
  zlim = c(-max_val, max_val),
  main = "Module-trait relationships adjusted for Sex, Age, and pmi \n Prefrontal Cortex"
)

# -----------------------------------------------------------------------------
# STEP 7: Moderation analysis (logistic): trait ~ scale(ME) * apoe_genotype + covariates
# -----------------------------------------------------------------------------
# NOTE: This only makes sense when trait is binary and modeled with binomial.
results_mod <- list()

for (trait in phenotypes) {
  for (module in names(MEs_DNAme)) {

    formula_str <- paste0(
      trait, " ~ scale(", module, ") * apoe_genotype + Female + Age + neurons_pfc + pmi"
    )

    model <- glm(
      as.formula(formula_str),
      data = cbind(MEs_DNAme, big_pheno_all),
      family = binomial(link = "logit")
    )

    sm <- summary(model)

    coef_name <- paste0("scale(", module, "):apoe_genotypeTRUE")

    if (coef_name %in% rownames(sm$coefficients)) {
      results_mod[[paste(module, trait, sep = "_")]] <-
        c(sm$coefficients[coef_name, ], module = module, trait = trait)
    }
  }
}

results_df_mod <- as.data.frame(do.call(rbind, results_mod))
results_df_mod[, 1:4] <- sapply(results_df_mod[, 1:4], as.numeric)

results_matrix_mod <- acast(results_df_mod, module ~ trait, value.var = "Estimate")
pvalue_matrix_mod  <- acast(results_df_mod, module ~ trait, value.var = "Pr(>|z|)")

results_matrix_mod <- results_matrix_mod[rownames(results_matrix_mod) != "MEgrey", , drop = FALSE]
pvalue_matrix_mod  <- pvalue_matrix_mod[rownames(pvalue_matrix_mod)  != "MEgrey", , drop = FALSE]

results_matrix_mod <- results_matrix_mod[, phenotypes, drop = FALSE]
pvalue_matrix_mod  <- pvalue_matrix_mod[, phenotypes, drop = FALSE]

row_hclust <- hclust(dist(results_matrix_mod, method = "euclidean"), method = "complete")
row_order  <- order.dendrogram(as.dendrogram(row_hclust))

results_matrix_mod <- results_matrix_mod[row_order, , drop = FALSE]
pvalue_matrix_mod  <- pvalue_matrix_mod[row_order,  , drop = FALSE]

textMatrix <- paste(signif(results_matrix_mod, 2), "\n(", signif(pvalue_matrix_mod, 1), ")", sep = "")
dim(textMatrix) <- dim(results_matrix_mod)

par(mar = c(11, 2, 4, 1))
max_val <- max(abs(range(results_matrix_mod, na.rm = TRUE)))

labeledHeatmap(
  Matrix = results_matrix_mod,
  xLabels = c("Reagan", "Clinical+Reagan", "APOE4", "CERAD", "Braak"),
  yLabels = rownames(results_matrix_mod),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 1.4,
  cex.lab  = 2,
  cex.main = 1.8,
  zlim = c(-max_val, max_val),
  main = "Module's Moderation Effect on the \n APOE4-Trait Relationship"
)

# -----------------------------------------------------------------------------
# STEP 8: Module enrichment (GO/KEGG) using missMethyl::gometh
# -----------------------------------------------------------------------------
# Example: test a single module 
cpgs_to_test <- names(net_PFC_DNAme$colors)[moduleColors_DNAme %in% c("pink")]

gst_go <- gometh(sig.cpg = cpgs_to_test, all.cpg = colnames(datExpr), array.type = "EPIC")
gst_go[gst_go$FDR < 0.1, c("TERM", "FDR")][order(gst_go$FDR[gst_go$FDR < 0.1]), ]

filtered_gst <- gst_go[gst_go$FDR <= 0.1, ]
sorted_gst   <- filtered_gst[order(filtered_gst$FDR), ]
top_gst      <- head(sorted_gst, 10)

library(ggplot2)
ggplot(top_gst, aes(x = reorder(TERM, FDR), y = FDR)) +
  geom_bar(stat = "identity", color = "black") +
  coord_flip() +
  labs(x = "GO Term", y = "FDR") +
  theme_minimal()

write.csv(top_gst, "./enrichment_moderator_positive_PFC.csv")

# KEGG enrichment example 
gst_kegg <- gometh(sig.cpg = cpgs, all.cpg = colnames(DNAm_all), array.type = "EPIC", collection = "KEGG")
gst_kegg[gst_kegg$FDR < 0.1, c("Description", "FDR")][order(gst_kegg$FDR[gst_kegg$FDR < 0.1]), ]

# check which module is CpG in
moduleColors_DNAme[which(names(net_PFC_DNAme$colors) %in% c("cg06329447"))]

# -----------------------------------------------------------------------------
# STEP 9: Module membership (MM) and gene significance (GS) for niareagansc_recode
# -----------------------------------------------------------------------------
niareagansc_recode <- as.data.frame(as.numeric(big_pheno_all$niareagansc_recode))
names(niareagansc_recode) <- "niareagansc_recode"

modNames <- substring(names(MEs_DNAme), 3)

geneModuleMembership <- as.data.frame(cor(datExpr, MEs_DNAme, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) <- paste("MM", modNames, sep = "")
names(MMPvalue)            <- paste("p.MM", modNames, sep = "")

geneTraitSignificance <- as.data.frame(cor(datExpr, niareagansc_recode, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) <- paste("GS.", names(niareagansc_recode), sep = "")
names(GSPvalue)             <- paste("p.GS.", names(niareagansc_recode), sep = "")

# -----------------------------------------------------------------------------
# STEP 10: Annotation and export CpG assignment file
# -----------------------------------------------------------------------------
annEPIC <- getAnnotation("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
probes <- colnames(datExpr)
probes2annEPIC <- match(probes, annEPIC$Name)

# How many probes missing annotation (expect 0)
sum(is.na(probes2annEPIC))

geneInfo0 <- data.frame(
  probes = probes,
  geneSymbol = annEPIC$UCSC_RefGene_Name[probes2annEPIC],
  RefGene_ID = annEPIC$UCSC_RefGene_Accession[probes2annEPIC],
  moduleColor = moduleColors_DNAme,
  geneTraitSignificance,
  GSPvalue,
  stringsAsFactors = FALSE
)

# Order modules by correlation with trait, then add MM columns in that order
modOrder <- order(-abs(cor(MEs_DNAme, niareagansc_recode, use = "p")))

for (mod in 1:ncol(geneModuleMembership)) {
  oldnames <- names(geneInfo0)
  geneInfo0 <- data.frame(geneInfo0,
                          geneModuleMembership[, modOrder[mod]],
                          MMPvalue[, modOrder[mod]])
  names(geneInfo0) <- c(oldnames,
                        paste("MM.", modNames[modOrder[mod]], sep = ""),
                        paste("p.MM.", modNames[modOrder[mod]], sep = ""))
}

geneOrder <- order(geneInfo0$moduleColor, -abs(geneInfo0$GS.niareagansc_recode))
geneInfo <- geneInfo0[geneOrder, ]

geneInfo
write.csv(geneInfo, file = "./cpg_assignment_PFC.csv")

# -----------------------------------------------------------------------------
# STEP 11 (Optional): Hub CpG enrichment using TOM
# -----------------------------------------------------------------------------
# Requires that you saved TOMs in blockwiseModules
# load("PFC_AD_anycog_TOM-block.1.RData")
#
# NOTE: If TOM is stored as a compact vector (lower triangle), we can use listToMatrix().

# load("PFC_AD_anycog_TOM-block.1.RData")
#
# listToMatrix <- function(TOM, n) {
#   tomMatrix <- matrix(0, n, n)
#   tomMatrix[lower.tri(tomMatrix)] <- TOM
#   tomMatrix <- tomMatrix + t(tomMatrix)
#   tomMatrix
# }
#
# tomMatrix <- listToMatrix(TOM, nGenes)
# rm(TOM)
#
# # Choose which modules to use
# selectedModules <- gsub("^ME", "", names(MEs_DNAme))  # or set manually: c("blue","pink",...)
#
# allHubCpGs <- character(0)
#
# for (module in selectedModules) {
#   probes_idx <- which(moduleColors_DNAme == module)
#   if (length(probes_idx) < 2) next
#
#   modTOM <- tomMatrix[probes_idx, probes_idx, drop = FALSE]
#   cpg_ids <- names(net_PFC_DNAme$colors)[moduleColors_DNAme == module]
#   dimnames(modTOM) <- list(cpg_ids, cpg_ids)
#
#   TOMsimilaritySum <- rowSums(modTOM)
#   hubCpGs <- cpg_ids[TOMsimilaritySum >= quantile(TOMsimilaritySum, 0.9, na.rm = TRUE)]
#   allHubCpGs <- c(allHubCpGs, hubCpGs)
# }
#
# # Enrichment on hub CpGs
# gst_hub_go <- gometh(sig.cpg = allHubCpGs, all.cpg = colnames(PFCSamples), array.type = "EPIC")
# gst_hub_go[order(gst_hub_go$FDR), c("TERM", "FDR")]
#
# gst_hub_kegg <- gometh(sig.cpg = allHubCpGs, all.cpg = colnames(PFCSamples), array.type = "EPIC", collection = "KEGG")
# gst_hub_kegg[order(gst_hub_kegg$FDR), c("Description", "FDR")]
#
# rm(tomMatrix)

# -----------------------------------------------------------------------------
# STEP 12 (Optional): check coefficient in a CpG of interest
# -----------------------------------------------------------------------------
# This uses DNAm_APOE
# summary(lm(DNAm_APOE[, "cg06108056"] ~ niareagansc_recode + Female + Age + pmi + neurons_pfc, data = big_pheno_APOE))
