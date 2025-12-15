# =============================================================================
# 40_multiomic_integration_and_final_analyses.R
#
# Goals:
#   A) Identify DNAme APOE-interaction CpGs tied to protein-gene list
#   B) Link CpG(s) of interest to proteins with APOE moderation
#   C) Relate DNAme module eigengenes to individual proteins and protein-group PC1s
#   D) Compare predictive power (Protein PCs vs DNAme eigengenes vs Combined)
#      stratified by APOE4
#   E) Module↔protein-group overlap enrichment (hypergeometric + missMethyl gsameth)
#
# Expected objects already in environment:
#   - significant_interactions (data.frame with Protein column like "X--GENE")
#   - combined_results (DNAme interaction results with Combined_Names, etc.)
#   - pheno_PFC (must have projid, CleanedSampleID, Female, Age, pmi, neurons_pfc, apoe_genotype, niareagansc_recode)
#   - DNAm_all (matrix/data.frame with rownames = projid, colnames = CpGs)
#   - proteins_PFC_scaled (matrix/data.frame with rownames = projid OR CleanedSampleID)
#   - datExpr (DNAme matrix for WGCNA eigengenes with rownames = projid; columns = CpGs used in WGCNA)
#   - moduleColors_DNAme (vector of module colors corresponding to columns of datExpr)
#   - net_PFC_DNAme (WGCNA network object; net_PFC_DNAme$colors for probe->module)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(broom)
  library(WGCNA)
  library(caret)
  library(pROC)
  library(ggplot2)
  library(minfi)
  library(missMethyl)
  library(org.Hs.eg.db)
})

# -----------------------------------------------------------------------------
# 0) Helpers
# -----------------------------------------------------------------------------
extract_gene <- function(colname) {
  if (!grepl("--", colname)) return(colname)
  strsplit(colname, "--", fixed = TRUE)[[1]][2]
}

sanitize_gene_vec <- function(x) {
  x <- as.character(x)
  x <- unlist(strsplit(x, ";", fixed = TRUE))
  x <- trimws(x)
  x <- x[nzchar(x)]
  unique(x)
}

ensure_ids_are_projid <- function() {
  stopifnot("projid" %in% colnames(pheno_PFC))
  rownames(pheno_PFC) <- pheno_PFC$projid

  if (!all(rownames(DNAm_all) %in% pheno_PFC$projid)) {
    message("Note: DNAm_all has projids not all present in pheno_PFC (OK).")
  }

  # proteins: if rownames are CleanedSampleID, remap to projid via pheno_PFC
  if (!all(rownames(proteins_PFC_scaled) %in% pheno_PFC$projid)) {
    if ("CleanedSampleID" %in% colnames(pheno_PFC)) {
      idx <- match(rownames(proteins_PFC_scaled), pheno_PFC$CleanedSampleID)
      if (sum(is.na(idx)) < length(idx)) {
        message("Remapping proteins_PFC_scaled rownames from CleanedSampleID -> projid")
        rownames(proteins_PFC_scaled) <- pheno_PFC$projid[idx]
      }
    }
  }
}

get_ad_factor <- function(df) {
  ad_bin <- as.integer(as.character(df$niareagansc_recode))
  if (any(is.na(ad_bin))) ad_bin <- as.integer(df$niareagansc_recode) - 1
  factor(ifelse(ad_bin == 1, "AD", "Control"), levels = c("Control", "AD"))
}

# -----------------------------------------------------------------------------
# 1) Gene list and DNAme interaction filtering
# -----------------------------------------------------------------------------
gene_list <- c(
  "GRM2","GSTK1","GLUD1","SYN3","LRRC47","SFXN1","CAP2","HPRT1",
  "PI4KA","DMXL2","ARFGEF3","KIF5C","OXR1","DPP3","VAMP1","FARSB",
  "CASKIN1","SARS2","GNAI2","HSPA8","VPS26B","VGF","SLC25A22","NRXN3","ME3",
  "GNAI3","NCOA7","GNAO1","AHNAK","FGG","HEBP1","APEX1","RAB4A","SLC12A5",
  "LRP1","BAG6","ARPC4","WDR1","EEF2","ALDH7A1","ALDOC","EIF4B","SSBP1",
  "NADK2","KBTBD11","IDH1","CRKL","RYR2","GRIPAP1","BSN","APEH"
)

# Filter combined_results to CpG rows annotated to any gene in gene_list
filtered_results <- combined_results %>%
  filter(sapply(strsplit(Combined_Names, ";", fixed = TRUE),
                function(gs) any(trimws(gs) %in% gene_list)))

# Optional: promoter-only subset
filtered_results_prom <- filtered_results %>%
  filter(Regulatory_Feature_Group %in% c("Promoter_Associated_Cell_type_specific",
                                        "Promoter_Associated"))

filtered_results_prom$interaction_p_adj <- p.adjust(
  filtered_results_prom$APOE_Interaction_Significance, method = "BH"
)

# Within-gene BH
gene_wise_results <- filtered_results_prom %>%
  separate_rows(Combined_Names, sep = ";") %>%
  mutate(Combined_Names = trimws(Combined_Names)) %>%
  filter(Combined_Names %in% gene_list) %>%
  group_by(Combined_Names) %>%
  mutate(gene_specific_p_adj = p.adjust(APOE_Interaction_Significance, method = "BH")) %>%
  ungroup()

# -----------------------------------------------------------------------------
# 2) Align IDs across omics (canonical = projid)
# -----------------------------------------------------------------------------
ensure_ids_are_projid()

common_projid <- Reduce(intersect, list(
  pheno_PFC$projid,
  rownames(DNAm_all),
  rownames(proteins_PFC_scaled),
  rownames(datExpr)
))

message("Common subjects across DNAm/proteomics/datExpr/pheno: ", length(common_projid))

pheno_common <- pheno_PFC[common_projid, ]
DNAm_common  <- DNAm_all[common_projid, , drop = FALSE]
prot_common  <- proteins_PFC_scaled[common_projid, , drop = FALSE]

# Create APOE4 and AD factor
pheno_common$APOE4 <- as.integer(as.character(pheno_common$apoe_genotype))
if (any(is.na(pheno_common$APOE4))) pheno_common$APOE4 <- as.integer(pheno_common$apoe_genotype) - 1
pheno_common$AD_status <- get_ad_factor(pheno_common)

# -----------------------------------------------------------------------------
# 3) Example CpG↔protein relationship with APOE moderation 
# -----------------------------------------------------------------------------
prot_col <- "ELAV4--ELAVL4"
cpg_col  <- "cg06329447"

if (prot_col %in% colnames(prot_common) && cpg_col %in% colnames(DNAm_common)) {
  df_cpg_prot <- data.frame(
    Protein = as.numeric(prot_common[, prot_col]),
    CpG = as.numeric(DNAm_common[, cpg_col]),
    APOE4 = pheno_common$APOE4,
    Female = pheno_common$Female,
    Age = pheno_common$Age,
    neurons_pfc = pheno_common$neurons_pfc,
    pmi = pheno_common$pmi
  ) %>% filter(complete.cases(.))

  fit_cpg_prot <- lm(Protein ~ CpG * APOE4 + Female + Age + neurons_pfc + pmi, data = df_cpg_prot)
  print(summary(fit_cpg_prot))
} else {
  message("Skipping CpG↔protein test: missing ", prot_col, " or ", cpg_col)
}

# -----------------------------------------------------------------------------
# 4) DNAme module eigengenes and protein associations
# -----------------------------------------------------------------------------
# Compute MEs from datExpr restricted to common_projid
datExpr_common <- datExpr[common_projid, , drop = FALSE]
MEs0 <- moduleEigengenes(datExpr_common, colors = moduleColors_DNAme)$eigengenes
MEs_DNAme_common <- orderMEs(MEs0)

# If you want only a subset, keep those that exist:
keep_MEs <- intersect(colnames(MEs_DNAme_common),
                      c("MElightgreen","MEsalmon","MEmidnightblue","MEpink","MEblue","MEtan"))
MEs_DNAme_common <- MEs_DNAme_common[, keep_MEs, drop = FALSE]

# Define interesting proteins as gene symbols, then map to columns
interesting_proteins <- gene_list

prot_genes <- sapply(colnames(prot_common), extract_gene)
selected_cols <- which(prot_genes %in% interesting_proteins)
selected_proteins <- prot_common[, selected_cols, drop = FALSE]
colnames(selected_proteins) <- prot_genes[selected_cols]

# 4A) ME ~ protein (per-protein regressions)
results_list <- list()
for (prot in intersect(interesting_proteins, colnames(selected_proteins))) {
  for (mod in colnames(MEs_DNAme_common)) {
    df <- data.frame(
      ME = as.numeric(MEs_DNAme_common[, mod]),
      protein = as.numeric(selected_proteins[, prot]),
      Female = pheno_common$Female,
      Age = pheno_common$Age,
      pmi = pheno_common$pmi,
      neurons_pfc = pheno_common$neurons_pfc
    ) %>% filter(complete.cases(.))

    fit <- lm(ME ~ protein + Female + Age + pmi + neurons_pfc, data = df)
    tt <- broom::tidy(fit)
    row <- tt[tt$term == "protein", ]
    if (nrow(row) == 1) {
      results_list[[paste(prot, mod, sep = "_")]] <- data.frame(
        Protein = prot, Module = mod,
        Estimate = row$estimate, StdError = row$std.error,
        statistic = row$statistic, p_value = row$p.value
      )
    }
  }
}

results_df <- bind_rows(results_list) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)

# -----------------------------------------------------------------------------
# 5) Protein-group PC1 ↔ DNAme modules
# -----------------------------------------------------------------------------
if (!exists("groups_list")) {
  groups_list <- list(
    Cell_activation_AND_Chemical_synaptic_transmission = c("LRP1","VGF","FGG","BAG6","HEBP1","SLC12A5","BSN"),
    Negative_regulation_of_adenylate_cyclase = c("GNAO1","GRM2","GNAI2","CASKIN1"),
    MET_receptor_recycling_AND_Lipid_metabolism = c("PI4KA","CRKL","RAB4A","GRIPAP1"),
    TCA_cycle = c("GLUD1","DLST","ME3"),
    Protein_translation = c("FARSB","SARS2"),
    Vesicle_mediated_transport = c("SYN3","VAMP1"),
    Protein_folding = c("HSPA8","DPP3"),
    misc = setdiff(gene_list, unique(unlist(.)))
  )
}

prot_genes_common <- sapply(colnames(prot_common), extract_gene)

protein_PC1_df <- data.frame(projid = common_projid)
for (g in names(groups_list)) {
  genes <- groups_list[[g]]
  cols <- which(prot_genes_common %in% genes)
  if (length(cols) < 2) next
  X <- prot_common[, cols, drop = FALSE]
  pca <- prcomp(X, center = TRUE, scale. = TRUE)
  protein_PC1_df[[g]] <- pca$x[,1]
}

# PC1 ~ ME
pc_me_list <- list()
for (g in setdiff(names(protein_PC1_df), "projid")) {
  for (mod in colnames(MEs_DNAme_common)) {
    df <- data.frame(
      ME = as.numeric(MEs_DNAme_common[, mod]),
      PC1 = as.numeric(protein_PC1_df[[g]]),
      Female = pheno_common$Female,
      Age = pheno_common$Age,
      pmi = pheno_common$pmi,
      neurons_pfc = pheno_common$neurons_pfc
    ) %>% filter(complete.cases(.))

    fit <- lm(ME ~ PC1 + Female + Age + pmi + neurons_pfc, data = df)
    tt <- broom::tidy(fit)
    row <- tt[tt$term == "PC1", ]
    if (nrow(row) == 1) {
      pc_me_list[[paste(g, mod, sep = "_")]] <- data.frame(
        Group = g, Module = mod,
        Estimate = row$estimate, StdError = row$std.error,
        statistic = row$statistic, p_value = row$p.value
      )
    }
  }
}

all_results_df <- bind_rows(pc_me_list) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)

# -----------------------------------------------------------------------------
# 6) Predictive power comparison (Protein PCs vs DNAme eigengenes vs Combined)
#    Stratified by APOE4
# -----------------------------------------------------------------------------
ctrl <- trainControl(
  method = "repeatedcv", number = 10, repeats = 10,
  summaryFunction = twoClassSummary, classProbs = TRUE,
  savePredictions = "final"
)

train_glm_cv <- function(df, outcome_col = "AD_status") {
  df <- df %>% filter(complete.cases(.))
  df[[outcome_col]] <- factor(df[[outcome_col]], levels = c("Control","AD"))
  train(
    as.formula(paste(outcome_col, "~ .")),
    data = df,
    method = "glm",
    family = binomial(),
    metric = "ROC",
    trControl = ctrl
  )
}

# Build predictor matrices aligned
PC_df <- protein_PC1_df %>% column_to_rownames("projid")
ME_df <- as.data.frame(MEs_DNAme_common)

# Make sure rownames align
PC_df <- PC_df[common_projid, , drop = FALSE]
ME_df <- ME_df[common_projid, , drop = FALSE]

cov_df <- pheno_common %>%
  dplyr::select(Female, Age, pmi, neurons_pfc)

# Model datasets
df_protein   <- data.frame(AD_status = pheno_common$AD_status, PC_df, cov_df)
df_dname     <- data.frame(AD_status = pheno_common$AD_status, ME_df, cov_df)
df_combined  <- data.frame(AD_status = pheno_common$AD_status, PC_df, ME_df, cov_df)

# Split indices
carriers_idx    <- pheno_common$APOE4 == 1
noncarriers_idx <- pheno_common$APOE4 == 0

set.seed(100)
# Non-carriers
m_prot_nc <- train_glm_cv(df_protein[noncarriers_idx, , drop = FALSE])
m_me_nc   <- train_glm_cv(df_dname[noncarriers_idx, , drop = FALSE])
m_com_nc  <- train_glm_cv(df_combined[noncarriers_idx, , drop = FALSE])

# Carriers
m_prot_c <- train_glm_cv(df_protein[carriers_idx, , drop = FALSE])
m_me_c   <- train_glm_cv(df_dname[carriers_idx, , drop = FALSE])
m_com_c  <- train_glm_cv(df_combined[carriers_idx, , drop = FALSE])

# All
m_prot_all <- train_glm_cv(df_protein)
m_me_all   <- train_glm_cv(df_dname)
m_com_all  <- train_glm_cv(df_combined)

# ROC plotting helper 
get_prob_col <- function(pred_df, positive = "AD") {
  if (positive %in% names(pred_df)) return(positive)
  levs <- unique(pred_df$obs)
  # fallback: last level column
  tail(levs, 1)
}

plot_roc_three <- function(m1, m2, m3, title) {
  p1col <- get_prob_col(m1$pred, "AD")
  p2col <- get_prob_col(m2$pred, "AD")
  p3col <- get_prob_col(m3$pred, "AD")

  r1 <- roc(m1$pred$obs, m1$pred[[p1col]], levels = c("Control","AD"), quiet = TRUE)
  r2 <- roc(m2$pred$obs, m2$pred[[p2col]], levels = c("Control","AD"), quiet = TRUE)
  r3 <- roc(m3$pred$obs, m3$pred[[p3col]], levels = c("Control","AD"), quiet = TRUE)

  plot(r1, col="red", main=title, legacy.axes=TRUE)
  plot(r2, col="blue", add=TRUE)
  plot(r3, col="green", add=TRUE)
  legend("bottomright",
         legend = c(
           paste0("Protein PCs AUC=", round(auc(r1), 2)),
           paste0("DNAme MEs AUC=", round(auc(r2), 2)),
           paste0("Combined AUC=", round(auc(r3), 2))
         ),
         col=c("red","blue","green"), lwd=2, bty="n")
}

par(mfrow=c(3,1))
plot_roc_three(m_prot_nc, m_me_nc, m_com_nc, "Non-carriers")
plot_roc_three(m_prot_c,  m_me_c,  m_com_c,  "Carriers")
plot_roc_three(m_prot_all,m_me_all,m_com_all,"All subjects")

# -----------------------------------------------------------------------------
# 7) Module gene enrichment for protein groups
# -----------------------------------------------------------------------------
# Probe list per module from the network colors
module_probes_list <- split(names(net_PFC_DNAme$colors), moduleColors_DNAme)

# EPIC annotation (hg19)
probe_annotation <- getAnnotation("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")

# Build module genes list
module_genes_list <- lapply(module_probes_list, function(probes) {
  genes <- probe_annotation$UCSC_RefGene_Name[match(probes, probe_annotation$Name)]
  sanitize_gene_vec(genes)
})

# Universe: all genes on the array
all_genes_split <- sanitize_gene_vec(probe_annotation$UCSC_RefGene_Name)

# Hypergeometric enrichment
enrichment_results <- bind_rows(lapply(names(module_genes_list), function(mod) {
  module_genes <- module_genes_list[[mod]]
  bind_rows(lapply(names(groups_list), function(pg) {
    group_genes <- groups_list[[pg]]
    overlap <- length(intersect(module_genes, group_genes))
    m <- length(group_genes)
    n <- length(all_genes_split) - m
    N <- length(module_genes)
    k <- overlap
    p_val <- phyper(k - 1, m, n, N, lower.tail = FALSE)
    data.frame(Module = mod, ProteinGroup = pg, overlap = overlap, p_value = p_val)
  }))
})) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)

# missMethyl gsameth enrichment (accounts for CpG number bias)
groups_list_entrez <- lapply(groups_list, function(genes) {
  entrez <- mapIds(org.Hs.eg.db, keys = genes,
                   column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  unique(na.omit(entrez))
})

collection <- groups_list_entrez
modules_of_interest <- intersect(names(module_probes_list),
                                 c("lightgreen","blue","tan","salmon","midnightblue","pink"))

results_enrichment <- list()
for (mod in modules_of_interest) {
  sig.cpg <- module_probes_list[[mod]]
  res <- gsameth(
    sig.cpg = sig.cpg,
    collection = collection,
    array.type = "EPIC",
    plot.bias = FALSE,
    prior.prob = FALSE,
    equiv.cpg = FALSE,
    fract.counts = FALSE,
    sig.genes = TRUE
  )
  results_enrichment[[mod]] <- res
}
