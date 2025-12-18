# =============================================================================
# 21_protein_interactions_and_mediation.R
#
# Expected objects in environment:
#   - proteins_PFC: data.frame/matrix [samples x proteins] with columns like "TAU--MAPT"
#   - pheno_PFC: data.frame with rows aligned to proteins_PFC rows (same samples/order),
#       containing: niareagansc_recode, apoe_genotype, Female, Age, pmi, neurons_pfc,
#       CleanedSampleID (optional), Batch (optional), braaksc, ceradsc, apoe_long (optional)
#
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ggeffects)
  library(ggrepel)
  library(ggpubr)
  library(gridExtra)
  library(grid)
  library(MASS)
  library(mediation)
  library(parallel)
  library(broom)
})

# -----------------------------------------------------------------------------
# 0) Sanity checks + phenotype harmonization
# -----------------------------------------------------------------------------
stopifnot(nrow(proteins_PFC) == nrow(pheno_PFC))

# AD (0/1) for binomial glm
pheno_PFC$AD_status <- as.integer(as.character(pheno_PFC$niareagansc_recode))
if (any(is.na(pheno_PFC$AD_status))) {
  pheno_PFC$AD_status <- as.integer(pheno_PFC$niareagansc_recode) - 1
}

# APOE4 carrier (0/1)
pheno_PFC$APOE4 <- as.integer(as.character(pheno_PFC$apoe_genotype))
if (any(is.na(pheno_PFC$APOE4))) {
  pheno_PFC$APOE4 <- as.integer(pheno_PFC$apoe_genotype) - 1
}

# Missing PMI
pheno_PFC$pmi[is.na(pheno_PFC$pmi)] <- mean(pheno_PFC$pmi, na.rm = TRUE)

# -----------------------------------------------------------------------------
# 1) Scale proteins (robust: median/IQR) 
# -----------------------------------------------------------------------------
robust_scale <- function(x) {
  i <- IQR(x, na.rm = TRUE)
  if (is.na(i) || i == 0) return(rep(NA_real_, length(x)))
  (x - median(x, na.rm = TRUE)) / i
}

proteins_PFC_scaled <- as.data.frame(lapply(proteins_PFC, robust_scale))

# -----------------------------------------------------------------------------
# 2) Protein-wise APOE4 × Protein interaction (logistic AD)
# -----------------------------------------------------------------------------

# Robust coefficient extractor by regex
extract_term <- function(coefs, pattern) {
  rn <- rownames(coefs)
  idx <- grep(pattern, rn)
  if (length(idx) == 0) return(c(Estimate=NA, `Std. Error`=NA, `Pr(>|z|)`=NA))
  coefs[idx[1], c("Estimate", "Std. Error", "Pr(>|z|)")]
}

fit_interaction_one_protein <- function(protein) {
  df <- data.frame(
    AD = pheno_PFC$AD_status,
    APOE4 = pheno_PFC$APOE4,
    Protein = as.numeric(proteins_PFC_scaled[[protein]]),
    Female = pheno_PFC$Female,
    Age = pheno_PFC$Age,
    pmi = pheno_PFC$pmi,
    Neurons = pheno_PFC$neurons_pfc
  )
  df <- df[complete.cases(df), ]

  fit <- glm(AD ~ APOE4 * Protein + Age + Female + pmi + Neurons,
             data = df, family = binomial())

  coefs <- summary(fit)$coefficients
  apoe  <- extract_term(coefs, "^APOE4$")
  prot  <- extract_term(coefs, "^Protein$")
  inter <- extract_term(coefs, "APOE4:Protein|Protein:APOE4")

  data.frame(
    Protein = protein,
    Beta_APOE = apoe["Estimate"], SE_APOE = apoe["Std. Error"], P_APOE = apoe["Pr(>|z|)"],
    Beta_Protein = prot["Estimate"], SE_Protein = prot["Std. Error"], P_Protein = prot["Pr(>|z|)"],
    Beta_Interaction = inter["Estimate"], SE_Interaction = inter["Std. Error"], P_Interaction = inter["Pr(>|z|)"],
    N = nrow(df),
    stringsAsFactors = FALSE
  )
}

# Choose proteins to test based on missingness
proteins_to_test <- names(proteins_PFC_scaled)
if (exists("missing_percent")) {
  proteins_to_test <- proteins_to_test[proteins_to_test %in% names(missing_percent)]
  proteins_to_test <- proteins_to_test[missing_percent[proteins_to_test] <= 0.2]
}

interaction_results <- bind_rows(lapply(proteins_to_test, function(p) {
  message("Fitting: ", p)
  tryCatch(fit_interaction_one_protein(p), error = function(e) NULL)
}))

interaction_results <- interaction_results %>%
  mutate(Adj_P_Interaction = p.adjust(P_Interaction, method = "BH")) %>%
  arrange(Adj_P_Interaction)


significant_interactions <- interaction_results %>%
  filter(!is.na(Adj_P_Interaction) & Adj_P_Interaction <= 0.1) 

# -----------------------------------------------------------------------------
# 3) Plots for significant interactions
#    (a) interaction curves (predicted prob AD)
#    (b) adjusted residual boxplots by APOE/AD group
# -----------------------------------------------------------------------------

create_interaction_plot <- function(protein_name) {
  df <- data.frame(
    AD = pheno_PFC$AD_status,
    APOE4 = factor(pheno_PFC$APOE4, levels = c(0,1), labels = c("APOE4-", "APOE4+")),
    Protein = as.numeric(proteins_PFC_scaled[[protein_name]]),
    Female = pheno_PFC$Female,
    Age = pheno_PFC$Age,
    pmi = pheno_PFC$pmi,
    Neurons = pheno_PFC$neurons_pfc
  )
  df <- df[complete.cases(df), ]

  fit <- glm(AD ~ APOE4 * Protein + Age + Female + pmi + Neurons,
             data = df, family = binomial())

  pred <- ggpredict(fit, terms = c("Protein [all]", "APOE4"), type = "fixed")

  coefs <- summary(fit)$coefficients
  inter <- extract_term(coefs, "APOE4\\+?:Protein|Protein:APOE4") # safe-ish
  p_int <- inter["Pr(>|z|)"]

  ggplot(pred, aes(x = x, y = predicted, color = group, fill = group)) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, color = NA) +
    labs(
      title = paste0(protein_name, "  |  interaction p=", ifelse(is.na(p_int), "NA", signif(p_int, 3))),
      x = "Protein (robust-scaled)",
      y = "Predicted P(AD)",
      color = "APOE4",
      fill = "APOE4"
    ) +
    theme_minimal() +
    theme(legend.position = "none", plot.title = element_text(size = 10))
}

create_boxplot_plot <- function(protein_name) {
  df <- data.frame(
    Protein = as.numeric(proteins_PFC_scaled[[protein_name]]),
    APOE4 = factor(pheno_PFC$APOE4, levels = c(0,1), labels = c("APOE4-", "APOE4+")),
    AD = factor(pheno_PFC$AD_status, levels = c(0,1), labels = c("AD-", "AD+")),
    Female = pheno_PFC$Female,
    Age = pheno_PFC$Age,
    pmi = pheno_PFC$pmi,
    Neurons = pheno_PFC$neurons_pfc,
    CleanedSampleID = if ("CleanedSampleID" %in% colnames(pheno_PFC)) pheno_PFC$CleanedSampleID else rownames(pheno_PFC)
  )
  df <- df[complete.cases(df), ]

  df$Group <- factor(paste(df$APOE4, df$AD, sep="/"),
                     levels = c("APOE4-/AD-","APOE4-/AD+","APOE4+/AD-","APOE4+/AD+"))

  # Covariate-adjust residuals
  adj <- lm(Protein ~ Female + Age + pmi + Neurons, data = df)
  df$Resid <- residuals(adj)

  # Outliers by 3*IQR rule 
  df <- df %>%
    group_by(Group) %>%
    mutate(
      Q1 = quantile(Resid, 0.25, na.rm = TRUE),
      Q3 = quantile(Resid, 0.75, na.rm = TRUE),
      IQRv = Q3 - Q1,
      lower = Q1 - 3 * IQRv,
      upper = Q3 + 3 * IQRv,
      is_outlier = Resid < lower | Resid > upper
    ) %>% ungroup()

  outlier_ids <- df$CleanedSampleID[df$is_outlier]

  p <- ggplot(df, aes(x = Group, y = Resid)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 0.6) +
    theme_minimal() +
    labs(title = protein_name, x = NULL, y = "Adjusted residuals") +
    theme(plot.title = element_text(size = 10),
          axis.text.x = element_text(size = 7, hjust = 0.5),
          legend.position = "none") +
    stat_compare_means(
      comparisons = list(c("APOE4-/AD-", "APOE4-/AD+"),
                         c("APOE4+/AD-", "APOE4+/AD+")),
      method = "wilcox.test",
      label = "p.signif"
    )

  list(plot = p, outliers = outlier_ids)
}

sig_proteins <- significant_interactions$Protein
sig_proteins <- sig_proteins[!(sig_proteins %in% exclude_list)]

interaction_plot_list <- list()
boxplot_list <- list()
outlier_tracker <- list()

for (p in sig_proteins) {
  interaction_plot_list[[p]] <- tryCatch(create_interaction_plot(p), error = function(e) NULL)
  bp <- tryCatch(create_boxplot_plot(p), error = function(e) NULL)
  if (!is.null(bp)) {
    boxplot_list[[p]] <- bp$plot
    outlier_tracker[[p]] <- bp$outliers
  }
}

# -----------------------------------------------------------------------------
# 4) Outlier frequency diagnostics (optional)
# -----------------------------------------------------------------------------
if (length(outlier_tracker) > 0 && "CleanedSampleID" %in% colnames(pheno_PFC)) {

  outlier_counts <- table(unlist(outlier_tracker))
  outlier_df <- data.frame(SubjectID = names(outlier_counts),
                           Frequency = as.numeric(outlier_counts)) %>%
    arrange(desc(Frequency))

  outlier_group_df <- outlier_df %>%
    left_join(pheno_PFC, by = c("SubjectID" = "CleanedSampleID")) %>%
    mutate(
      Batch = if ("Batch" %in% colnames(pheno_PFC)) factor(Batch) else NA,
      Female = factor(Female, levels = c(0,1), labels = c("Male","Female"))
    )

  if ("Batch" %in% colnames(outlier_group_df)) {
    print(ggplot(outlier_group_df, aes(x = Batch, y = Frequency)) + geom_boxplot() + theme_minimal())
  }
  print(ggplot(outlier_group_df, aes(x = Age, y = Frequency)) + geom_point(alpha=0.6) + geom_smooth(method="loess") + theme_minimal())
  print(ggplot(outlier_group_df, aes(x = neurons_pfc, y = Frequency)) + geom_point(alpha=0.6) + geom_smooth(method="lm") + theme_minimal())
  print(ggplot(outlier_group_df, aes(x = Female, y = Frequency)) + geom_boxplot() + theme_minimal())
}

# -----------------------------------------------------------------------------
# 5) Save multi-page PDF of interaction+boxplots (grouping optional)
# -----------------------------------------------------------------------------
if (length(interaction_plot_list) > 0 && length(boxplot_list) > 0) {

  # Simple paging: 3 proteins per row, each protein gets 2 panels (interaction + box)
  proteins_per_page <- 6  # 6 proteins = 12 panels
  plotted <- intersect(names(interaction_plot_list), names(boxplot_list))
  plotted <- plotted[!sapply(plotted, function(p) is.null(interaction_plot_list[[p]]) || is.null(boxplot_list[[p]]))]

  if (length(plotted) > 0) {
    pdf("~/Desktop/interaction_plots_by_cluster.pdf", width = 11, height = 8.5)

    pages <- split(plotted, ceiling(seq_along(plotted) / proteins_per_page))
    for (pg in pages) {
      grid.newpage()
      grobs <- list()
      for (p in pg) {
        grobs <- c(grobs, list(interaction_plot_list[[p]], boxplot_list[[p]]))
      }
      grid.draw(arrangeGrob(grobs = grobs, ncol = 3)) # 3 columns of panels
    }

    dev.off()
  }
}

# -----------------------------------------------------------------------------
# 6) Mediation: APOE4 -> Protein -> AD 
# -----------------------------------------------------------------------------
perform_mediation_one_protein <- function(protein_name, sims = 1000) {
  df <- data.frame(
    protein = as.numeric(proteins_PFC_scaled[[protein_name]]),
    APOE4 = pheno_PFC$APOE4,
    AD = pheno_PFC$AD_status,
    Female = pheno_PFC$Female,
    Age = pheno_PFC$Age,
    pmi = pheno_PFC$pmi,
    neurons_pfc = pheno_PFC$neurons_pfc
  )
  df <- df[complete.cases(df), ]

  med_model <- lm(protein ~ APOE4 + Female + Age + pmi + neurons_pfc, data = df)

  # Keep your “moderated” version in outcome model if desired:
  out_model <- glm(AD ~ APOE4 * protein + Female + Age + pmi + neurons_pfc,
                   family = binomial(), data = df)

  mediate(med_model, out_model,
          treat = "APOE4", mediator = "protein",
          boot = TRUE, sims = sims)
}

# Choose proteins for mediation 
proteins_for_mediation <- sig_proteins

run_mediation_parallel <- function(protein_list, sims = 1000, ncores = max(1, detectCores()-2)) {
  cl <- makeCluster(ncores)
  clusterExport(cl, c("proteins_PFC_scaled", "pheno_PFC", "perform_mediation_one_protein", "sims"), envir = environment())
  clusterEvalQ(cl, { library(mediation); library(stats) })
  res <- parLapply(cl, protein_list, function(p) {
    tryCatch(perform_mediation_one_protein(p, sims = sims), error = function(e) NULL)
  })
  stopCluster(cl)
  setNames(res, protein_list)
}

protein_mediation <- run_mediation_parallel(proteins_for_mediation, sims = 1000)

mediation_table <- bind_rows(lapply(names(protein_mediation), function(p) {
  mr <- protein_mediation[[p]]
  if (is.null(mr)) return(NULL)
  data.frame(
    Protein = p,
    ACME = mr$d.avg,
    ACME_p = mr$d.avg.p,
    ADE = mr$z.avg,
    ADE_p = mr$z.avg.p,
    Total_Effect = mr$tau.coef,
    Total_Effect_p = mr$tau.p,
    Prop_Mediated = mr$n.avg,
    Prop_Mediated_p = mr$n.avg.p,
    stringsAsFactors = FALSE
  )
})) %>%
  mutate(ACME_adj = p.adjust(ACME_p, "BH")) %>%
  arrange(ACME_p)

# -----------------------------------------------------------------------------
# 7) Neuropath associations (binary Braak/CERAD)
# -----------------------------------------------------------------------------
pheno_PFC$braaksc_num <- ifelse(pheno_PFC$braaksc %in% c(3,4), 1,
                               ifelse(pheno_PFC$braaksc %in% c(1,2), 0, NA))
pheno_PFC$ceradsc_num <- ifelse(pheno_PFC$ceradsc == 2, 1,
                               ifelse(pheno_PFC$ceradsc == 1, 0, NA))

test_protein_neuro <- function(protein_name, outcome_col, with_interaction = FALSE) {
  df <- data.frame(
    outcome = pheno_PFC[[outcome_col]],
    Protein = as.numeric(proteins_PFC_scaled[[protein_name]]),
    APOE4 = pheno_PFC$APOE4,
    Female = pheno_PFC$Female,
    Age = pheno_PFC$Age,
    pmi = pheno_PFC$pmi,
    Neurons = pheno_PFC$neurons_pfc
  )
  df <- df[complete.cases(df), ]

  if (with_interaction) {
    f <- outcome ~ Protein * APOE4 + Female + Age + pmi + Neurons
  } else {
    f <- outcome ~ Protein + Female + Age + pmi + Neurons
  }

  fit <- glm(f, family = binomial(), data = df)
  broom::tidy(fit)
}

neuro_proteins <- unique(c(sig_proteins,
                           mediation_table$Protein[mediation_table$ACME_adj < 0.1]))

neuro_out <- bind_rows(lapply(neuro_proteins, function(p) {
  br <- tryCatch(test_protein_neuro(p, "braaksc_num", with_interaction = FALSE), error = function(e) NULL)
  ce <- tryCatch(test_protein_neuro(p, "ceradsc_num", with_interaction = FALSE), error = function(e) NULL)
  if (is.null(br) && is.null(ce)) return(NULL)
  bind_rows(
    if (!is.null(br)) mutate(br, ProteinName = p, Outcome = "Braak") else NULL,
    if (!is.null(ce)) mutate(ce, ProteinName = p, Outcome = "CERAD") else NULL
  )
}))

# -----------------------------------------------------------------------------
# 8) Enrichment (clusterProfiler)
# -----------------------------------------------------------------------------
# Example for interaction proteins using gene SYMBOL after "--"

# library(clusterProfiler)
# library(org.Hs.eg.db)
#
# sig_symbols <- sapply(strsplit(sig_proteins, "--"), function(x) if (length(x) > 1) x[2] else NA)
# sig_symbols <- unique(na.omit(sig_symbols))
#
# universe_symbols <- names(proteins_PFC_scaled)
# universe_symbols <- sapply(strsplit(universe_symbols, "--"), function(x) if (length(x) > 1) x[2] else NA)
# universe_symbols <- unique(na.omit(universe_symbols))
#
# ego <- enrichGO(
#   gene = sig_symbols,
#   OrgDb = org.Hs.eg.db,
#   ont = "BP",
#   universe = universe_symbols,
#   pAdjustMethod = "BH",
#   pvalueCutoff = 0.15,
#   keyType = "SYMBOL",
#   readable = TRUE
# )
