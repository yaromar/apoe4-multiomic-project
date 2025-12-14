# =============================================================================
# 11_proteomics_normalization_and_integration.R
#
# PURPOSE
# -------
# Starting from:
#   - peptide_intensities_df              (sample x ProteinName table; has SampleID/Batch/Region)
#   - rt_calibration_intensities_df       (RTC wide table; has Sample + RT peptide intensities + metadata)
#   - sample_metadata_df (optional)
#   - big_pheno                           (phenotype table keyed by projid)
#   - ID_mapping.csv                      (maps tube labels -> projid)
#
# This script:
#   1) Keeps only proteomics samples that have RTC calibration.
#   2) Optionally computes RT-based normalization factors.
#   3) Collapses protein entries to gene symbols extracted from "GN=" in ProteinName.
#   4) Splits into regions (PFC/ST/CBM).
#   5) Integrates ID mapping + phenotype (projid).
#   6) Saves cleaned data objects.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
})

setwd("~/Desktop/research/ROSMAP/proteomics/")


# -----------------------------------------------------------------------------
# 1) Drop known problematic samples 
# -----------------------------------------------------------------------------
# - 338P appears in both batch 2 and 3; RTC only batch 2.
# - 234P_r in batch 3 and 234P in batch 6; RTC only batch 6.
#
bad_samples <- c("338P", "234P_r")

peptide_intensities_df_grouped <- peptide_intensities_df %>%
  filter(!(SampleID %in% bad_samples))

# -----------------------------------------------------------------------------
# 2) Keep only samples that have RT calibration entries
# -----------------------------------------------------------------------------
# RTC table uses column 'Sample' as the sample code (e.g., "68P" / "195C" etc).
stopifnot("Sample" %in% colnames(rt_calibration_intensities_df))

peptide_intensities_df_grouped <- peptide_intensities_df_grouped %>%
  filter(SampleID %in% rt_calibration_intensities_df$Sample)

rownames(peptide_intensities_df_grouped) <- peptide_intensities_df_grouped$SampleID

# Recompute Region
peptide_intensities_df_grouped$Region <- ifelse(grepl("C", peptide_intensities_df_grouped$SampleID), "CBM",
                                        ifelse(grepl("P", peptide_intensities_df_grouped$SampleID), "PFC",
                                        ifelse(grepl("S", peptide_intensities_df_grouped$SampleID), "ST", NA)))

# Sample metadata table 
sample_info <- data.frame(
  SampleID = peptide_intensities_df_grouped$SampleID,
  Batch    = peptide_intensities_df_grouped$Batch,
  Region   = peptide_intensities_df_grouped$Region,
  stringsAsFactors = FALSE
)
rownames(sample_info) <- sample_info$SampleID
sample_info <- sample_info[peptide_intensities_df_grouped$SampleID, ]

# -----------------------------------------------------------------------------
# 3) Match RTC calibration rows to proteomics samples (aligned order)
# -----------------------------------------------------------------------------
temp_rt_calibration_intensities_df <- rt_calibration_intensities_df[
  match(peptide_intensities_df_grouped$SampleID, rt_calibration_intensities_df$Sample),
]

stopifnot(all(temp_rt_calibration_intensities_df$Sample == peptide_intensities_df_grouped$SampleID))

# Replace zeros with NA 
peptide_intensities_df_grouped[peptide_intensities_df_grouped == 0] <- NA

# -----------------------------------------------------------------------------
# 4) Define control peptides for normalization 
# -----------------------------------------------------------------------------
control_peptide_columns <- c(
  "Quant..Intensity_ELASGLSFPVGFK",
  "Quant..Intensity_GISNEGQNASIK",
  "Quant..Intensity_IGDYAGIK",
  "Quant..Intensity_LTILEELR",
  "Quant..Intensity_SAAGAFGPELSR"
)

# Keep only the ones that actually exist in the RTC table
control_peptide_columns <- intersect(control_peptide_columns, colnames(temp_rt_calibration_intensities_df))

if (length(control_peptide_columns) == 0) {
  warning("No control peptide columns found in temp_rt_calibration_intensities_df; skipping RT-based normalization.")
}

# Coerce to numeric 
temp_rt_calibration_intensities_df[control_peptide_columns] <- lapply(
  temp_rt_calibration_intensities_df[control_peptide_columns],
  function(x) as.numeric(as.character(x))
)

# Compute control sums 
if (length(control_peptide_columns) > 0) {
  temp_rt_calibration_intensities_df$control_sum <- rowSums(
    temp_rt_calibration_intensities_df[, control_peptide_columns],
    na.rm = TRUE
  )
} else {
  temp_rt_calibration_intensities_df$control_sum <- NA_real_
}

# -----------------------------------------------------------------------------
# 5) OPTIONAL: apply RT normalization (default = OFF )
# -----------------------------------------------------------------------------
apply_rt_normalization <- FALSE  # set TRUE when you want to use RTC scaling

peptide_intensities_df_grouped_norm <- peptide_intensities_df_grouped

normalized_columns <- setdiff(
  names(peptide_intensities_df_grouped_norm),
  c("SampleID", "Batch", "Region")
)

if (apply_rt_normalization && length(control_peptide_columns) > 0) {
  # Example: divide by control_sum 
  peptide_intensities_df_grouped_norm[normalized_columns] <-
    peptide_intensities_df_grouped_norm[normalized_columns] / temp_rt_calibration_intensities_df$control_sum
}

# -----------------------------------------------------------------------------
# 6) Collapse ProteinName -> gene symbol using GN=... extraction 
# -----------------------------------------------------------------------------
# Remove metadata columns before pivoting
expr_only <- peptide_intensities_df_grouped_norm %>%
  select(-SampleID, -Batch, -Region)

# Add row_id to preserve sample order through long->wide transforms
expr_only$row_id <- seq_len(nrow(expr_only))

# STRICT vs PERMISSIVE:
# - strict: drop proteins without GN= 
# - permissive: keep them by assigning a fallback ID (e.g. "UNMAPPED_<protein>")
strict_gn_only <- TRUE

collapsed_to_gene <- expr_only %>%
  pivot_longer(cols = -row_id, names_to = "ProteinName", values_to = "value") %>%
  mutate(
    Gene = str_extract(ProteinName, "GN=([^ ]+)")
  ) %>%
  mutate(
    Gene = sub("^GN=", "", Gene),
    Gene = ifelse(is.na(Gene) & !strict_gn_only,
                  paste0("UNMAPPED_", ProteinName),
                  Gene)
  ) %>%
  filter(!is.na(Gene)) %>%
  group_by(row_id, Gene) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Gene, values_from = value, id_cols = row_id) %>%
  select(-row_id)

# Turn summed zeros into NA 
collapsed_to_gene[collapsed_to_gene == 0] <- NA

# Reattach metadata
peptide_intensities_df_grouped_norm <- collapsed_to_gene
peptide_intensities_df_grouped_norm$SampleID <- sample_info$SampleID
peptide_intensities_df_grouped_norm$Batch    <- sample_info$Batch
peptide_intensities_df_grouped_norm$Region   <- sample_info$Region

# Rowname assignment 
rownames(peptide_intensities_df_grouped_norm) <- peptide_intensities_df_grouped_norm$SampleID

# -----------------------------------------------------------------------------
# 7) Split by brain region
# -----------------------------------------------------------------------------
proteins_CBM <- peptide_intensities_df_grouped_norm[grepl("C", sample_info$SampleID), ]
proteins_PFC <- peptide_intensities_df_grouped_norm[grepl("P", sample_info$SampleID), ]
proteins_ST  <- peptide_intensities_df_grouped_norm[grepl("S", sample_info$SampleID), ]

# -----------------------------------------------------------------------------
# 8) Integrate ID mapping + phenotype (projid)
# -----------------------------------------------------------------------------
ID_mapping <- read.csv("~/Desktop/research/ROSMAP/proteomics/ID_mapping.csv", stringsAsFactors = FALSE)

# Clean tube labels 
sample_info2 <- sample_info %>%
  mutate(CleanedSampleID = gsub("_r.*$", "", SampleID))

# Merge mapping (NOTE: merge reorders; we re-align after)
sample_info2 <- merge(
  sample_info2,
  ID_mapping,
  by.x = "CleanedSampleID",
  by.y = "Tubel.ID.labels",
  all.x = TRUE
)

# Restore rownames and original ordering
sample_info2$SampleID <- as.character(sample_info2$SampleID)
rownames(sample_info2) <- sample_info2$SampleID
sample_info2 <- sample_info2[sample_info$SampleID, , drop = FALSE]

# Ensure projid as character
sample_info2$projid <- as.character(sample_info2$projid)

# ---- pheno_all for all samples (SampleID-aligned)
pheno_all <- merge(
  sample_info2[sample_info2$SampleID %in% rownames(peptide_intensities_df_grouped_norm), ],
  big_pheno,
  by = "projid",
  all.x = TRUE
)
pheno_all <- pheno_all[match(rownames(peptide_intensities_df_grouped_norm), pheno_all$SampleID), ]

# ---- region-specific phenotype merges (projid rownames for protein matrices)
# PFC
rownames(proteins_PFC) <- proteins_PFC$SampleID
pheno_PFC <- merge(
  sample_info2[sample_info2$SampleID %in% rownames(proteins_PFC), ],
  big_pheno,
  by = "projid",
  all.x = TRUE
)
pheno_PFC <- pheno_PFC[match(rownames(proteins_PFC), pheno_PFC$SampleID), ]
rownames(proteins_PFC) <- pheno_PFC$projid

# ST
rownames(proteins_ST) <- proteins_ST$SampleID
pheno_ST <- merge(
  sample_info2[sample_info2$SampleID %in% rownames(proteins_ST), ],
  big_pheno,
  by = "projid",
  all.x = TRUE
)
pheno_ST <- pheno_ST[match(rownames(proteins_ST), pheno_ST$SampleID), ]
rownames(proteins_ST) <- pheno_ST$projid

# CBM
rownames(proteins_CBM) <- proteins_CBM$SampleID
pheno_CBM <- merge(
  sample_info2[sample_info2$SampleID %in% rownames(proteins_CBM), ],
  big_pheno,
  by = "projid",
  all.x = TRUE
)
pheno_CBM <- pheno_CBM[match(rownames(proteins_CBM), pheno_CBM$SampleID), ]
rownames(proteins_CBM) <- pheno_CBM$projid

# -----------------------------------------------------------------------------
# 9) Cleanup + save outputs 
# -----------------------------------------------------------------------------
rm(ID_mapping)

save(
  list = c("proteins_CBM", "proteins_PFC", "proteins_ST",
           "pheno_PFC", "pheno_ST", "pheno_CBM", "sample_info2"),
  file = "cleaned_proteins_not_aggregated.RData"
)

save(
  list = c("peptide_intensities_df_grouped_norm",
           "temp_rt_calibration_intensities_df",
           "control_peptide_columns",
           "pheno_all"),
  file = "peptide_intensities_pre_combination.RData"
)
