# =============================================================================
# 10_proteomics_import_and_cleaning.R
#
# PURPOSE
# -------
# 1) Read Levine DIA proteomics "Samples Table" CSVs across batches 1–6.
#    Build a sample x protein wide matrix of intensities, with Batch label.
#    Do quick missingness and abundance QC; infer brain region from SampleID.
#
# 2) Read Levine RTC "Individual Peptide Result" CSVs across batches 1–6.
#    Parse sample identifiers from mzML-like strings; build peptide-wide table
#    and extract sample metadata (Sample, SampleNumber, Batch_RTC, Region).
#
# EXPECTED INPUT FILE STRUCTURE
# -----------------------------
# Protein-level tables:
#   ./proteomics/ROSMAP Proteomics/Levine DIA Result_Batch-0X/.../*.csv
#    - must contain column 'Protein Name'
#    - intensity is assumed to be the last column
#
# RTC peptide result tables:
#   ./proteomics/ROSMAP Proteomics/Levine_Batch-0X_RTC Individual Peptide Result/*.csv
#    - must contain column 'Sample' with mzML filename string
#    - columns: Sequence, Quant..Intensity, RT.Start..min., RT.Center..min., RT.Stop..min.
#
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(ggplot2)
  library(stringr)
})

setwd("~/Desktop/research/ROSMAP/")

# -----------------------------------------------------------------------------
# PART A: Protein-level DIA "Samples Table" import (sample x protein)
# -----------------------------------------------------------------------------

# ---- Helper: read a single "Samples Table" csv and return a tidy long table
# Output columns: ProteinName, Intensity, SampleID
process_sample_file <- function(file_path) {
  # Read the CSV
  data <- read_csv(file_path, show_col_types = FALSE)

  # Extract sample identifier from filename
  sample_id <- tools::file_path_sans_ext(basename(file_path))
  sample_id <- gsub("Samples Table of Levine DIA ", "", sample_id)

  # Keep only the two needed fields: 'Protein Name' + last column as Intensity
  data <- data %>%
    select(ProteinName = `Protein Name`, Intensity = last_col()) %>%
    mutate(
      Intensity = as.numeric(Intensity),
      SampleID  = sample_id
    )

  data
}

# ---- Helper: import one batch directory and return wide sample x protein table
read_protein_batch <- function(batch_id, batch_dir) {
  files <- list.files(
    path = batch_dir,
    pattern = "\\.csv$",
    full.names = TRUE,
    recursive = TRUE
  )

  # Read all sample files in batch into long form
  combined <- map_df(files, process_sample_file)

  # Remove rows without intensity
  combined <- combined[!is.na(combined$Intensity), ]

  combined_no_agg <- combined %>%
    group_by(SampleID, ProteinName) %>%
    ungroup()

  # Pivot to wide: rows = samples, columns = proteins
  wide <- combined_no_agg %>%
    pivot_wider(
      names_from  = ProteinName,
      values_from = Intensity,
      values_fill = list(Intensity = NA)
    )

  wide$Batch <- as.character(batch_id)
  wide
}

# ---- Batch directories 
batch_dirs <- list(
  `1` = "./proteomics/ROSMAP Proteomics/Levine DIA Result_Batch-01/Samples Table of Levine DIA Batch-01",
  `2` = "./proteomics/ROSMAP Proteomics/Levine DIA Result_Batch-02/Samples Table of Levine DIA Batch-02",
  `3` = "./proteomics/ROSMAP Proteomics/Levine DIA Result_Batch-03/Samples Table of Levine DIA_Batch-03",
  `4` = "./proteomics/ROSMAP Proteomics/Levine DIA Result_Batch-04/Samples Table of Levine DIA_Batch-04",
  `5` = "./proteomics/ROSMAP Proteomics/Levine DIA Result_Batch-05/Samples Table of Levine DIA_Batch-05",
  `6` = "./proteomics/ROSMAP Proteomics/Levine DIA Result_Batch-06/Samples Table of Levine DIA_Batch-06"
)

# ---- Read all batches, bind to one protein intensity matrix
final_by_batch <- imap(batch_dirs, ~ read_protein_batch(batch_id = .y, batch_dir = .x))
peptide_intensities_df <- bind_rows(final_by_batch)

# One-off SampleID cleanup
peptide_intensities_df$SampleID <- gsub("^Samples Table of Levine Dia 254S$", "254S", peptide_intensities_df$SampleID)

# Remove proteins with extreme missingness 
# NOTE: This assumes ~1040 samples
peptide_intensities_df <- peptide_intensities_df[, -which(colSums(is.na(peptide_intensities_df)) >= 1040)]

# Remove fully missing samples
peptide_intensities_df <- peptide_intensities_df %>%
  filter(rowSums(is.na(.)) < ncol(.))

# Infer brain region from SampleID 
peptide_intensities_df$Region <- ifelse(grepl("C", peptide_intensities_df$SampleID), "CBM",
                                ifelse(grepl("P", peptide_intensities_df$SampleID), "PFC",
                                ifelse(grepl("S", peptide_intensities_df$SampleID), "ST", NA)))

# -----------------------------------------------------------------------------
# Protein-level QC (optional)
# -----------------------------------------------------------------------------
missing_per_column <- peptide_intensities_df %>%
  summarise(across(where(is.numeric), ~ sum(is.na(.))))

hist(as.numeric(missing_per_column), main = "Missing values per protein", xlab = "# missing")

missing_per_row_n <- peptide_intensities_df %>%
  mutate(missing_count_n = rowSums(!is.na(.))) %>%
  select(SampleID, missing_count_n)

ggplot(missing_per_row_n, aes(x = missing_count_n)) +
  theme_minimal() +
  geom_histogram(binwidth = 50, fill = "blue", color = "black") +
  labs(
    title = "Histogram of Unique Proteins per Sample",
    x = "Number of Unique Proteins",
    y = "Frequency"
  )

# Total intensity distribution by batch 
total_intensity_all <- rowSums(
  peptide_intensities_df[, !(colnames(peptide_intensities_df) %in% c("SampleID", "Batch", "Region"))],
  na.rm = TRUE
)

hist(
  rowSums(
    peptide_intensities_df[peptide_intensities_df$Batch == "4",
                           !(colnames(peptide_intensities_df) %in% c("SampleID", "Batch", "Region"))],
    na.rm = TRUE
  ),
  breaks = 30,
  xlim = range(total_intensity_all),
  main = "Abundance of proteins -- batch 4",
  xlab = "Total Intensity per Sample"
)

# Abundance by region 
hist(
  rowSums(
    peptide_intensities_df[peptide_intensities_df$Region == "CBM",
                           !(colnames(peptide_intensities_df) %in% c("SampleID", "Batch", "Region"))],
    na.rm = TRUE
  ),
  breaks = 30,
  xlim = range(total_intensity_all),
  main = "Abundance of proteins -- CBM",
  xlab = "Total Intensity per Sample"
)

# -----------------------------------------------------------------------------
# PART B: RTC calibration peptide result import
# -----------------------------------------------------------------------------

# ---- Helper: parse mzML-like sample strings into SampleNumber + Sample code
parse_rtc_sample_fields <- function(sample_str) {
  pattern <- "HFX(\\d+)-(\\d+)_?.*?_(\\d+[CPS](_r\\d?)?)_?.*?.mzML"
  m <- str_match(sample_str, pattern)

  # sample_number = instrument + run 
  sample_number <- paste0(m[, 2], m[, 3])
  sample_code   <- m[, 4]

  list(SampleNumber = sample_number, Sample = sample_code)
}

# ---- Helper: process one RTC peptide CSV
process_rtc_file <- function(file_path) {
  data <- read.csv(file_path)

  # 'Sample' column is expected to contain mzML filename string(s)
  matches <- str_match(data$Sample, "HFX(\\d+)-(\\d+)_?.*?_(\\d+[CPS](_r\\d?)?)_?.*?.mzML")

  sample_number <- apply(matches[, 2:3], 1, function(x) paste0(x, collapse = ""))
  sample_code   <- matches[, 4]

  data$SampleNumber <- sample_number
  data$Sample <- sample_code

  data %>%
    select(Sample, Sequence, Quant..Intensity, RT.Start..min., RT.Center..min., RT.Stop..min., SampleNumber)
}

# ---- Batch RTC directories + batch-specific file slicing 
rtc_batch_dirs <- list(
  `1` = list(dir = "./proteomics/ROSMAP Proteomics/Levine_Batch_01_RTC Individual Peptide Result", slice = 3:16),
  `2` = list(dir = "./proteomics/ROSMAP Proteomics/Levine_Batch-02_RTC Individual Peptide Result", slice = 2:15),
  `3` = list(dir = "./proteomics/ROSMAP Proteomics/Levine_Batch-03_RTC Individual Peptide Result", slice = 2:15),
  `4` = list(dir = "./proteomics/ROSMAP Proteomics/Levine_Batch-04_RTC Individual Peptide Result", slice = 2:15),
  `5` = list(dir = "./proteomics/ROSMAP Proteomics/Levine_Batch-05_RTC Individual Peptide Result", slice = 2:15),
  `6` = list(dir = "./proteomics/ROSMAP Proteomics/Levine_Batch-06_RTC Individual Peptide Result", slice = 2:15)
)

read_rtc_batch <- function(batch_id, info) {
  file_paths <- list.files(info$dir, pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)

  file_paths <- file_paths[info$slice]

  dat <- map_df(file_paths, process_rtc_file)
  dat$Batch_RTC <- as.character(batch_id)
  dat
}

rtc_list <- imap(rtc_batch_dirs, ~ read_rtc_batch(batch_id = .y, info = .x))
rt_calibration_intensities_df <- bind_rows(rtc_list)

# -----------------------------------------------------------------------------
# RTC cleaning before pivot
# -----------------------------------------------------------------------------
# Identify duplicates 
rt_calibration_intensities_df %>%
  group_by(Sample, Sequence) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n > 1)

# Drop NA samples
rt_calibration_intensities_df <- rt_calibration_intensities_df %>% filter(!is.na(Sample))

# Drop "END OF FILE" pseudo-sequence rows
rt_calibration_intensities_df <- rt_calibration_intensities_df %>% filter(Sequence != "END OF FILE")

# Drop exact duplicates
rt_calibration_intensities_df <- rt_calibration_intensities_df %>% distinct()

# Re-check duplicates post-distinct 
rt_calibration_intensities_df %>%
  summarise(n = n(), .by = c(Sample, Sequence)) %>%
  filter(n > 1)

# -----------------------------------------------------------------------------
# Pivot RTC table to wide: one row per Sample, columns per peptide sequence
# -----------------------------------------------------------------------------
rt_calibration_intensities_df <- rt_calibration_intensities_df %>%
  pivot_wider(
    names_from = Sequence,
    values_from = c(Quant..Intensity, RT.Start..min., RT.Center..min., RT.Stop..min., Batch_RTC, SampleNumber),
    names_glue = "{.value}_{Sequence}",
    id_cols = Sample,
    values_fill = list(
      Quant..Intensity = NA,
      RT.Start..min. = NA,
      RT.Center..min. = NA,
      RT.Stop..min. = NA,
      Batch_RTC = NA,
      SampleNumber = NA
    )
  )

# Add Region column inferred from Sample code 
rt_calibration_intensities_df$Region <- ifelse(grepl("C", rt_calibration_intensities_df$Sample), "CBM",
                                       ifelse(grepl("P", rt_calibration_intensities_df$Sample), "PFC",
                                       ifelse(grepl("S", rt_calibration_intensities_df$Sample), "ST", NA)))

# -----------------------------------------------------------------------------
# Extract SampleNumber + Batch from wide table
# -----------------------------------------------------------------------------:
#   SampleNumber <- SampleNumber_SFANQPLEVVYSK
#   Batch <- Batch_RTC_ELASGLSFPVGFK

# Find any SampleNumber_* column and any Batch_RTC_* column to use as representative.
sampleNumber_cols <- grep("^SampleNumber_", names(rt_calibration_intensities_df), value = TRUE)
batch_cols        <- grep("^Batch_RTC_",    names(rt_calibration_intensities_df), value = TRUE)

# Prefer the original specific sequences if present, else take the first available
preferred_sn <- "SampleNumber_SFANQPLEVVYSK"
preferred_b  <- "Batch_RTC_ELASGLSFPVGFK"

sn_col <- if (preferred_sn %in% sampleNumber_cols) preferred_sn else sampleNumber_cols[1]
b_col  <- if (preferred_b  %in% batch_cols)        preferred_b  else batch_cols[1]

rt_calibration_intensities_df$SampleNumber <- rt_calibration_intensities_df[[sn_col]]
rt_calibration_intensities_df$Batch        <- rt_calibration_intensities_df[[b_col]]

# Build sample metadata table
sample_metadata_df <- rt_calibration_intensities_df %>%
  select(Sample, SampleNumber, Batch, Region) %>%
  distinct()

# -----------------------------------------------------------------------------
# Safe cleanup 
# -----------------------------------------------------------------------------
to_rm <- c(
  "rtc_list", "final_by_batch",
  "missing_per_column", "missing_per_row_n", "total_intensity_all"
)
rm(list = intersect(to_rm, ls()))
