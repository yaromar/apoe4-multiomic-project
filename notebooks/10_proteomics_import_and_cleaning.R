# Yaro
# Proteomics + RT calibration import + QC scaffolding
#
# Purpose:
#   1) Read per-sample "Samples Table of Levine DIA ..." CSVs across batches.
#   2) Build a SampleID × ProteinName intensity matrix with Batch info.
#   3) Read RT calibration "Individual Peptide Result" CSVs across batches.
#   4) Parse sample code from mzML-style Sample strings into a clean Sample key.
#   5) Pivot RT calibration to wide (Sample × {Quant/RT/Batch/etc}_PEPTIDESEQ).
#   6) Build sample metadata table (Sample, SampleNumber, Batch, Region).
#
# Notes:
#   - This script currently keeps ProteinName (not collapsing to GeneName).
#   - It does not yet merge proteins with RTC (but sets up consistent keys to do so).

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(stringr)
  library(ggplot2)
})

setwd("~/Desktop/research/ROSMAP/")
set.seed(0)

# -----------------------------
# 1) DIA protein tables import
# -----------------------------

batch_dirs <- list(
  "1" = "./proteomics/ROSMAP Proteomics/Levine DIA Result_Batch-01/Samples Table of Levine DIA Batch-01",
  "2" = "./proteomics/ROSMAP Proteomics/Levine DIA Result_Batch-02/Samples Table of Levine DIA Batch-02",
  "3" = "./proteomics/ROSMAP Proteomics/Levine DIA Result_Batch-03/Samples Table of Levine DIA_Batch-03",
  "4" = "./proteomics/ROSMAP Proteomics/Levine DIA Result_Batch-04/Samples Table of Levine DIA_Batch-04",
  "5" = "./proteomics/ROSMAP Proteomics/Levine DIA Result_Batch-05/Samples Table of Levine DIA Batch-05",
  "6" = "./proteomics/ROSMAP Proteomics/Levine DIA Result_Batch-06/Samples Table of Levine DIA Batch-06"
)

extract_sample_id_from_filename <- function(file_path) {
  sample_id <- tools::file_path_sans_ext(basename(file_path))
  sample_id <- gsub("Samples Table of Levine DIA ", "", sample_id)
  sample_id
}

process_sample_file <- function(file_path, batch) {
  df <- readr::read_csv(file_path, show_col_types = FALSE)

  # Assumes:
  #   - "Protein Name" exists
  #   - last column is the intensity column
  df <- df %>%
    transmute(
      ProteinName = .data[["Protein Name"]],
      Intensity = suppressWarnings(as.numeric(dplyr::last_col())),
      SampleID = extract_sample_id_from_filename(file_path),
      Batch = as.character(batch)
    ) %>%
    filter(!is.na(Intensity), !is.na(ProteinName))

  df
}

protein_long <- imap_dfr(batch_dirs, function(dir_path, batch_id) {
  files <- list.files(dir_path, pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)
  if (length(files) == 0) stop(paste("No CSV files in:", dir_path))

  map_dfr(files, ~process_sample_file(.x, batch = batch_id))
})

# If there are duplicate rows for the same SampleID × ProteinName, collapse by sum
protein_long <- protein_long %>%
  group_by(SampleID, ProteinName, Batch) %>%
  summarise(Intensity = sum(Intensity, na.rm = TRUE), .groups = "drop")

# Fix one known SampleID typo (only if it exists)
protein_long$SampleID[protein_long$SampleID == "Samples Table of Levine Dia 254S"] <- "254S"

# Build Sample × Protein matrix; keep SampleID + Batch as identifiers.
peptide_intensities_df <- protein_long %>%
  pivot_wider(
    id_cols = c(SampleID, Batch),
    names_from = ProteinName,
    values_from = Intensity,
    values_fill = NA
  )

# Remove proteins missing in (almost) everyone
# NOTE: this uses a dynamic threshold rather than a hard-coded 1040.
na_per_col <- colSums(is.na(peptide_intensities_df))
# don’t apply to metadata columns
na_per_col[names(na_per_col) %in% c("SampleID", "Batch")] <- 0

# Example policy: drop proteins present in < 2 samples
keep_cols <- names(na_per_col)[na_per_col < (nrow(peptide_intensities_df) - 1)]
peptide_intensities_df <- peptide_intensities_df %>%
  select(all_of(c("SampleID", "Batch", keep_cols[!keep_cols %in% c("SampleID", "Batch")])))

# Remove completely missing samples (all protein columns NA)
peptide_intensities_df <- peptide_intensities_df %>%
  filter(rowSums(is.na(across(-c(SampleID, Batch)))) < (ncol(.) - 2))

# Add Region derived from SampleID (C/P/S)
peptide_intensities_df <- peptide_intensities_df %>%
  mutate(
    Region = case_when(
      grepl("C", SampleID) ~ "CBM",
      grepl("P", SampleID) ~ "PFC",
      grepl("S", SampleID) ~ "ST",
      TRUE ~ NA_character_
    )
  )

# QC: missingness per protein
missing_per_column <- peptide_intensities_df %>%
  summarise(across(where(is.numeric), ~sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "Protein", values_to = "MissingCount")

ggplot(missing_per_column, aes(x = MissingCount)) +
  geom_histogram(bins = 50) +
  theme_minimal() +
  labs(title = "Missing values per protein column", x = "# missing", y = "Count of proteins")

# QC: number of non-missing proteins per sample
nonmissing_per_sample <- peptide_intensities_df %>%
  transmute(SampleID, NonMissingProteins = rowSums(!is.na(across(where(is.numeric)))))

ggplot(nonmissing_per_sample, aes(x = NonMissingProteins)) +
  geom_histogram(binwidth = 50) +
  theme_minimal() +
  labs(title = "Unique proteins per sample", x = "# non-missing proteins", y = "Frequency")

# -----------------------------------------
# 2) RT calibration peptide results import
# -----------------------------------------

rtc_dirs <- list(
  "1" = "./proteomics/ROSMAP Proteomics/Levine_Batch_01_RTC Individual Peptide Result",
  "2" = "./proteomics/ROSMAP Proteomics/Levine_Batch-02_RTC Individual Peptide Result",
  "3" = "./proteomics/ROSMAP Proteomics/Levine_Batch-03_RTC Individual Peptide Result",
  "4" = "./proteomics/ROSMAP Proteomics/Levine_Batch-04_RTC Individual Peptide Result",
  "5" = "./proteomics/ROSMAP Proteomics/Levine_Batch-05_RTC Individual Peptide Result",
  "6" = "./proteomics/ROSMAP Proteomics/Levine_Batch-06_RTC Individual Peptide Result"
)

parse_rtc_sample_fields <- function(sample_vec) {
  # Expected sample string looks like:
  #   HFX20-2876_Levine_68P_r_20201104183859.mzML
  # or HFX19-1468_Levine_195C.mzML
  #
  # We want:
  #   InstrumentNumber: 20
  #   SampleNumber: 2876  (or concatenated 20 + 2876 if you really want that)
  #   SampleCode: 68P / 195C / etc (including _r or _r1)
  pat <- "HFX(\\d+)-(\\d+)_?.*?_(\\d+[CPS](?:_r\\d?)?)_?.*?\\.mzML"

  m <- stringr::str_match(sample_vec, pat)
  tibble(
    instrument = m[, 2],
    run_number = m[, 3],
    Sample = m[, 4]
  )
}

process_rtc_file <- function(file_path, batch) {
  df <- read.csv(file_path)

  # Require expected columns
  stopifnot(all(c("Sample", "Sequence", "Quant..Intensity",
                  "RT.Start..min.", "RT.Center..min.", "RT.Stop..min.") %in% names(df)))

  parsed <- parse_rtc_sample_fields(df$Sample)

  df <- df %>%
    mutate(
      Sample = parsed$Sample,
      SampleNumber = paste0(parsed$instrument, parsed$run_number),
      Batch_RTC = as.character(batch)
    ) %>%
    filter(!is.na(Sample), Sequence != "END OF FILE") %>%
    select(Sample, Sequence, Quant..Intensity, RT.Start..min., RT.Center..min., RT.Stop..min., SampleNumber, Batch_RTC) %>%
    distinct()

  df
}

rt_calibration_long <- imap_dfr(rtc_dirs, function(dir_path, batch_id) {
  files <- list.files(dir_path, pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)
  if (length(files) == 0) stop(paste("No CSV files in:", dir_path))

  # IMPORTANT: do NOT slice files[2:15] etc. unless you have a known reason.
  map_dfr(files, ~process_rtc_file(.x, batch = batch_id))
})

# Sanity: check duplicates (should be 1 row per Sample×Sequence after distinct)
dup_check <- rt_calibration_long %>%
  count(Sample, Sequence) %>%
  filter(n > 1)

if (nrow(dup_check) > 0) {
  message("WARNING: duplicates remain for some Sample×Sequence pairs.")
  print(dup_check)
}

# Pivot RTC to wide: each peptide sequence becomes a set of columns
rt_calibration_wide <- rt_calibration_long %>%
  pivot_wider(
    id_cols = Sample,
    names_from = Sequence,
    values_from = c(Quant..Intensity, RT.Start..min., RT.Center..min., RT.Stop..min., Batch_RTC, SampleNumber),
    names_glue = "{.value}_{Sequence}",
    values_fill = NA
  )

# Add Region derived from Sample
rt_calibration_wide <- rt_calibration_wide %>%
  mutate(
    Region = case_when(
      grepl("C", Sample) ~ "CBM",
      grepl("P", Sample) ~ "PFC",
      grepl("S", Sample) ~ "ST",
      TRUE ~ NA_character_
    )
  )

# Pick one representative SampleNumber and Batch column (choose a peptide you KNOW is always present)
# Here we use ELASGLSFPVGFK and SFANQPLEVVYSK only if they exist.
if ("SampleNumber_SFANQPLEVVYSK" %in% names(rt_calibration_wide)) {
  rt_calibration_wide <- rt_calibration_wide %>%
    mutate(SampleNumber = .data$SampleNumber_SFANQPLEVVYSK)
}

if ("Batch_RTC_ELASGLSFPVGFK" %in% names(rt_calibration_wide)) {
  rt_calibration_wide <- rt_calibration_wide %>%
    mutate(Batch = .data$Batch_RTC_ELASGLSFPVGFK)
}

sample_metadata_df <- rt_calibration_wide %>%
  select(Sample, SampleNumber, Batch, Region) %>%
  distinct()

print(head(sample_metadata_df))
