# =========================================================
# 01_prepare_dna_quantification_data.R
# Preparation and cleaning of DNA quantification data
# =========================================================

# -----------------------------
# 1. Libraries
# -----------------------------
library(readxl)
library(dplyr)
library(stringr)
library(readr)
library(here)

# -----------------------------
# 2. Paths
# -----------------------------
input_file  <- here("analysis", "01_data", "dna_quantification.xlsx")

out_tables  <- here("analysis", "02_dna_quantification", "tables")
out_objects <- here("analysis", "02_dna_quantification", "objects")

dir.create(out_tables,  showWarnings = FALSE, recursive = TRUE)
dir.create(out_objects, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 3. Load data
# -----------------------------
# Assumption:
# The first sheet contains the curated table with columns such as:
# IIM ID, Novogene ID, Species, Individual, Swabbing, Tissue,
# Extraction date, Kit, [ng/µL] Nanodrop, Ratio 260/280,
# Ratio 260/230, [ng/µL] Qubit, Sequencing

dna_raw <- read_excel(input_file, sheet = 1)

# -----------------------------
# 4. Basic checks
# -----------------------------
required_cols <- c(
  "IIM ID",
  "Novogene ID",
  "Species",
  "Individual",
  "Swabbing",
  "Tissue",
  "Extraction date",
  "Kit",
  "[ng/µL] Nanodrop",
  "Ratio 260/280",
  "Ratio 260/230",
  "[ng/µL] Qubit",
  "Sequencing"
)

missing_cols <- setdiff(required_cols, names(dna_raw))

if (length(missing_cols) > 0) {
  stop(
    "The following required columns are missing from the input file: ",
    paste(missing_cols, collapse = ", ")
  )
}

# -----------------------------
# 5. Clean and prepare data
# -----------------------------
dna_clean <- dna_raw %>%
  mutate(
    Species = str_trim(as.character(Species)),
    Tissue  = str_trim(as.character(Tissue)),
    Kit     = str_trim(as.character(Kit)),
    Sequencing = str_trim(as.character(Sequencing)),
    
    # Factors with controlled order
    Species = factor(Species, levels = c("Octopus", "Turbot")),
    Tissue  = factor(Tissue, levels = c("Skin", "Mantle", "Gut")),
    Kit     = factor(Kit, levels = c("MasterPure", "PowerFecal", "PowerSoil", "PureLink")),
    
    # Numeric / date conversion
    Individual = as.integer(Individual),
    Swabbing   = as.integer(Swabbing),
    `Extraction date` = as.Date(`Extraction date`),
    
    `[ng/µL] Nanodrop` = as.numeric(`[ng/µL] Nanodrop`),
    `Ratio 260/280`    = as.numeric(`Ratio 260/280`),
    `Ratio 260/230`    = as.numeric(`Ratio 260/230`),
    
    # Keep original Qubit as character
    Qubit_raw = as.character(`[ng/µL] Qubit`),
    
    # Qubit status
    Qubit_status = case_when(
      is.na(Qubit_raw) ~ "Missing",
      str_to_lower(Qubit_raw) %in% c("out of range", "oor") ~ "Out of range",
      TRUE ~ "Measured"
    ),
    
    # Numeric Qubit value only when measured
    Qubit_ng_uL = case_when(
      Qubit_status == "Measured" ~ suppressWarnings(as.numeric(Qubit_raw)),
      TRUE ~ NA_real_
    ),
    
    # Optional version coding "Out of range" as zero
    Qubit_ng_uL_zero = case_when(
      Qubit_status == "Measured" ~ suppressWarnings(as.numeric(Qubit_raw)),
      Qubit_status == "Out of range" ~ 0,
      TRUE ~ NA_real_
    ),
    
    # Common quality flags
    pass_260_280_17_20 = dplyr::between(`Ratio 260/280`, 1.7, 2.0),
    pass_260_280_15_20 = dplyr::between(`Ratio 260/280`, 1.5, 2.0),
    pass_260_230_15    = `Ratio 260/230` >= 1.5,
    pass_260_230_175   = `Ratio 260/230` >= 1.75,
    
    # Sequencing flag
    sequencing_selected = case_when(
      str_to_upper(Sequencing) == "Y" ~ TRUE,
      str_to_upper(Sequencing) == "N" ~ FALSE,
      TRUE ~ NA
    )
  ) %>%
  arrange(Species, Tissue, Individual, Kit, Swabbing)

# -----------------------------
# 6. Design/structure checks
# -----------------------------

# Full count by biological replicate and kit
design_check <- dna_clean %>%
  count(Species, Tissue, Individual, Kit, name = "n_rows") %>%
  arrange(Species, Tissue, Individual, Kit)

# Coverage summary by matrix and kit
coverage_summary <- dna_clean %>%
  distinct(Species, Tissue, Individual, Kit) %>%
  count(Species, Tissue, Kit, name = "n_individuals") %>%
  arrange(Species, Tissue, Kit)

# Count by Qubit status
qubit_status_overview <- dna_clean %>%
  count(Qubit_status, name = "n") %>%
  mutate(prop = n / sum(n))

# Quick missingness overview
missingness_summary <- tibble(
  variable = c(
    "Nanodrop",
    "Ratio 260/280",
    "Ratio 260/230",
    "Qubit numeric"
  ),
  n_missing = c(
    sum(is.na(dna_clean$`[ng/µL] Nanodrop`)),
    sum(is.na(dna_clean$`Ratio 260/280`)),
    sum(is.na(dna_clean$`Ratio 260/230`)),
    sum(is.na(dna_clean$Qubit_ng_uL))
  )
)

# -----------------------------
# 7. Export tables
# -----------------------------
write_csv(
  dna_clean,
  file.path(out_tables, "dna_quantification_clean.csv")
)

write_csv(
  design_check,
  file.path(out_tables, "dna_quantification_design_check.csv")
)

write_csv(
  coverage_summary,
  file.path(out_tables, "dna_quantification_coverage_summary.csv")
)

write_csv(
  qubit_status_overview,
  file.path(out_tables, "dna_quantification_qubit_status_overview.csv")
)

write_csv(
  missingness_summary,
  file.path(out_tables, "dna_quantification_missingness_summary.csv")
)

# -----------------------------
# 8. Export objects
# -----------------------------
saveRDS(
  dna_clean,
  file.path(out_objects, "dna_quantification_clean.rds")
)

saveRDS(
  list(
    dna_clean = dna_clean,
    design_check = design_check,
    coverage_summary = coverage_summary,
    qubit_status_overview = qubit_status_overview,
    missingness_summary = missingness_summary
  ),
  file.path(out_objects, "dna_quantification_preparation_objects.rds")
)

# -----------------------------
# 9. Console messages
# -----------------------------
message("DNA quantification data preparation complete.")
message("Clean dataset saved to: analysis/02_dna_quantification/objects/")
message("Auxiliary tables saved to: analysis/02_dna_quantification/tables/")