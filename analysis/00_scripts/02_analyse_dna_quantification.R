# =========================================================
# 02_analyse_dna_quantification.R
# Analysis, summaries and plots for DNA quantification data
# =========================================================

# -----------------------------
# 1. Libraries
# -----------------------------
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(here)
library(purrr)
library(stringr)
library(tibble)
library(scales)

# -----------------------------
# 2. Paths
# -----------------------------
input_rds   <- here("analysis", "02_dna_quantification", "objects", "dna_quantification_clean.rds")

out_tables  <- here("analysis", "02_dna_quantification", "tables")
out_plots   <- here("analysis", "02_dna_quantification", "plots")
out_objects <- here("analysis", "02_dna_quantification", "objects")

dir.create(out_tables,  showWarnings = FALSE, recursive = TRUE)
dir.create(out_plots,   showWarnings = FALSE, recursive = TRUE)
dir.create(out_objects, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 3. Helper functions
# -----------------------------
safe_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
safe_sd   <- function(x) if (sum(!is.na(x)) < 2) NA_real_ else sd(x, na.rm = TRUE)
safe_med  <- function(x) if (all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)
safe_min  <- function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
safe_max  <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)

cv_percent <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 2) return(NA_real_)
  m <- mean(x)
  s <- sd(x)
  if (isTRUE(all.equal(m, 0))) return(NA_real_)
  100 * s / m
}

format_metric_label <- function(x) {
  recode(
    x,
    nanodrop_mean      = "Nanodrop DNA concentration",
    ratio_260_280_mean = "A260/A280 ratio",
    ratio_260_230_mean = "A260/A230 ratio",
    qubit_mean         = "Qubit DNA concentration",
    .default = x
  )
}

# Friedman test within each Species x Tissue using Individual as block
run_friedman <- function(df, value_col) {
  value_sym <- rlang::sym(value_col)
  
  tmp <- df %>%
    select(Species, Tissue, Individual, Kit, !!value_sym) %>%
    rename(value = !!value_sym) %>%
    group_by(Species, Tissue, Individual, Kit) %>%
    summarise(value = safe_mean(value), .groups = "drop")
  
  groups <- split(tmp, list(tmp$Species, tmp$Tissue), drop = TRUE)
  
  res <- lapply(groups, function(g) {
    wide <- g %>%
      select(Individual, Kit, value) %>%
      pivot_wider(names_from = Kit, values_from = value)
    
    wide_complete <- wide %>%
      filter(if_all(-Individual, ~ !is.na(.x)))
    
    n_blocks <- nrow(wide_complete)
    n_kits   <- ncol(wide_complete) - 1
    
    if (n_blocks < 2 || n_kits < 2) {
      return(tibble(
        metric = value_col,
        species = unique(g$Species)[1],
        tissue = unique(g$Tissue)[1],
        n_blocks = n_blocks,
        statistic = NA_real_,
        p_value = NA_real_
      ))
    }
    
    long_complete <- wide_complete %>%
      pivot_longer(-Individual, names_to = "Kit", values_to = "value")
    
    ft <- friedman.test(value ~ Kit | Individual, data = long_complete)
    
    tibble(
      metric = value_col,
      species = unique(g$Species)[1],
      tissue = unique(g$Tissue)[1],
      n_blocks = n_blocks,
      statistic = unname(ft$statistic),
      p_value = ft$p.value
    )
  })
  
  bind_rows(res)
}

# Pairwise paired Wilcoxon tests within each Species x Tissue
run_pairwise_wilcox <- function(df, value_col) {
  value_sym <- rlang::sym(value_col)
  kit_levels <- levels(df$Kit)
  pairs <- combn(kit_levels, 2, simplify = FALSE)
  
  tmp <- df %>%
    select(Species, Tissue, Individual, Kit, !!value_sym) %>%
    rename(value = !!value_sym) %>%
    group_by(Species, Tissue, Individual, Kit) %>%
    summarise(value = safe_mean(value), .groups = "drop")
  
  groups <- tmp %>%
    group_by(Species, Tissue) %>%
    group_split()
  
  out <- lapply(groups, function(g) {
    lapply(pairs, function(p) {
      sub <- g %>%
        filter(Kit %in% p) %>%
        select(Individual, Kit, value) %>%
        pivot_wider(names_from = Kit, values_from = value)
      
      if (!all(p %in% names(sub))) {
        return(tibble(
          metric = value_col,
          species = unique(g$Species)[1],
          tissue = unique(g$Tissue)[1],
          group1 = p[1],
          group2 = p[2],
          n_pairs = 0,
          statistic = NA_real_,
          p_value = NA_real_
        ))
      }
      
      sub_complete <- sub %>%
        filter(!is.na(.data[[p[1]]]), !is.na(.data[[p[2]]]))
      
      n_pairs <- nrow(sub_complete)
      
      if (n_pairs < 2) {
        return(tibble(
          metric = value_col,
          species = unique(g$Species)[1],
          tissue = unique(g$Tissue)[1],
          group1 = p[1],
          group2 = p[2],
          n_pairs = n_pairs,
          statistic = NA_real_,
          p_value = NA_real_
        ))
      }
      
      wt <- wilcox.test(
        x = sub_complete[[p[1]]],
        y = sub_complete[[p[2]]],
        paired = TRUE,
        exact = FALSE
      )
      
      tibble(
        metric = value_col,
        species = unique(g$Species)[1],
        tissue = unique(g$Tissue)[1],
        group1 = p[1],
        group2 = p[2],
        n_pairs = n_pairs,
        statistic = unname(wt$statistic),
        p_value = wt$p.value
      )
    }) %>% bind_rows()
  })
  
  bind_rows(out) %>%
    group_by(metric, species, tissue) %>%
    mutate(p_adj_bh = p.adjust(p_value, method = "BH")) %>%
    ungroup()
}

# -----------------------------
# 4. Load clean data
# -----------------------------
dna_clean <- readRDS(input_rds)

# -----------------------------
# 5. Create analysis-ready dataset
# -----------------------------
dna_analysis <- dna_clean %>%
  mutate(
    Nanodrop_analysis = if_else(`[ng/µL] Nanodrop` > 0, `[ng/µL] Nanodrop`, NA_real_),
    
    nanodrop_nonpositive_flag = `[ng/µL] Nanodrop` <= 0,
    ratio_260_280_invalid_flag = is.na(`Ratio 260/280`) | `Ratio 260/280` <= 0,
    ratio_260_230_invalid_flag = is.na(`Ratio 260/230`) | `Ratio 260/230` < 0,
    
    qubit_detected_flag = Qubit_status %in% c("Measured", "Out of range"),
    qubit_measured_flag = Qubit_status == "Measured"
  ) %>%
  arrange(Species, Tissue, Individual, Kit, Swabbing)

write_csv(
  dna_analysis,
  file.path(out_tables, "dna_quantification_analysis_ready.csv")
)

saveRDS(
  dna_analysis,
  file.path(out_objects, "dna_quantification_analysis_ready.rds")
)

# -----------------------------
# 6. Descriptive summaries
# -----------------------------

# 6.1 Summary by kit
summary_by_kit_raw <- dna_analysis %>%
  group_by(Kit) %>%
  summarise(
    n_samples = n(),
    
    nanodrop_mean   = safe_mean(Nanodrop_analysis),
    nanodrop_sd     = safe_sd(Nanodrop_analysis),
    nanodrop_median = safe_med(Nanodrop_analysis),
    nanodrop_min    = safe_min(Nanodrop_analysis),
    nanodrop_max    = safe_max(Nanodrop_analysis),
    
    ratio_260_280_mean   = safe_mean(`Ratio 260/280`),
    ratio_260_280_sd     = safe_sd(`Ratio 260/280`),
    ratio_260_280_median = safe_med(`Ratio 260/280`),
    
    ratio_260_230_mean   = safe_mean(`Ratio 260/230`),
    ratio_260_230_sd     = safe_sd(`Ratio 260/230`),
    ratio_260_230_median = safe_med(`Ratio 260/230`),
    
    qubit_n_measured     = sum(Qubit_status == "Measured"),
    qubit_n_out_of_range = sum(Qubit_status == "Out of range"),
    qubit_n_missing      = sum(Qubit_status == "Missing"),
    
    qubit_mean   = safe_mean(Qubit_ng_uL),
    qubit_sd     = safe_sd(Qubit_ng_uL),
    qubit_median = safe_med(Qubit_ng_uL),
    qubit_min    = safe_min(Qubit_ng_uL),
    qubit_max    = safe_max(Qubit_ng_uL),
    .groups = "drop"
  )

summary_by_kit <- summary_by_kit_raw %>%
  rename(
    `Number of samples` = n_samples,
    `Mean DNA concentration (Nanodrop, ng/µL)` = nanodrop_mean,
    `SD DNA concentration (Nanodrop, ng/µL)` = nanodrop_sd,
    `Median DNA concentration (Nanodrop, ng/µL)` = nanodrop_median,
    `Minimum DNA concentration (Nanodrop, ng/µL)` = nanodrop_min,
    `Maximum DNA concentration (Nanodrop, ng/µL)` = nanodrop_max,
    `Mean A260/A280 ratio` = ratio_260_280_mean,
    `SD A260/A280 ratio` = ratio_260_280_sd,
    `Median A260/A280 ratio` = ratio_260_280_median,
    `Mean A260/A230 ratio` = ratio_260_230_mean,
    `SD A260/A230 ratio` = ratio_260_230_sd,
    `Median A260/A230 ratio` = ratio_260_230_median,
    `Qubit: number measured` = qubit_n_measured,
    `Qubit: number out of range` = qubit_n_out_of_range,
    `Qubit: number missing` = qubit_n_missing,
    `Mean DNA concentration (Qubit, ng/µL)` = qubit_mean,
    `SD DNA concentration (Qubit, ng/µL)` = qubit_sd,
    `Median DNA concentration (Qubit, ng/µL)` = qubit_median,
    `Minimum DNA concentration (Qubit, ng/µL)` = qubit_min,
    `Maximum DNA concentration (Qubit, ng/µL)` = qubit_max
  )

write_csv(summary_by_kit, file.path(out_tables, "summary_by_kit.csv"))

# 6.2 Summary by kit x species x tissue
summary_by_kit_species_tissue_raw <- dna_analysis %>%
  group_by(Kit, Species, Tissue) %>%
  summarise(
    n_samples = n(),
    
    nanodrop_mean   = safe_mean(Nanodrop_analysis),
    nanodrop_sd     = safe_sd(Nanodrop_analysis),
    nanodrop_median = safe_med(Nanodrop_analysis),
    nanodrop_min    = safe_min(Nanodrop_analysis),
    nanodrop_max    = safe_max(Nanodrop_analysis),
    
    ratio_260_280_mean   = safe_mean(`Ratio 260/280`),
    ratio_260_280_sd     = safe_sd(`Ratio 260/280`),
    ratio_260_280_median = safe_med(`Ratio 260/280`),
    
    ratio_260_230_mean   = safe_mean(`Ratio 260/230`),
    ratio_260_230_sd     = safe_sd(`Ratio 260/230`),
    ratio_260_230_median = safe_med(`Ratio 260/230`),
    
    qubit_n_measured     = sum(Qubit_status == "Measured"),
    qubit_n_out_of_range = sum(Qubit_status == "Out of range"),
    qubit_n_missing      = sum(Qubit_status == "Missing"),
    
    qubit_mean   = safe_mean(Qubit_ng_uL),
    qubit_sd     = safe_sd(Qubit_ng_uL),
    qubit_median = safe_med(Qubit_ng_uL),
    .groups = "drop"
  ) %>%
  arrange(Species, Tissue, Kit)

summary_by_kit_species_tissue <- summary_by_kit_species_tissue_raw %>%
  rename(
    `Number of samples` = n_samples,
    `Mean DNA concentration (Nanodrop, ng/µL)` = nanodrop_mean,
    `SD DNA concentration (Nanodrop, ng/µL)` = nanodrop_sd,
    `Median DNA concentration (Nanodrop, ng/µL)` = nanodrop_median,
    `Minimum DNA concentration (Nanodrop, ng/µL)` = nanodrop_min,
    `Maximum DNA concentration (Nanodrop, ng/µL)` = nanodrop_max,
    `Mean A260/A280 ratio` = ratio_260_280_mean,
    `SD A260/A280 ratio` = ratio_260_280_sd,
    `Median A260/A280 ratio` = ratio_260_280_median,
    `Mean A260/A230 ratio` = ratio_260_230_mean,
    `SD A260/A230 ratio` = ratio_260_230_sd,
    `Median A260/A230 ratio` = ratio_260_230_median,
    `Qubit: number measured` = qubit_n_measured,
    `Qubit: number out of range` = qubit_n_out_of_range,
    `Qubit: number missing` = qubit_n_missing,
    `Mean DNA concentration (Qubit, ng/µL)` = qubit_mean,
    `SD DNA concentration (Qubit, ng/µL)` = qubit_sd,
    `Median DNA concentration (Qubit, ng/µL)` = qubit_median
  )

write_csv(summary_by_kit_species_tissue, file.path(out_tables, "summary_by_kit_species_tissue.csv"))

# 6.3 Summary of purity pass rates
purity_pass_by_kit_raw <- dna_analysis %>%
  group_by(Kit) %>%
  summarise(
    n_samples = n(),
    prop_260_280_17_20 = mean(pass_260_280_17_20, na.rm = TRUE),
    prop_260_280_15_20 = mean(pass_260_280_15_20, na.rm = TRUE),
    prop_260_230_15    = mean(pass_260_230_15, na.rm = TRUE),
    prop_260_230_175   = mean(pass_260_230_175, na.rm = TRUE),
    .groups = "drop"
  )

purity_pass_by_kit <- purity_pass_by_kit_raw %>%
  rename(
    `Number of samples` = n_samples,
    `Proportion within A260/A280 = 1.7-2.0` = prop_260_280_17_20,
    `Proportion within A260/A280 = 1.5-2.0` = prop_260_280_15_20,
    `Proportion with A260/A230 >= 1.5` = prop_260_230_15,
    `Proportion with A260/A230 >= 1.75` = prop_260_230_175
  )

write_csv(purity_pass_by_kit, file.path(out_tables, "purity_pass_by_kit.csv"))

purity_pass_by_kit_species_tissue_raw <- dna_analysis %>%
  group_by(Kit, Species, Tissue) %>%
  summarise(
    n_samples = n(),
    prop_260_280_17_20 = mean(pass_260_280_17_20, na.rm = TRUE),
    prop_260_280_15_20 = mean(pass_260_280_15_20, na.rm = TRUE),
    prop_260_230_15    = mean(pass_260_230_15, na.rm = TRUE),
    prop_260_230_175   = mean(pass_260_230_175, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Species, Tissue, Kit)

purity_pass_by_kit_species_tissue <- purity_pass_by_kit_species_tissue_raw %>%
  rename(
    `Number of samples` = n_samples,
    `Proportion within A260/A280 = 1.7-2.0` = prop_260_280_17_20,
    `Proportion within A260/A280 = 1.5-2.0` = prop_260_280_15_20,
    `Proportion with A260/A230 >= 1.5` = prop_260_230_15,
    `Proportion with A260/A230 >= 1.75` = prop_260_230_175
  )

write_csv(
  purity_pass_by_kit_species_tissue,
  file.path(out_tables, "purity_pass_by_kit_species_tissue.csv")
)

# 6.4 Qubit status by kit x species x tissue
qubit_status_by_kit_species_tissue <- dna_analysis %>%
  group_by(Kit, Species, Tissue, Qubit_status) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Kit, Species, Tissue) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  arrange(Species, Tissue, Kit, Qubit_status) %>%
  rename(
    `Number of samples` = n,
    `Proportion of samples` = prop
  )

write_csv(
  qubit_status_by_kit_species_tissue,
  file.path(out_tables, "qubit_status_by_kit_species_tissue.csv")
)

# 6.5 Sequencing selection summary
sequencing_by_kit_species_tissue <- dna_analysis %>%
  group_by(Kit, Species, Tissue, Sequencing) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(Species, Tissue, Kit, Sequencing) %>%
  rename(`Number of samples` = n)

write_csv(
  sequencing_by_kit_species_tissue,
  file.path(out_tables, "sequencing_by_kit_species_tissue.csv")
)

# -----------------------------
# 7. Individual-level paired dataset
# -----------------------------
dna_individual_raw <- dna_analysis %>%
  group_by(Species, Tissue, Individual, Kit) %>%
  summarise(
    n_rows = n(),
    nanodrop_mean = safe_mean(Nanodrop_analysis),
    ratio_260_280_mean = safe_mean(`Ratio 260/280`),
    ratio_260_230_mean = safe_mean(`Ratio 260/230`),
    qubit_mean = safe_mean(Qubit_ng_uL),
    qubit_status_mode = names(sort(table(Qubit_status), decreasing = TRUE))[1],
    .groups = "drop"
  ) %>%
  arrange(Species, Tissue, Individual, Kit)

dna_individual <- dna_individual_raw %>%
  rename(
    `Number of technical rows` = n_rows,
    `Mean DNA concentration (Nanodrop, ng/µL)` = nanodrop_mean,
    `Mean A260/A280 ratio` = ratio_260_280_mean,
    `Mean A260/A230 ratio` = ratio_260_230_mean,
    `Mean DNA concentration (Qubit, ng/µL)` = qubit_mean,
    `Most frequent Qubit status` = qubit_status_mode
  )

write_csv(
  dna_individual,
  file.path(out_tables, "dna_quantification_individual.csv")
)

saveRDS(
  dna_individual_raw,
  file.path(out_objects, "dna_quantification_individual.rds")
)

# 7.1 Individual-level summary by kit x matrix
summary_biological_replicates_raw <- dna_individual_raw %>%
  group_by(Kit, Species, Tissue) %>%
  summarise(
    n_individuals = n(),
    nanodrop_mean   = safe_mean(nanodrop_mean),
    nanodrop_sd     = safe_sd(nanodrop_mean),
    nanodrop_median = safe_med(nanodrop_mean),
    
    ratio_260_280_mean = safe_mean(ratio_260_280_mean),
    ratio_260_280_sd   = safe_sd(ratio_260_280_mean),
    
    ratio_260_230_mean = safe_mean(ratio_260_230_mean),
    ratio_260_230_sd   = safe_sd(ratio_260_230_mean),
    
    qubit_mean = safe_mean(qubit_mean),
    qubit_sd   = safe_sd(qubit_mean),
    .groups = "drop"
  ) %>%
  arrange(Species, Tissue, Kit)

summary_biological_replicates <- summary_biological_replicates_raw %>%
  rename(
    `Number of individuals` = n_individuals,
    `Mean DNA concentration (Nanodrop, ng/µL)` = nanodrop_mean,
    `SD DNA concentration (Nanodrop, ng/µL)` = nanodrop_sd,
    `Median DNA concentration (Nanodrop, ng/µL)` = nanodrop_median,
    `Mean A260/A280 ratio` = ratio_260_280_mean,
    `SD A260/A280 ratio` = ratio_260_280_sd,
    `Mean A260/A230 ratio` = ratio_260_230_mean,
    `SD A260/A230 ratio` = ratio_260_230_sd,
    `Mean DNA concentration (Qubit, ng/µL)` = qubit_mean,
    `SD DNA concentration (Qubit, ng/µL)` = qubit_sd
  )

write_csv(
  summary_biological_replicates,
  file.path(out_tables, "summary_biological_replicates.csv")
)

# -----------------------------
# 8. Technical reproducibility
# -----------------------------
technical_cv_raw <- dna_analysis %>%
  group_by(Species, Tissue, Individual, Kit) %>%
  summarise(
    n_rows = n(),
    nanodrop_cv = cv_percent(Nanodrop_analysis),
    ratio_260_280_cv = cv_percent(`Ratio 260/280`),
    ratio_260_230_cv = cv_percent(`Ratio 260/230`),
    qubit_cv = cv_percent(Qubit_ng_uL),
    .groups = "drop"
  )

technical_cv <- technical_cv_raw %>%
  rename(
    `Number of technical rows` = n_rows,
    `CV Nanodrop DNA concentration (%)` = nanodrop_cv,
    `CV A260/A280 ratio (%)` = ratio_260_280_cv,
    `CV A260/A230 ratio (%)` = ratio_260_230_cv,
    `CV Qubit DNA concentration (%)` = qubit_cv
  )

write_csv(
  technical_cv,
  file.path(out_tables, "technical_cv_by_individual_kit.csv")
)

technical_cv_summary_raw <- technical_cv_raw %>%
  group_by(Kit, Species, Tissue) %>%
  summarise(
    n_groups = n(),
    nanodrop_cv_mean = safe_mean(nanodrop_cv),
    ratio_260_280_cv_mean = safe_mean(ratio_260_280_cv),
    ratio_260_230_cv_mean = safe_mean(ratio_260_230_cv),
    qubit_cv_mean = safe_mean(qubit_cv),
    .groups = "drop"
  ) %>%
  arrange(Species, Tissue, Kit)

technical_cv_summary <- technical_cv_summary_raw %>%
  rename(
    `Number of individual-kit groups` = n_groups,
    `Mean CV Nanodrop DNA concentration (%)` = nanodrop_cv_mean,
    `Mean CV A260/A280 ratio (%)` = ratio_260_280_cv_mean,
    `Mean CV A260/A230 ratio (%)` = ratio_260_230_cv_mean,
    `Mean CV Qubit DNA concentration (%)` = qubit_cv_mean
  )

write_csv(
  technical_cv_summary,
  file.path(out_tables, "technical_cv_summary.csv")
)

# -----------------------------
# 9. Inferential statistics
# -----------------------------
metrics_to_test <- c(
  "nanodrop_mean",
  "ratio_260_280_mean",
  "ratio_260_230_mean",
  "qubit_mean"
)

friedman_results_raw <- bind_rows(
  lapply(metrics_to_test, function(m) run_friedman(dna_individual_raw, m))
) %>%
  mutate(metric = format_metric_label(metric))

friedman_results <- friedman_results_raw %>%
  rename(
    Metric = metric,
    Species = species,
    Tissue = tissue,
    `Number of paired individuals` = n_blocks,
    `Friedman statistic` = statistic,
    `P-value` = p_value
  )

write_csv(
  friedman_results,
  file.path(out_tables, "friedman_results_by_matrix.csv")
)

pairwise_paired_wilcox_raw <- bind_rows(
  lapply(metrics_to_test, function(m) run_pairwise_wilcox(dna_individual_raw, m))
) %>%
  mutate(metric = format_metric_label(metric))

pairwise_paired_wilcox <- pairwise_paired_wilcox_raw %>%
  rename(
    Metric = metric,
    Species = species,
    Tissue = tissue,
    `Comparison 1` = group1,
    `Comparison 2` = group2,
    `Number of paired individuals` = n_pairs,
    `Wilcoxon statistic` = statistic,
    `P-value` = p_value,
    `BH-adjusted P-value` = p_adj_bh
  )

write_csv(
  pairwise_paired_wilcox,
  file.path(out_tables, "pairwise_paired_wilcox_by_matrix.csv")
)

# -----------------------------
# 10. Plots
# -----------------------------
theme_set(theme_bw())

# 10.1 Nanodrop by kit
p_nanodrop <- ggplot(
  dna_analysis,
  aes(x = Kit, y = Nanodrop_analysis, fill = Kit)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.7) +
  facet_grid(Species ~ Tissue, scales = "free_y") +
  labs(
    title = "Nanodrop DNA concentration by extraction kit",
    x = "Extraction kit",
    y = "DNA concentration (ng/µL)"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  filename = file.path(out_plots, "nanodrop_by_kit_species_tissue.png"),
  plot = p_nanodrop,
  width = 10, height = 6, dpi = 300
)

# 10.2 A260/A280
p_260_280 <- ggplot(
  dna_analysis,
  aes(x = Kit, y = `Ratio 260/280`, fill = Kit)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.7) +
  geom_hline(yintercept = c(1.7, 2.0), linetype = "dashed", color = "red") +
  facet_grid(Species ~ Tissue) +
  labs(
    title = "A260/A280 ratio by extraction kit",
    x = "Extraction kit",
    y = "A260/A280"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  filename = file.path(out_plots, "ratio_260_280_by_kit_species_tissue.png"),
  plot = p_260_280,
  width = 10, height = 6, dpi = 300
)

# 10.3 A260/A230
p_260_230 <- ggplot(
  dna_analysis,
  aes(x = Kit, y = `Ratio 260/230`, fill = Kit)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.7) +
  geom_hline(yintercept = c(1.5, 1.75), linetype = "dashed", color = "red") +
  facet_grid(Species ~ Tissue) +
  labs(
    title = "A260/A230 ratio by extraction kit",
    x = "Extraction kit",
    y = "A260/A230"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  filename = file.path(out_plots, "ratio_260_230_by_kit_species_tissue.png"),
  plot = p_260_230,
  width = 10, height = 6, dpi = 300
)

# 10.4 Qubit status
p_qubit_status <- ggplot(
  qubit_status_by_kit_species_tissue,
  aes(x = Kit, y = `Proportion of samples`, fill = Qubit_status)
) +
  geom_col(position = "fill") +
  facet_grid(Species ~ Tissue) +
  scale_y_continuous(labels = percent_format()) +
  labs(
    title = "Qubit measurement status by extraction kit",
    x = "Extraction kit",
    y = "Percentage of samples",
    fill = "Qubit status"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  filename = file.path(out_plots, "qubit_status_by_kit_species_tissue.png"),
  plot = p_qubit_status,
  width = 10, height = 6, dpi = 300
)

# 10.5 Qubit measured only
dna_qubit_measured <- dna_analysis %>%
  filter(Qubit_status == "Measured")

p_qubit_measured <- ggplot(
  dna_qubit_measured,
  aes(x = Kit, y = Qubit_ng_uL, fill = Kit)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.7) +
  facet_grid(Species ~ Tissue, scales = "free_y") +
  labs(
    title = "Qubit DNA concentration by extraction kit (measured samples only)",
    x = "Extraction kit",
    y = "Qubit DNA concentration (ng/µL)"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  filename = file.path(out_plots, "qubit_measured_by_kit_species_tissue.png"),
  plot = p_qubit_measured,
  width = 10, height = 6, dpi = 300
)

# 10.6 Purity pass rates
purity_long <- purity_pass_by_kit_species_tissue %>%
  pivot_longer(
    cols = c(
      `Proportion within A260/A280 = 1.7-2.0`,
      `Proportion within A260/A280 = 1.5-2.0`,
      `Proportion with A260/A230 >= 1.5`,
      `Proportion with A260/A230 >= 1.75`
    ),
    names_to = "metric",
    values_to = "proportion"
  ) %>%
  mutate(
    metric = recode(
      metric,
      `Proportion within A260/A280 = 1.7-2.0` = "A260/A280 in 1.7-2.0",
      `Proportion within A260/A280 = 1.5-2.0` = "A260/A280 in 1.5-2.0",
      `Proportion with A260/A230 >= 1.5` = "A260/A230 ≥ 1.5",
      `Proportion with A260/A230 >= 1.75` = "A260/A230 ≥ 1.75"
    )
  )

p_purity_pass <- ggplot(
  purity_long,
  aes(x = Kit, y = proportion, fill = Kit)
) +
  geom_col() +
  facet_grid(metric ~ Species + Tissue) +
  scale_y_continuous(labels = percent_format()) +
  labs(
    title = "Proportion of samples passing common purity thresholds",
    x = "Extraction kit",
    y = "Proportion"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  filename = file.path(out_plots, "purity_pass_rates_by_kit_species_tissue.png"),
  plot = p_purity_pass,
  width = 12, height = 8, dpi = 300
)

# 10.7 Paired spaghetti plot: Nanodrop at biological replicate level
p_spaghetti_nanodrop <- ggplot(
  dna_individual_raw,
  aes(x = Kit, y = nanodrop_mean, group = Individual, color = factor(Individual))
) +
  geom_line(alpha = 0.6) +
  geom_point(size = 2) +
  facet_grid(Species ~ Tissue, scales = "free_y") +
  labs(
    title = "Individual-level paired Nanodrop concentration across extraction kits",
    x = "Extraction kit",
    y = "Mean Nanodrop concentration (ng/µL)",
    color = "Individual"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  filename = file.path(out_plots, "paired_nanodrop_spaghetti.png"),
  plot = p_spaghetti_nanodrop,
  width = 10, height = 6, dpi = 300
)

# -----------------------------
# 11. Save master objects list
# -----------------------------
analysis_objects <- list(
  dna_analysis = dna_analysis,
  summary_by_kit = summary_by_kit,
  summary_by_kit_species_tissue = summary_by_kit_species_tissue,
  purity_pass_by_kit = purity_pass_by_kit,
  purity_pass_by_kit_species_tissue = purity_pass_by_kit_species_tissue,
  qubit_status_by_kit_species_tissue = qubit_status_by_kit_species_tissue,
  sequencing_by_kit_species_tissue = sequencing_by_kit_species_tissue,
  dna_individual = dna_individual_raw,
  summary_biological_replicates = summary_biological_replicates,
  technical_cv = technical_cv,
  technical_cv_summary = technical_cv_summary,
  friedman_results = friedman_results,
  pairwise_paired_wilcox = pairwise_paired_wilcox
)

saveRDS(
  analysis_objects,
  file.path(out_objects, "dna_quantification_analysis_objects.rds")
)

# -----------------------------
# 12. Console summary for interpretation
# -----------------------------
summary_for_console <- summary_by_kit %>%
  left_join(
    purity_pass_by_kit %>%
      select(
        Kit,
        `Proportion within A260/A280 = 1.7-2.0`,
        `Proportion with A260/A230 >= 1.5`
      ),
    by = "Kit"
  ) %>%
  mutate(
    `Proportion Qubit measured` = `Qubit: number measured` / `Number of samples`
  ) %>%
  select(
    Kit,
    `Number of samples`,
    `Mean DNA concentration (Nanodrop, ng/µL)`,
    `Mean A260/A280 ratio`,
    `Mean A260/A230 ratio`,
    `Mean DNA concentration (Qubit, ng/µL)`,
    `Proportion Qubit measured`,
    `Proportion within A260/A280 = 1.7-2.0`,
    `Proportion with A260/A230 >= 1.5`
  ) %>%
  arrange(desc(`Mean DNA concentration (Nanodrop, ng/µL)`))

friedman_sig <- friedman_results %>%
  filter(!is.na(`P-value`), `P-value` < 0.05) %>%
  arrange(Metric, Species, Tissue)

wilcox_sig <- pairwise_paired_wilcox %>%
  filter(!is.na(`BH-adjusted P-value`), `BH-adjusted P-value` < 0.05) %>%
  arrange(Metric, Species, Tissue, `BH-adjusted P-value`)

cat("\n====================================================\n")
cat("DNA QUANTIFICATION: MAIN SUMMARY BY KIT\n")
cat("====================================================\n\n")
print(summary_for_console, n = Inf, width = Inf)

cat("\n====================================================\n")
cat("SIGNIFICANT FRIEDMAN TESTS (P < 0.05)\n")
cat("====================================================\n\n")
if (nrow(friedman_sig) == 0) {
  cat("No significant Friedman tests detected.\n")
} else {
  print(friedman_sig, n = Inf, width = Inf)
}

cat("\n====================================================\n")
cat("SIGNIFICANT PAIRWISE WILCOXON TESTS (BH-adjusted P < 0.05)\n")
cat("====================================================\n\n")
if (nrow(wilcox_sig) == 0) {
  cat("No significant pairwise Wilcoxon comparisons detected.\n")
} else {
  print(wilcox_sig, n = Inf, width = Inf)
}

cat("\n====================================================\n")
cat("TOP-LEVEL INTERPRETATION AID\n")
cat("====================================================\n\n")

top_nanodrop <- summary_for_console %>%
  slice_max(order_by = `Mean DNA concentration (Nanodrop, ng/µL)`, n = 1, with_ties = FALSE)

top_qubit <- summary_for_console %>%
  filter(!is.na(`Mean DNA concentration (Qubit, ng/µL)`)) %>%
  slice_max(order_by = `Mean DNA concentration (Qubit, ng/µL)`, n = 1, with_ties = FALSE)

top_260280 <- summary_for_console %>%
  filter(!is.na(`Proportion within A260/A280 = 1.7-2.0`)) %>%
  slice_max(order_by = `Proportion within A260/A280 = 1.7-2.0`, n = 1, with_ties = FALSE)

top_260230 <- summary_for_console %>%
  filter(!is.na(`Proportion with A260/A230 >= 1.5`)) %>%
  slice_max(order_by = `Proportion with A260/A230 >= 1.5`, n = 1, with_ties = FALSE)

cat("Highest mean Nanodrop concentration:\n")
print(top_nanodrop, n = Inf, width = Inf)

cat("\nHighest mean Qubit concentration:\n")
print(top_qubit, n = Inf, width = Inf)

cat("\nBest proportion within A260/A280 = 1.7-2.0:\n")
print(top_260280, n = Inf, width = Inf)

cat("\nBest proportion with A260/A230 >= 1.5:\n")
print(top_260230, n = Inf, width = Inf)

message("DNA quantification analysis complete.")