# 08_aim1_representativeness.R
# Purpose: Aim 1 -- Quantify representativeness of Stanford ADRC sample
#          relative to the ACS catchment population (Alameda, San Mateo,
#          Santa Clara counties).
#
# Generates:
#   - Standardized mean difference (SMD) comparison: ADRC vs. BRFSS
#   - Love plots (before/after weighting)
#   - Propensity score overlap/positivity diagnostics
#   - Table 1 (descriptive characteristics, ADRC vs. BRFSS)
#
# INPUTS:
#   - processed_data_dir/05_imputed_with_weights.rds
#   - ACS 5-year estimates (2019-2023) for 3 counties (to be loaded separately)
#
# OUTPUTS:
#   - results_dir/08_table1.rds
#   - results_dir/08_love_plot.pdf
#   - results_dir/08_pscore_overlap.pdf
#   - results_dir/08_positivity_diagnostics.rds
#
# ============================================================================

source("config.R")
library(twang)
library(survey)
# library(tidycensus)  # uncomment when ACS API key is configured


# ============================================================================
# 1. LOAD DATA
# ============================================================================

# For Aim 1 representativeness: use ALL-AGES harmonized data (no age filter)
# to show the full picture of how ADRC compares to the catchment community.
# The 5-category race_eth_h is used here (not the 4-category collapse).
harmonized_all <- readRDS(file.path(processed_data_dir,
                                     "03_adrc_brfss_harmonized_all_ages.rds"))
adrc_all  <- harmonized_all %>% filter(adrc_h == 1L)
brfss_all <- harmonized_all %>% filter(adrc_h == 0L)

# Also load weighted data for love plot (age 50+ only, post-imputation)
imputed_weighted <- readRDS(file.path(processed_data_dir, "05_imputed_with_weights.rds"))
imp1 <- imputed_weighted %>% filter(imp_h == 1)
adrc  <- imp1 %>% filter(adrc_h == 1L)
brfss <- imp1 %>% filter(adrc_h == 0L)


# ============================================================================
# 2. ACS DATA (3-COUNTY, 2019-2023 5-YEAR ESTIMATES)
# ============================================================================
# TODO: Fetch via tidycensus::get_acs() or load pre-downloaded file.
#
# Key ACS tables:
#   B01001 - Sex by age
#   B03002 - Hispanic/Latino origin by race
#   B15003 - Educational attainment (population 25+)
#   B12001 - Marital status (population 15+)
#   B19013 - Median household income
#
# County FIPS: 06001 (Alameda), 06081 (San Mateo), 06085 (Santa Clara)

message("NOTE: ACS comparison data not yet loaded.")


# ============================================================================
# 3. TABLE 1: DESCRIPTIVE CHARACTERISTICS
# ============================================================================

covariates <- c("age_h", "male_h", "race_eth_h", "education_h", "marital_h",
                "diabetes_h", "depression_h", "mi_h", "angina_h",
                "stroke_h", "hypertension_h")

# ADRC unweighted summary (aggregate counts only)
adrc_summary <- adrc %>%
  summarise(
    n = n(),
    age_mean = mean(age_h, na.rm = TRUE),
    age_sd   = sd(age_h, na.rm = TRUE),
    across(c(male_h, diabetes_h, depression_h, mi_h, angina_h,
             stroke_h, hypertension_h),
           list(mean = ~mean(., na.rm = TRUE)),
           .names = "{.col}_mean")
  )

# BRFSS survey-weighted summary
brfss_design <- svydesign(ids = ~1, weights = ~brfss_sampwt_h, data = brfss)

# TODO: Format Table 1 with small cell suppression (n < 11) and save


# ============================================================================
# 4. LOVE PLOT (SMDs BEFORE/AFTER WEIGHTING)
# ============================================================================

# Compute SMDs for all covariates including factor levels (dummy-coded)
compute_full_smd <- function(data, weight_var = NULL) {
  adrc_d  <- data %>% filter(adrc_h == 1)
  brfss_d <- data %>% filter(adrc_h == 0)
  brfss_w <- brfss_d$brfss_sampwt_h
  adrc_w  <- if (is.null(weight_var)) rep(1, nrow(adrc_d)) else adrc_d[[weight_var]]

  # Helper: compute SMD for a single numeric vector
  one_smd <- function(va, vb) {
    ma   <- weighted.mean(va, adrc_w, na.rm = TRUE)
    mb   <- weighted.mean(vb, brfss_w, na.rm = TRUE)
    sd_b <- sqrt(weighted.mean((vb - mb)^2, brfss_w, na.rm = TRUE))
    if (sd_b > 0) (ma - mb) / sd_b else NA_real_
  }

  results <- list()

  # Continuous / binary covariates
  for (v in c("age_h", "male_h", "diabetes_h", "depression_h",
              "mi_h", "angina_h", "stroke_h", "hypertension_h")) {
    results[[v]] <- one_smd(adrc_d[[v]], brfss_d[[v]])
  }

  # Factor covariates: compute SMD for each level
  for (v in c("race_eth4_h", "education_h", "marital_h")) {
    for (lev in levels(data[[v]])) {
      label <- paste0(v, ": ", lev)
      results[[label]] <- one_smd(
        as.numeric(adrc_d[[v]] == lev),
        as.numeric(brfss_d[[v]] == lev)
      )
    }
  }

  tibble(covariate = names(results), smd = unlist(results))
}

smd_unweighted <- compute_full_smd(imp1, weight_var = NULL) %>%
  rename(unweighted = smd)
smd_weighted <- compute_full_smd(imp1, weight_var = "sw3_final") %>%
  rename(weighted = smd)

smd_df <- smd_unweighted %>%
  left_join(smd_weighted, by = "covariate") %>%
  pivot_longer(cols = c(unweighted, weighted),
               names_to = "method", values_to = "smd") %>%
  # Clean up labels for display
  mutate(
    label = covariate %>%
      str_replace("_h$", "") %>%
      str_replace("race_eth4_h: ", "Race: ") %>%
      str_replace("education_h: ", "Edu: ") %>%
      str_replace("marital_h: ", "Marital: ") %>%
      str_replace("age$", "Age") %>%
      str_replace("male$", "Male") %>%
      str_replace("diabetes$", "Diabetes") %>%
      str_replace("depression$", "Depression") %>%
      str_replace("mi$", "Heart attack / MI") %>%
      str_replace("angina$", "Angina / CHD") %>%
      str_replace("stroke$", "Stroke") %>%
      str_replace("hypertension$", "Hypertension"),
    # Group for faceting / ordering
    group = case_when(
      str_detect(covariate, "race_eth4") ~ "Race/Ethnicity",
      str_detect(covariate, "education") ~ "Education",
      str_detect(covariate, "marital")   ~ "Marital Status",
      covariate %in% c("age_h", "male_h") ~ "Demographics",
      TRUE ~ "Health Conditions"
    ),
    group = factor(group, levels = c("Demographics", "Race/Ethnicity",
                                      "Education", "Marital Status",
                                      "Health Conditions"))
  )

# Order covariates by absolute unweighted SMD within each group
smd_order <- smd_df %>%
  filter(method == "unweighted") %>%
  arrange(group, abs(smd)) %>%
  pull(label)
smd_df$label <- factor(smd_df$label, levels = smd_order)

love_plot <- ggplot(smd_df, aes(x = smd, y = label,
                                 color = method, shape = method)) +
  geom_point(size = 2.5) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black") +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed",
             color = "red", alpha = 0.5) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  scale_color_manual(
    values = c("unweighted" = "grey50", "weighted" = "steelblue"),
    labels = c("unweighted" = "Unweighted", "weighted" = "IOSW Weighted")
  ) +
  scale_shape_manual(
    values = c("unweighted" = 1, "weighted" = 16),
    labels = c("unweighted" = "Unweighted", "weighted" = "IOSW Weighted")
  ) +
  labs(x = "Standardized Mean Difference (ADRC \u2212 BRFSS) / SD(BRFSS)",
       y = NULL,
       title = "Covariate Balance: Stanford ADRC vs. 3-County BRFSS Target Population",
       subtitle = "Dashed red lines indicate \u00b10.25 SMD threshold",
       color = NULL, shape = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    strip.text.y = element_text(angle = 0, hjust = 0, face = "bold"),
    panel.grid.major.y = element_line(color = "grey92"),
    plot.title = element_text(size = 12, face = "bold")
  )

ggsave(file.path(results_dir, "08_love_plot.pdf"), love_plot,
       width = 9, height = 10)


# ============================================================================
# 5. PROPENSITY SCORE OVERLAP / POSITIVITY DIAGNOSTICS
# ============================================================================

pscore_plot <- ggplot(imp1, aes(x = p3, fill = factor(adrc_h),
                                 color = factor(adrc_h))) +
  geom_density(alpha = 0.4) +
  scale_fill_manual(values = c("0" = "grey70", "1" = "steelblue"),
                    labels = c("BRFSS", "ADRC")) +
  scale_color_manual(values = c("0" = "grey50", "1" = "steelblue4"),
                     labels = c("BRFSS", "ADRC")) +
  labs(x = "Propensity Score: P(ADRC | X)",
       y = "Density",
       title = "Propensity Score Overlap",
       fill = NULL, color = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(file.path(results_dir, "08_pscore_overlap.pdf"), pscore_plot,
       width = 8, height = 5)

positivity <- adrc %>%
  summarise(
    n                = n(),
    min_pscore       = min(p3, na.rm = TRUE),
    max_pscore       = max(p3, na.rm = TRUE),
    pct_pscore_lt01  = mean(p3 < 0.01, na.rm = TRUE) * 100,
    pct_pscore_gt99  = mean(p3 > 0.99, na.rm = TRUE) * 100,
    max_weight       = max(sw3, na.rm = TRUE),
    pct_weight_gt10  = mean(sw3 > 10, na.rm = TRUE) * 100
  )

saveRDS(positivity, file.path(results_dir, "08_positivity_diagnostics.rds"))

message("Positivity diagnostics:")
message(paste(capture.output(positivity), collapse = "\n"))
message("\nSaved: love plot, pscore overlap, positivity diagnostics")
