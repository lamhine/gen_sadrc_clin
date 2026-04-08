# 09_results_summary.R
# Purpose: Aim 3 -- Generate publication-ready figures, tables, and summaries.
#
# INPUTS:
#   - results_dir/07_aipsw_results.rds
#   - results_dir/07_aipsw_bootests.rds
#   - results_dir/08_positivity_diagnostics.rds
#
# OUTPUTS:
#   - results_dir/09_forest_estimator_comparison.pdf
#   - results_dir/09_forest_race_stratified.pdf
#   - results_dir/09_results_table.rds
#
# ============================================================================

source("config.R")


# ============================================================================
# 1. LOAD RESULTS
# ============================================================================

results    <- readRDS(file.path(results_dir, "07_aipsw_results.rds"))
positivity <- readRDS(file.path(results_dir, "08_positivity_diagnostics.rds"))


# ============================================================================
# 2. FOREST PLOT: ESTIMATOR COMPARISON
# ============================================================================

overall <- results %>%
  filter(estimand %in% c("aipsw", "iosw", "gcomp", "unweighted")) %>%
  mutate(
    estimand = factor(estimand,
                      levels = c("unweighted", "iosw", "gcomp", "aipsw"),
                      labels = c("Unweighted ADRC", "IOSW only",
                                 "Outcome model only", "AIPSW (DR)")),
    outcome = factor(outcome,
                     levels = c("mci_h", "ad_h", "adrd_h", "moca_h", "cdrsb_h"),
                     labels = c("MCI", "AD", "ADRD", "MoCA", "CDR-SB"))
  )

p_estimator <- ggplot(overall, aes(x = estimate, y = estimand,
                                    xmin = lci, xmax = uci,
                                    color = estimand)) +
  geom_pointrange(size = 0.8) +
  facet_wrap(~outcome, scales = "free_x", ncol = 2) +
  scale_color_manual(values = c(
    "Unweighted ADRC"    = "grey50",
    "IOSW only"          = "forestgreen",
    "Outcome model only" = "darkorange",
    "AIPSW (DR)"         = "steelblue"
  )) +
  labs(x = "Estimate (95% CI)", y = NULL,
       title = "Generalized Estimates: Estimator Comparison",
       color = NULL) +
  theme_minimal() +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold"))

ggsave(file.path(results_dir, "09_forest_estimator_comparison.pdf"),
       p_estimator, width = 10, height = 6)


# ============================================================================
# 3. FOREST PLOT: RACE/ETHNICITY-STRATIFIED AIPSW
# ============================================================================

race_results <- results %>%
  filter(grepl("^aipsw_", estimand)) %>%
  mutate(
    race = str_replace(estimand, "^aipsw_", "") %>%
      str_replace_all("_", " ") %>%
      str_to_title() %>%
      factor(levels = c("White Nh", "Hispanic", "Asian Nh", "Nh Other")),
    outcome = factor(outcome,
                     levels = c("mci_h", "ad_h", "adrd_h", "moca_h", "cdrsb_h"),
                     labels = c("MCI", "AD", "ADRD", "MoCA", "CDR-SB"))
  ) %>%
  filter(!is.na(race))

p_race <- ggplot(race_results, aes(x = estimate, y = race,
                                    xmin = lci, xmax = uci,
                                    color = race)) +
  geom_pointrange(size = 0.8) +
  facet_wrap(~outcome, scales = "free_x", ncol = 2) +
  labs(x = "AIPSW Estimate (95% CI)", y = NULL,
       title = "Generalized Estimates by Race/Ethnicity (AIPSW)",
       color = NULL) +
  theme_minimal() +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold"))

ggsave(file.path(results_dir, "09_forest_race_stratified.pdf"),
       p_race, width = 10, height = 6)


# ============================================================================
# 4. PREVALENCE RATIOS AND DIFFERENCES (vs. White NH reference)
# ============================================================================
# Compute from bootstrap estimates stored in 07_aipsw_bootests.rds.
# For each bootstrap replicate:
#   PR = est_race / est_white_nh
#   PD = est_race - est_white_nh
# Then take mean + 2.5th/97.5th percentiles for CIs.

bootests <- readRDS(file.path(results_dir, "07_aipsw_bootests.rds"))

# Race column names in the bootstrap matrix
race_cols <- c("aipsw_white_nh", "aipsw_hispanic", "aipsw_asian_nh", "aipsw_nh_other")
race_labels <- c("White NH", "Hispanic", "Asian NH", "NH Other")
ref_col <- "aipsw_white_nh"

pr_pd_list <- list()

for (out_name in names(bootests)) {
  bmat <- bootests[[out_name]]

  # For each non-reference race, compute PR and PD vs White NH
  for (i in seq_along(race_cols)) {
    rc <- race_cols[i]
    if (rc == ref_col) next  # skip reference group

    boot_ref  <- bmat[, ref_col]
    boot_race <- bmat[, rc]

    # Prevalence ratio
    boot_pr <- boot_race / boot_ref
    # Remove Inf/NaN from cases where ref = 0
    boot_pr[!is.finite(boot_pr)] <- NA_real_

    # Prevalence/mean difference
    boot_pd <- boot_race - boot_ref

    pr_pd_list[[paste0(out_name, "_", rc, "_pr")]] <- tibble(
      outcome  = out_name,
      race     = race_labels[i],
      measure  = ifelse(out_name %in% c("mci_h", "ad_h", "adrd_h"),
                        "Prevalence Ratio", "Mean Ratio"),
      estimate = mean(boot_pr, na.rm = TRUE),
      lci      = quantile(boot_pr, 0.025, na.rm = TRUE),
      uci      = quantile(boot_pr, 0.975, na.rm = TRUE)
    )

    pr_pd_list[[paste0(out_name, "_", rc, "_pd")]] <- tibble(
      outcome  = out_name,
      race     = race_labels[i],
      measure  = ifelse(out_name %in% c("mci_h", "ad_h", "adrd_h"),
                        "Prevalence Difference", "Mean Difference"),
      estimate = mean(boot_pd, na.rm = TRUE),
      lci      = quantile(boot_pd, 0.025, na.rm = TRUE),
      uci      = quantile(boot_pd, 0.975, na.rm = TRUE)
    )
  }
}

pr_pd_results <- bind_rows(pr_pd_list)

# Format for display
pr_pd_table <- pr_pd_results %>%
  mutate(
    est_ci = sprintf("%.3f (%.3f, %.3f)", estimate, lci, uci),
    outcome_label = factor(outcome,
      levels = c("mci_h", "ad_h", "adrd_h", "moca_h", "cdrsb_h"),
      labels = c("MCI", "AD", "ADRD", "MoCA", "CDR-SB"))
  ) %>%
  arrange(outcome_label, measure, race)

saveRDS(pr_pd_results, file.path(results_dir, "09_pr_pd_results.rds"))

message("\nPrevalence/Mean Ratios and Differences (vs. White NH reference):")
message(paste(capture.output(print(pr_pd_table %>% select(outcome_label, race, measure, est_ci),
                                    n = Inf)), collapse = "\n"))

# ---- Forest plot: PR/PD ----
p_pr <- pr_pd_results %>%
  filter(grepl("Ratio", measure)) %>%
  mutate(
    outcome_label = factor(outcome,
      levels = c("mci_h", "ad_h", "adrd_h", "moca_h", "cdrsb_h"),
      labels = c("MCI", "AD", "ADRD", "MoCA", "CDR-SB")),
    race = factor(race, levels = c("Hispanic", "Asian NH", "NH Other"))
  ) %>%
  ggplot(aes(x = estimate, y = race, xmin = lci, xmax = uci, color = race)) +
  geom_pointrange(size = 0.8) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  facet_wrap(~outcome_label, scales = "free_x", ncol = 2) +
  labs(x = "Ratio vs. White NH (95% CI)", y = NULL,
       title = "AIPSW Race/Ethnicity Ratios (vs. White NH Reference)",
       color = NULL) +
  theme_minimal() +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold"))

ggsave(file.path(results_dir, "09_forest_pr.pdf"), p_pr, width = 10, height = 5)


# ============================================================================
# 5. SUMMARY RESULTS TABLE
# ============================================================================

results_table <- results %>%
  mutate(est_ci = sprintf("%.3f (%.3f, %.3f)", estimate, lci, uci)) %>%
  select(outcome, estimand, est_ci)

saveRDS(results_table, file.path(results_dir, "09_results_table.rds"))

message("Saved all results figures and tables to: ", results_dir)
