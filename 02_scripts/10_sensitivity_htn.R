# 10_sensitivity_htn.R
# Purpose: Sensitivity analysis -- Compare AIPSW estimates WITH vs WITHOUT
#          hypertension in the participation model.
#
# RATIONALE (Decision Log D05/Flag I):
#   Hypertension has ~76% structural missingness in BRFSS because the HTN
#   question (RFHYPE6/RFHYPE5) is only asked in odd-numbered survey years.
#   Even after multiple imputation with survey_year as an auxiliary predictor,
#   there is concern about the quality of imputed HTN values.
#
#   This sensitivity analysis fits the full pipeline WITHOUT hypertension
#   and compares results to the main analysis to assess the impact.
#
# APPROACH:
#   1. Re-fit participation model (script 05 Model 3) dropping hypertension_h
#      and its race interaction
#   2. Re-fit outcome models (script 06) dropping hypertension_h
#   3. Compute AIPSW point estimates (no bootstrap -- just compare direction)
#   4. Display side-by-side comparison with main results
#
# INPUTS:
#   - processed_data_dir/04_imputed_stacked.rds
#   - results_dir/07_aipsw_results.rds (main analysis for comparison)
#
# OUTPUTS:
#   - results_dir/10_sensitivity_htn_comparison.rds
#   - results_dir/10_sensitivity_htn_comparison.pdf
#
# ============================================================================

source("config.R")
library(survey)

# ============================================================================
# 1. LOAD DATA
# ============================================================================

imputed <- readRDS(file.path(processed_data_dir, "04_imputed_stacked.rds"))
main_results <- readRDS(file.path(results_dir, "07_aipsw_results.rds"))


# ============================================================================
# 2. PARTICIPATION MODEL WITHOUT HYPERTENSION
# ============================================================================
# Same as script 05 Model 3 but dropping hypertension_h and race:hypertension

participation_formula_no_htn <- adrc_h ~ race_eth4_h + male_h + age_h +
  education_h + marital_h + diabetes_h + stroke_h + depression_h +
  race_eth4_h:male_h + race_eth4_h:age_h +
  race_eth4_h:education_h + race_eth4_h:diabetes_h +
  race_eth4_h:depression_h

# Outcome model formula also without hypertension
outcome_formula_no_htn <- ~ race_eth4_h + male_h + age_h + education_h +
  marital_h + diabetes_h + stroke_h + depression_h + mi_h + angina_h +
  race_eth4_h:age_h + race_eth4_h:male_h + race_eth4_h:education_h

imputations <- sort(unique(imputed$imp_h))
weighted_list <- list()

message("Fitting participation models WITHOUT hypertension...")

for (imp in imputations) {
  imp_data <- imputed %>% filter(imp_h == imp)

  combined_design <- svydesign(ids = ~1, weights = ~brfss_sampwt_h,
                               data = imp_data)

  m_no_htn <- svyglm(participation_formula_no_htn,
                      family = quasibinomial,
                      design = combined_design)

  imp_data$p_nohtn <- as.numeric(predict(m_no_htn, newdata = imp_data,
                                          type = "response"))

  # Stabilized IOSW
  p_adrc_pop <- sum(imp_data$brfss_sampwt_h[imp_data$adrc_h == 1L]) /
                sum(imp_data$brfss_sampwt_h)

  imp_data <- imp_data %>%
    mutate(
      w_nohtn  = (1 - p_nohtn) / p_nohtn,
      sw_nohtn = w_nohtn * (p_adrc_pop / (1 - p_adrc_pop)),
      sw_nohtn_final = case_when(
        adrc_h == 1L ~ sw_nohtn,
        adrc_h == 0L ~ 1
      )
    )

  weighted_list[[imp]] <- imp_data
}

imputed_nohtn <- bind_rows(weighted_list)


# ============================================================================
# 3. OUTCOME MODELS WITHOUT HYPERTENSION
# ============================================================================

message("Fitting outcome models WITHOUT hypertension...")

pred_list <- list()

for (imp in imputations) {
  imp_data <- imputed_nohtn %>% filter(imp_h == imp)
  adrc <- imp_data %>% filter(adrc_h == 1L)

  # MCI
  fit_mci <- glm(update(outcome_formula_no_htn, mci_h ~ .),
                 family = binomial, data = adrc)
  imp_data$pred_mci_nh <- predict(fit_mci, newdata = imp_data, type = "response")

  # AD
  fit_ad <- glm(update(outcome_formula_no_htn, ad_h ~ .),
                family = binomial, data = adrc)
  imp_data$pred_ad_nh <- predict(fit_ad, newdata = imp_data, type = "response")

  # MoCA
  fit_moca <- lm(update(outcome_formula_no_htn, moca_h ~ .), data = adrc)
  imp_data$pred_moca_nh <- predict(fit_moca, newdata = imp_data)

  # CDR-SB (two-part)
  adrc$cdrsb_pos <- as.integer(adrc$cdrsb_h > 0)
  fit_p1 <- glm(update(outcome_formula_no_htn, cdrsb_pos ~ .),
                family = binomial, data = adrc)
  pred_p1 <- predict(fit_p1, newdata = imp_data, type = "response")

  adrc_pos <- adrc %>% filter(cdrsb_h > 0)
  if (nrow(adrc_pos) >= 10) {
    fit_p2 <- lm(update(outcome_formula_no_htn, cdrsb_h ~ .), data = adrc_pos)
    pred_p2 <- pmax(predict(fit_p2, newdata = imp_data), 0)
  } else {
    pred_p2 <- mean(adrc_pos$cdrsb_h, na.rm = TRUE)
  }
  imp_data$pred_cdrsb_nh <- pred_p1 * pred_p2

  pred_list[[imp]] <- imp_data
}

outcome_nohtn <- bind_rows(pred_list)


# ============================================================================
# 4. COMPUTE AIPSW ESTIMATES (point estimates across imputations)
# ============================================================================

message("Computing AIPSW estimates without hypertension...")

outcomes <- list(
  list(name = "mci_h",   pred = "pred_mci_nh",   family = "binomial"),
  list(name = "ad_h",    pred = "pred_ad_nh",     family = "binomial"),
  list(name = "moca_h",  pred = "pred_moca_nh",   family = "gaussian"),
  list(name = "cdrsb_h", pred = "pred_cdrsb_nh",  family = "gaussian")
)

nohtn_results <- list()

for (out in outcomes) {
  imp_estimates <- numeric(length(imputations))

  for (i in seq_along(imputations)) {
    imp <- imputations[i]
    imp_data <- outcome_nohtn %>% filter(imp_h == imp)

    adrc  <- imp_data %>% filter(adrc_h == 1L)
    brfss <- imp_data %>% filter(adrc_h == 0L)

    Y        <- adrc[[out$name]]
    mu_adrc  <- adrc[[out$pred]]
    mu_brfss <- brfss[[out$pred]]
    brfss_wt <- brfss$brfss_sampwt_h
    adrc_sw  <- adrc$sw_nohtn

    # G-computation in target
    est_gcomp <- sum(mu_brfss * brfss_wt) / sum(brfss_wt)

    # IOSW-weighted residual correction
    residuals_adrc <- Y - mu_adrc
    correction <- weighted.mean(residuals_adrc, adrc_sw)

    imp_estimates[i] <- est_gcomp + correction
  }

  nohtn_results[[out$name]] <- tibble(
    outcome  = out$name,
    estimand = "aipsw_no_htn",
    estimate = mean(imp_estimates, na.rm = TRUE)
  )
}

nohtn_df <- bind_rows(nohtn_results)


# ============================================================================
# 5. COMPARE WITH MAIN RESULTS
# ============================================================================

main_aipsw <- main_results %>%
  filter(estimand == "aipsw") %>%
  select(outcome, estimate_main = estimate, lci_main = lci, uci_main = uci)

comparison <- nohtn_df %>%
  left_join(main_aipsw, by = "outcome") %>%
  mutate(
    abs_diff    = abs(estimate - estimate_main),
    pct_diff    = abs_diff / abs(estimate_main) * 100,
    within_ci   = estimate >= lci_main & estimate <= uci_main,
    outcome_label = factor(outcome,
      levels = c("mci_h", "ad_h", "moca_h", "cdrsb_h"),
      labels = c("MCI Prevalence", "AD Prevalence",
                  "MoCA Score", "CDR-SB Score"))
  )

saveRDS(comparison, file.path(results_dir, "10_sensitivity_htn_comparison.rds"))

message("\n========== HTN SENSITIVITY ANALYSIS ==========")
message("Comparison: AIPSW with HTN (main) vs. without HTN (sensitivity)\n")
for (i in seq_len(nrow(comparison))) {
  r <- comparison[i, ]
  message(sprintf("%-16s  Main: %.3f (%.3f, %.3f)  |  No HTN: %.3f  |  Diff: %.3f (%.1f%%)  |  Within main CI: %s",
    r$outcome_label, r$estimate_main, r$lci_main, r$uci_main,
    r$estimate, r$abs_diff, r$pct_diff,
    ifelse(r$within_ci, "YES", "NO")))
}


# ============================================================================
# 6. VISUALIZATION
# ============================================================================

plot_data <- bind_rows(
  main_results %>%
    filter(estimand == "aipsw") %>%
    mutate(model = "With HTN (main)"),
  nohtn_df %>%
    mutate(lci = NA_real_, uci = NA_real_, model = "Without HTN (sensitivity)")
) %>%
  mutate(
    outcome_label = factor(outcome,
      levels = c("mci_h", "ad_h", "moca_h", "cdrsb_h"),
      labels = c("MCI", "AD", "MoCA", "CDR-SB")),
    model = factor(model, levels = c("With HTN (main)", "Without HTN (sensitivity)"))
  )

p_htn <- ggplot(plot_data, aes(x = estimate, y = model,
                                xmin = lci, xmax = uci,
                                color = model, shape = model)) +
  geom_pointrange(size = 0.8, na.rm = TRUE) +
  facet_wrap(~outcome_label, scales = "free_x", ncol = 2) +
  scale_color_manual(values = c("With HTN (main)" = "steelblue",
                                 "Without HTN (sensitivity)" = "darkorange")) +
  scale_shape_manual(values = c("With HTN (main)" = 16,
                                 "Without HTN (sensitivity)" = 17)) +
  labs(x = "AIPSW Estimate", y = NULL,
       title = "Sensitivity Analysis: Effect of Hypertension in Models",
       subtitle = "Point estimates with vs. without hypertension (76% structurally missing in BRFSS)",
       color = NULL, shape = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold"))

ggsave(file.path(results_dir, "10_sensitivity_htn_comparison.pdf"),
       p_htn, width = 10, height = 6)

message("\nSaved: 10_sensitivity_htn_comparison.rds and .pdf")
