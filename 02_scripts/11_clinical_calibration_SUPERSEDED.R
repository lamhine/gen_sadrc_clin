# 11_clinical_calibration.R
# Purpose: Secondary analysis -- Correct for clinical enrichment in ADRC
#          using external ADRD prevalence from Matthews et al. (2019).
#
# The ADRC over-represents cognitively impaired participants (~52% ADRD)
# relative to the general 65+ population (~10%). The primary AIPSW analysis
# (script 07) corrects for demographic selection but not clinical enrichment.
#
# This script implements a model-based calibration:
#   1. Fit P(ADRD | race, age, sex) in ADRC
#   2. Recalibrate intercept so predicted prevalence matches Matthews et al.
#      in the BRFSS target population
#   3. Compute clinical enrichment weights for ADRC participants
#   4. Re-estimate ADRD prevalence with combined demographic + clinical weights
#
# Two calibrated estimators:
#   A. G-computation (calibrated): predict P_pop(ADRD|X) in BRFSS, take
#      survey-weighted mean. Clean but ignores ADRC residual information.
#   B. IOSW × clinical weight: multiply demographic IOSW by outcome-based
#      calibration weight. Uses ADRC outcomes directly.
#
# Reference: Matthews KA et al. (2019). Racial and ethnic estimates of
#   Alzheimer's disease and related dementias in the United States (2015-2060)
#   in adults aged ≥65 years. Alzheimers Dement. 15(1):17-24. PMC6333531.
#
# INPUTS:
#   - processed_data_dir/05_imputed_with_weights.rds
#
# OUTPUTS:
#   - results_dir/11_calibrated_adrd.rds
#   - results_dir/11_calibration_comparison.pdf
#
# ============================================================================

source("config.R")
library(survey)


# ============================================================================
# 1. MATTHEWS ET AL. (2019) EXTERNAL PREVALENCE TARGETS
# ============================================================================
# Table 2 / Figure 1: ADRD prevalence by demographic group (2015 estimates)
# Based on Medicare FFS claims using CMS Chronic Conditions Warehouse algorithm.
#
# Note: Matthews uses national Medicare FFS data. Our target is 3-county Bay
# Area adults 65+. We assume national race/age/sex-specific prevalences are
# reasonable approximations for this catchment area.

matthews_race <- tribble(
  ~race,       ~prev_matthews,
  "White NH",   0.103,
  "Hispanic",   0.122,
  "Asian NH",   0.084,   # "Asian/Pacific Islander" in Matthews
  "NH Other",   0.138    # Black NH (13.8%) is primary component of our NH Other
)

matthews_age <- tribble(
  ~age_group, ~prev_matthews,
  "65-74",     0.036,
  "75-84",     0.136,
  "85+",       0.346
)

matthews_sex <- tribble(
  ~sex,     ~prev_matthews,
  "Female",  0.122,
  "Male",    0.086
)


# ============================================================================
# 2. LOAD DATA
# ============================================================================

imputed_weighted <- readRDS(file.path(processed_data_dir, "05_imputed_with_weights.rds"))
imputations <- sort(unique(imputed_weighted$imp_h))


# ============================================================================
# 3. MODEL-BASED CALIBRATION (LOOP OVER IMPUTATIONS)
# ============================================================================
#
# Strategy:
#   a) Fit logistic P(ADRD | race, age, sex) in ADRC -- the "ADRC model"
#   b) The population model uses the same coefficients but shifts the intercept
#      by δ so that E_target[expit(η + δ)] = π_Matthews
#   c) Clinical enrichment weight for each ADRC participant:
#        ADRD=1: w_clin = P_pop(X) / P_adrc(X)     (down-weight cases)
#        ADRD=0: w_clin = (1-P_pop(X)) / (1-P_adrc(X))  (up-weight controls)
#   d) Combined weight = IOSW_demographic × w_clinical

calibrated_list <- list()

for (imp in imputations) {

  imp_data <- imputed_weighted %>% filter(imp_h == imp)
  adrc  <- imp_data %>% filter(adrc_h == 1L)
  brfss <- imp_data %>% filter(adrc_h == 0L)

  # ---- 3a. Fit ADRD model in ADRC ----
  # Parsimonious model: race + age + sex (matches Matthews stratification)
  adrc_mod <- glm(adrd_h ~ race_eth4_h + age_h + male_h,
                  family = binomial, data = adrc)

  # Linear predictor for all participants
  eta <- predict(adrc_mod, newdata = imp_data, type = "link")
  p_adrc_all <- plogis(eta)

  # ---- 3b. Compute target prevalence ----
  # Weight Matthews race-specific prevalences by BRFSS race distribution
  brfss_race_wt <- brfss %>%
    group_by(race_eth4_h) %>%
    summarise(wt = sum(brfss_sampwt_h), .groups = "drop") %>%
    mutate(prop = wt / sum(wt))

  target_prev <- brfss_race_wt %>%
    left_join(matthews_race, by = c("race_eth4_h" = "race")) %>%
    summarise(overall = sum(prop * prev_matthews, na.rm = TRUE)) %>%
    pull(overall)

  # ---- 3c. Find intercept shift δ via uniroot ----
  # Solve: E_BRFSS[expit(η_brfss + δ)] = target_prev
  eta_brfss <- eta[imp_data$adrc_h == 0L]
  wt_brfss  <- brfss$brfss_sampwt_h

  cal_fn <- function(delta) {
    weighted.mean(plogis(eta_brfss + delta), wt_brfss) - target_prev
  }

  delta <- uniroot(cal_fn, interval = c(-15, 5), tol = 1e-10)$root

  # ---- 3d. Calibrated predictions ----
  p_pop_all  <- plogis(eta + delta)
  p_pop_adrc <- p_pop_all[imp_data$adrc_h == 1L]
  p_adrc_adrc <- p_adrc_all[imp_data$adrc_h == 1L]

  # Verify calibration hit target
  calib_check <- weighted.mean(p_pop_all[imp_data$adrc_h == 0L], wt_brfss)

  # ---- 3e. Clinical enrichment weights (ADRC only) ----
  Y <- adrc$adrd_h

  w_clinical <- ifelse(
    Y == 1L,
    p_pop_adrc / p_adrc_adrc,             # < 1: down-weight cases
    (1 - p_pop_adrc) / (1 - p_adrc_adrc)  # > 1: up-weight controls
  )

  # Guard against numerical edge cases
  w_clinical[!is.finite(w_clinical)] <- 1
  w_clinical <- pmax(w_clinical, 0.01)

  # ---- 3f. Compute estimates ----
  iosw <- adrc$sw3_final
  w_combined <- iosw * w_clinical

  # Overall estimates
  est_unweighted       <- mean(Y, na.rm = TRUE)
  est_iosw_only        <- weighted.mean(Y, iosw)
  est_calibrated_iosw  <- weighted.mean(Y, w_combined)
  est_gcomp_calibrated <- weighted.mean(p_pop_all[imp_data$adrc_h == 0L], wt_brfss)

  # Race-stratified estimates (g-computation with calibrated model)
  gcomp_race <- setNames(numeric(4), c("White NH", "Hispanic", "Asian NH", "NH Other"))
  calib_iosw_race <- gcomp_race

  for (r in names(gcomp_race)) {
    # G-comp: calibrated predictions in BRFSS target
    idx_b <- which(brfss$race_eth4_h == r)
    if (length(idx_b) >= 5) {
      gcomp_race[r] <- weighted.mean(
        p_pop_all[imp_data$adrc_h == 0L][idx_b],
        wt_brfss[idx_b]
      )
    } else {
      gcomp_race[r] <- NA_real_
    }

    # IOSW × clinical: weighted ADRC outcomes
    idx_a <- which(adrc$race_eth4_h == r)
    if (length(idx_a) >= 5) {
      calib_iosw_race[r] <- weighted.mean(Y[idx_a], w_combined[idx_a])
    } else {
      calib_iosw_race[r] <- NA_real_
    }
  }

  calibrated_list[[imp]] <- tibble(
    imp_h = imp,
    delta = delta,
    target_prev = target_prev,
    calib_check = calib_check,
    est_unweighted       = est_unweighted,
    est_iosw_only        = est_iosw_only,
    est_calibrated_iosw  = est_calibrated_iosw,
    est_gcomp_calibrated = est_gcomp_calibrated,
    # Race-stratified g-comp (calibrated)
    gcomp_white_nh = gcomp_race["White NH"],
    gcomp_hispanic = gcomp_race["Hispanic"],
    gcomp_asian_nh = gcomp_race["Asian NH"],
    gcomp_nh_other = gcomp_race["NH Other"],
    # Race-stratified IOSW × clinical
    calib_white_nh = calib_iosw_race["White NH"],
    calib_hispanic = calib_iosw_race["Hispanic"],
    calib_asian_nh = calib_iosw_race["Asian NH"],
    calib_nh_other = calib_iosw_race["NH Other"],
    # Diagnostics
    w_clin_mean_case = mean(w_clinical[Y == 1L]),
    w_clin_mean_ctrl = mean(w_clinical[Y == 0L]),
    w_clin_max = max(w_clinical),
    enrichment_ratio = est_iosw_only / est_gcomp_calibrated
  )

  if (imp %% 10 == 0 || imp == 1) {
    message(sprintf(
      "Imp %2d: δ=%.2f  target=%.3f  unwtd=%.3f  IOSW=%.3f  calib_IOSW=%.3f  gcomp=%.3f",
      imp, delta, target_prev, est_unweighted, est_iosw_only,
      est_calibrated_iosw, est_gcomp_calibrated
    ))
  }
}

calibrated_results <- bind_rows(calibrated_list)


# ============================================================================
# 4. COMBINE ACROSS IMPUTATIONS (RUBIN'S RULES)
# ============================================================================
# Note: Without within-imputation bootstrap, CIs reflect only between-
# imputation variance. These will be narrower than full bootstrap CIs.
# For final analysis, wrap in bootstrap (as in script 07).

estimate_cols <- c(
  "est_unweighted", "est_iosw_only", "est_calibrated_iosw", "est_gcomp_calibrated",
  "gcomp_white_nh", "gcomp_hispanic", "gcomp_asian_nh", "gcomp_nh_other",
  "calib_white_nh", "calib_hispanic", "calib_asian_nh", "calib_nh_other"
)

M <- length(imputations)

combined <- map_dfr(estimate_cols, function(col) {
  vals <- calibrated_results[[col]]
  theta_bar <- mean(vals, na.rm = TRUE)
  B <- var(vals, na.rm = TRUE)
  # Without within-imputation variance, use (1 + 1/M) * B as lower bound
  Total_var <- (1 + 1/M) * B
  se <- sqrt(Total_var)
  tibble(
    estimand = col,
    estimate = theta_bar,
    se = se,
    lci = theta_bar - 1.96 * se,
    uci = theta_bar + 1.96 * se
  )
})


# ============================================================================
# 5. FORMAT AND DISPLAY
# ============================================================================

labels_map <- c(
  est_unweighted       = "Unweighted ADRC",
  est_iosw_only        = "IOSW only (demographic adjustment)",
  est_calibrated_iosw  = "IOSW x clinical calibration",
  est_gcomp_calibrated = "G-comp (calibrated model)",
  gcomp_white_nh       = "  G-comp calibrated: White NH",
  gcomp_hispanic       = "  G-comp calibrated: Hispanic",
  gcomp_asian_nh       = "  G-comp calibrated: Asian NH",
  gcomp_nh_other       = "  G-comp calibrated: NH Other",
  calib_white_nh       = "  Calib IOSW: White NH",
  calib_hispanic       = "  Calib IOSW: Hispanic",
  calib_asian_nh       = "  Calib IOSW: Asian NH",
  calib_nh_other       = "  Calib IOSW: NH Other"
)

comparison <- combined %>%
  mutate(
    est_ci = sprintf("%.3f (%.3f, %.3f)", estimate, lci, uci),
    label = labels_map[estimand],
    label = factor(label, levels = labels_map)
  )

message("\n", strrep("=", 60))
message("ADRD PREVALENCE: CLINICAL CALIBRATION RESULTS")
message(strrep("=", 60))
message("\nMatthews et al. (2019) external targets:")
message("  White: 10.3%  |  Hispanic: 12.2%  |  Asian/PI: 8.4%  |  Black: 13.8%")
message("  Weighted target for this catchment: ",
        round(mean(calibrated_results$target_prev) * 100, 1), "%")
message("\nMean intercept shift δ = ",
        round(mean(calibrated_results$delta), 2),
        " (ADRC model intercept shifted down to match population prevalence)")
message("\nMean clinical weights:")
message("  Cases (ADRD=1):    ", round(mean(calibrated_results$w_clin_mean_case), 3),
        "  (down-weighted)")
message("  Controls (ADRD=0): ", round(mean(calibrated_results$w_clin_mean_ctrl), 3),
        "  (up-weighted)")
message("\nEstimated ADRD Prevalence (65+):")
message(paste(capture.output(
  print(comparison %>% select(label, est_ci), n = Inf)
), collapse = "\n"))


# ============================================================================
# 6. COMPARISON FIGURE
# ============================================================================

# Panel A: Overall estimator comparison
overall_df <- combined %>%
  filter(estimand %in% c("est_unweighted", "est_iosw_only",
                          "est_calibrated_iosw", "est_gcomp_calibrated")) %>%
  mutate(
    label = factor(
      c("Unweighted ADRC", "IOSW (demographic only)",
        "IOSW x clinical calibration", "G-comp (calibrated)"),
      levels = rev(c("Unweighted ADRC", "IOSW (demographic only)",
                      "IOSW x clinical calibration", "G-comp (calibrated)"))
    )
  )

p_overall <- ggplot(overall_df, aes(x = estimate, y = label,
                                     xmin = lci, xmax = uci)) +
  geom_pointrange(size = 0.8, color = "steelblue") +
  geom_vline(xintercept = mean(calibrated_results$target_prev),
             linetype = "dashed", color = "red", alpha = 0.7) +
  annotate("text",
           x = mean(calibrated_results$target_prev) + 0.01,
           y = 0.6, label = "Matthews target",
           color = "red", hjust = 0, size = 3) +
  scale_x_continuous(labels = scales::percent_format(),
                     limits = c(0, 0.65)) +
  labs(x = "ADRD Prevalence (65+)", y = NULL,
       title = "ADRD Prevalence Estimates: Effect of Clinical Calibration") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12))

# Panel B: Race-stratified (g-comp calibrated vs Matthews targets)
race_df <- combined %>%
  filter(grepl("^gcomp_", estimand)) %>%
  mutate(
    race = c("White NH", "Hispanic", "Asian NH", "NH Other"),
    race = factor(race, levels = c("White NH", "Hispanic", "Asian NH", "NH Other"))
  ) %>%
  left_join(matthews_race, by = "race")

p_race <- ggplot(race_df, aes(x = estimate, y = race,
                               xmin = lci, xmax = uci)) +
  geom_pointrange(size = 0.8, color = "steelblue") +
  geom_point(aes(x = prev_matthews), shape = 4, size = 3, color = "red") +
  scale_x_continuous(labels = scales::percent_format(),
                     limits = c(0, 0.25)) +
  labs(x = "ADRD Prevalence (65+)", y = NULL,
       title = "Race/Ethnicity-Stratified (G-comp Calibrated)",
       subtitle = "Blue = calibrated estimate, Red X = Matthews et al. target") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12))

# Combine panels
p_combined <- gridExtra::grid.arrange(p_overall, p_race, ncol = 1, heights = c(1, 1))

ggsave(file.path(results_dir, "11_calibration_comparison.pdf"),
       p_combined, width = 9, height = 8)


# ============================================================================
# 7. SAVE
# ============================================================================

saveRDS(list(
  calibrated_by_imp = calibrated_results,
  combined          = combined,
  comparison        = comparison,
  matthews_targets  = matthews_race,
  diagnostics = tibble(
    mean_delta           = mean(calibrated_results$delta),
    mean_target_prev     = mean(calibrated_results$target_prev),
    mean_enrichment_ratio = mean(calibrated_results$enrichment_ratio),
    mean_w_clin_case     = mean(calibrated_results$w_clin_mean_case),
    mean_w_clin_ctrl     = mean(calibrated_results$w_clin_mean_ctrl)
  )
), file.path(results_dir, "11_calibrated_adrd.rds"))

message("\nSaved: 11_calibrated_adrd.rds")
message("Saved: 11_calibration_comparison.pdf")
