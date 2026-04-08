# 10_sensitivity_trimmed_weights.R
# Purpose: Sensitivity analysis comparing AIPSW estimates using untrimmed
#          vs 99th-percentile-trimmed stabilized IOSW weights.
#
# This addresses FLAG-P: Does trimming extreme weights change conclusions?
#
# Approach: Re-compute AIPSW point estimates across all 40 imputations
#           using trimmed weights, then compare to untrimmed estimates via
#           Rubin's rules for point estimates. No bootstrap (too expensive);
#           the question is whether trimming shifts the point estimates
#           materially, not whether the CIs change.
#
# INPUTS:
#   - processed_data_dir/06_outcome_predictions.rds
#   - results_dir/07_aipsw_results.rds (for comparison)
#
# OUTPUTS:
#   - results_dir/10_sensitivity_trimmed.rds
#   - results_dir/10_sensitivity_trimmed.csv
#
# ============================================================================

source("config.R")

# ============================================================================
# 1. LOAD DATA
# ============================================================================

outcome_data <- readRDS(file.path(processed_data_dir, "06_outcome_predictions.rds"))
aipsw_untrimmed <- readRDS(file.path(results_dir, "07_aipsw_results.rds"))

# ============================================================================
# 2. DEFINE ESTIMATOR (mirrors script 07 but adds trim_pct argument)
# ============================================================================

participation_formula <- adrc_h ~ race_eth4_h + male_h + age_h + education_h +
  marital_h + diabetes_h + hypertension_h + stroke_h +
  depression_h +
  race_eth4_h:male_h + race_eth4_h:age_h +
  race_eth4_h:education_h +
  race_eth4_h:diabetes_h + race_eth4_h:hypertension_h +
  race_eth4_h:depression_h

outcome_formula_base <- ~ race_eth4_h + male_h + age_h + education_h +
  marital_h + diabetes_h + hypertension_h + stroke_h +
  depression_h + mi_h + angina_h +
  race_eth4_h:age_h + race_eth4_h:male_h + race_eth4_h:education_h

compute_estimates_trimmed <- function(data, outcome, family = "binomial",
                                      two_part = FALSE, trim_pct = 0.99) {

  adrc  <- data %>% filter(adrc_h == 1L)
  brfss <- data %>% filter(adrc_h == 0L)
  Y     <- adrc[[outcome]]

  # Participation model
  n_brfss <- sum(data$adrc_h == 0L)
  brfss_wt_sum <- sum(data$brfss_sampwt_h[data$adrc_h == 0L])
  data$boot_wt <- ifelse(
    data$adrc_h == 1L, 1,
    data$brfss_sampwt_h / brfss_wt_sum * n_brfss
  )

  p_mod <- glm(participation_formula, family = binomial,
               data = data, weights = boot_wt)
  p_hat <- predict(p_mod, newdata = data, type = "response")

  p_adrc_sample <- sum(data$adrc_h == 1L) / nrow(data)
  w <- (1 - p_hat) / p_hat
  sw <- w * (p_adrc_sample / (1 - p_adrc_sample))
  adrc_sw <- sw[data$adrc_h == 1L]

  # ---- TRIM ----
  trim_threshold <- quantile(adrc_sw, trim_pct, na.rm = TRUE)
  adrc_sw_trimmed <- pmin(adrc_sw, trim_threshold)

  # Outcome model
  if (!two_part) {
    fam <- if (family == "binomial") binomial else gaussian
    o_mod <- glm(update(outcome_formula_base, as.formula(paste(outcome, "~ ."))),
                 family = fam, data = adrc)
    mu_hat <- predict(o_mod, newdata = data, type = "response")
  } else {
    adrc$cdrsb_pos <- as.integer(adrc$cdrsb_h > 0)
    p1_mod <- glm(update(outcome_formula_base, cdrsb_pos ~ .),
                  family = binomial, data = adrc)
    pred_p1 <- predict(p1_mod, newdata = data, type = "response")

    adrc_pos <- adrc %>% filter(cdrsb_h > 0)
    if (nrow(adrc_pos) >= 10) {
      p2_mod <- glm(update(outcome_formula_base, cdrsb_h ~ .),
                     family = Gamma(link = "log"), data = adrc_pos)
      pred_p2 <- pmin(predict(p2_mod, newdata = data, type = "response"), 18)
    } else {
      pred_p2 <- mean(adrc_pos$cdrsb_h, na.rm = TRUE)
    }
    mu_hat <- pred_p1 * pred_p2
  }

  mu_adrc  <- mu_hat[data$adrc_h == 1L]
  mu_brfss <- mu_hat[data$adrc_h == 0L]

  brfss_wt <- brfss$brfss_sampwt_h
  N_target <- sum(brfss_wt)

  # g-comp (same for trimmed and untrimmed — doesn't use IOSW)
  est_gcomp <- sum(mu_brfss * brfss_wt) / N_target

  # AIPSW with UNTRIMMED weights
  residuals_adrc <- Y - mu_adrc
  correction_untrimmed <- weighted.mean(residuals_adrc, adrc_sw)
  est_aipsw_untrimmed  <- est_gcomp + correction_untrimmed

  # AIPSW with TRIMMED weights
  correction_trimmed <- weighted.mean(residuals_adrc, adrc_sw_trimmed)
  est_aipsw_trimmed  <- est_gcomp + correction_trimmed

  # IOSW with trimmed and untrimmed
  est_iosw_untrimmed <- weighted.mean(Y, adrc_sw)
  est_iosw_trimmed   <- weighted.mean(Y, adrc_sw_trimmed)

  # Race-stratified AIPSW (trimmed)
  races <- c("White NH", "Hispanic", "Asian NH", "NH Other")
  race_trimmed <- setNames(numeric(length(races)),
                           paste0("aipsw_", gsub(" ", "_", tolower(races))))
  race_untrimmed <- race_trimmed

  adrc_race  <- data$race_eth4_h[data$adrc_h == 1L]
  brfss_race <- data$race_eth4_h[data$adrc_h == 0L]

  for (r in races) {
    r_in_adrc  <- which(adrc_race == r)
    r_in_brfss <- which(brfss_race == r)

    if (length(r_in_adrc) < 5 || length(r_in_brfss) < 5) {
      race_trimmed[paste0("aipsw_", gsub(" ", "_", tolower(r)))] <- NA_real_
      race_untrimmed[paste0("aipsw_", gsub(" ", "_", tolower(r)))] <- NA_real_
      next
    }

    wt_r    <- brfss_wt[r_in_brfss]
    gcomp_r <- sum(mu_brfss[r_in_brfss] * wt_r) / sum(wt_r)

    corr_r_untrimmed <- weighted.mean(residuals_adrc[r_in_adrc], adrc_sw[r_in_adrc])
    corr_r_trimmed   <- weighted.mean(residuals_adrc[r_in_adrc], adrc_sw_trimmed[r_in_adrc])

    key <- paste0("aipsw_", gsub(" ", "_", tolower(r)))
    race_untrimmed[key] <- gcomp_r + corr_r_untrimmed
    race_trimmed[key]   <- gcomp_r + corr_r_trimmed
  }

  list(
    untrimmed = c(aipsw = est_aipsw_untrimmed, iosw = est_iosw_untrimmed,
                  gcomp = est_gcomp, race_untrimmed),
    trimmed   = c(aipsw = est_aipsw_trimmed, iosw = est_iosw_trimmed,
                  gcomp = est_gcomp, race_trimmed),
    trim_threshold = trim_threshold,
    n_trimmed = sum(adrc_sw > trim_threshold),
    max_weight_before = max(adrc_sw),
    max_weight_after  = max(adrc_sw_trimmed)
  )
}


# ============================================================================
# 3. RUN ACROSS ALL IMPUTATIONS AND OUTCOMES
# ============================================================================

outcomes <- list(
  list(name = "mci_h",   family = "binomial", two_part = FALSE),
  list(name = "ad_h",    family = "binomial", two_part = FALSE),
  list(name = "adrd_h",  family = "binomial", two_part = FALSE),
  list(name = "moca_h",  family = "gaussian", two_part = FALSE),
  list(name = "cdrsb_h", family = "gaussian", two_part = TRUE)
)

imputations <- sort(unique(outcome_data$imp_h))
comparison_results <- list()

for (out in outcomes) {
  message("Processing: ", out$name)

  est_untrimmed <- list()
  est_trimmed   <- list()
  trim_info     <- list()

  for (imp in imputations) {
    imp_data <- outcome_data %>%
      filter(imp_h == imp) %>%
      filter(!is.na(.data[[out$name]]) | adrc_h == 0L)

    res <- tryCatch(
      compute_estimates_trimmed(imp_data, outcome = out$name,
                                family = out$family,
                                two_part = out$two_part,
                                trim_pct = 0.99),
      error = function(e) {
        message("  Error in imputation ", imp, ": ", e$message)
        NULL
      }
    )

    if (!is.null(res)) {
      est_untrimmed[[imp]] <- res$untrimmed
      est_trimmed[[imp]]   <- res$trimmed
      trim_info[[imp]]     <- c(threshold = res$trim_threshold,
                                n_trimmed = res$n_trimmed,
                                max_before = res$max_weight_before,
                                max_after  = res$max_weight_after)
    }
  }

  # Pool across imputations (simple mean — point estimates only)
  mat_untrimmed <- do.call(rbind, est_untrimmed)
  mat_trimmed   <- do.call(rbind, est_trimmed)
  mat_trim_info <- do.call(rbind, trim_info)

  est_names <- names(est_untrimmed[[1]])

  pooled_untrimmed <- colMeans(mat_untrimmed, na.rm = TRUE)
  pooled_trimmed   <- colMeans(mat_trimmed, na.rm = TRUE)
  mean_trim_info   <- colMeans(mat_trim_info, na.rm = TRUE)

  comparison_results[[out$name]] <- tibble(
    outcome = out$name,
    estimand = est_names,
    untrimmed = pooled_untrimmed,
    trimmed = pooled_trimmed,
    abs_diff = abs(pooled_trimmed - pooled_untrimmed),
    pct_diff = abs(pooled_trimmed - pooled_untrimmed) / abs(pooled_untrimmed) * 100,
    avg_trim_threshold = mean_trim_info["threshold"],
    avg_n_trimmed = mean_trim_info["n_trimmed"],
    avg_max_wt_before = mean_trim_info["max_before"],
    avg_max_wt_after = mean_trim_info["max_after"]
  )
}

sensitivity_table <- bind_rows(comparison_results)


# ============================================================================
# 4. DISPLAY AND SAVE
# ============================================================================

message("\n=== SENSITIVITY ANALYSIS: TRIMMED vs UNTRIMMED WEIGHTS ===")
message("99th percentile trimming of stabilized IOSW\n")

# Summary
summary_table <- sensitivity_table %>%
  filter(estimand == "aipsw") %>%
  transmute(
    Outcome = outcome,
    Untrimmed = sprintf("%.3f", untrimmed),
    Trimmed = sprintf("%.3f", trimmed),
    `Abs Diff` = sprintf("%.4f", abs_diff),
    `% Diff` = sprintf("%.1f%%", pct_diff),
    `Avg Max Wt (before)` = sprintf("%.1f", avg_max_wt_before),
    `Avg Max Wt (after)` = sprintf("%.1f", avg_max_wt_after),
    `Avg N Trimmed` = sprintf("%.1f", avg_n_trimmed)
  )

message("Overall AIPSW estimates:")
message(paste(capture.output(print(summary_table, n = Inf, width = 120)), collapse = "\n"))

# Race-stratified
message("\nRace-stratified AIPSW estimates:")
race_table <- sensitivity_table %>%
  filter(grepl("^aipsw_", estimand)) %>%
  transmute(
    Outcome = outcome,
    Race = estimand,
    Untrimmed = sprintf("%.3f", untrimmed),
    Trimmed = sprintf("%.3f", trimmed),
    `Abs Diff` = sprintf("%.4f", abs_diff),
    `% Diff` = sprintf("%.1f%%", pct_diff)
  )
message(paste(capture.output(print(race_table, n = Inf, width = 120)), collapse = "\n"))

# Conclusion
max_pct_overall <- max(sensitivity_table$pct_diff[sensitivity_table$estimand == "aipsw"], na.rm = TRUE)
max_pct_race <- max(sensitivity_table$pct_diff[grepl("^aipsw_", sensitivity_table$estimand)], na.rm = TRUE)

message(sprintf("\nMax %% change (overall): %.1f%%", max_pct_overall))
message(sprintf("Max %% change (race-specific): %.1f%%", max_pct_race))

if (max_pct_overall < 5) {
  message("CONCLUSION: Trimming has minimal impact on overall estimates (<5% change).")
} else {
  message("CAUTION: Trimming changes overall estimates by >5%. Consider investigating extreme weights.")
}

if (max_pct_race < 10) {
  message("CONCLUSION: Trimming has minimal impact on race-specific estimates (<10% change).")
} else {
  message("CAUTION: Trimming changes some race-specific estimates by >10%. Check small-N groups.")
}

saveRDS(sensitivity_table, file.path(results_dir, "10_sensitivity_trimmed.rds"))
write.csv(sensitivity_table, file.path(results_dir, "10_sensitivity_trimmed.csv"),
          row.names = FALSE)

message("\nSaved: 10_sensitivity_trimmed.rds")
message("Saved: 10_sensitivity_trimmed.csv")
