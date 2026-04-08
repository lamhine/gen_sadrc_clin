# 05_weight_development.R
# Purpose: Develop stabilized inverse odds of selection weights (sIOSW) to
#          balance ADRC and BRFSS cohorts. Iterative propensity score modeling
#          with race/ethnicity interactions, assessed via SMDs.
#
# This is the PARTICIPATION MODEL component of the doubly robust estimator.
# Weights from this script feed into both the IOSW-only sensitivity analysis
# and the AIPSW primary analysis.
#
# PARTICIPATION MODEL COVARIATES (from authoritative harmonization map):
#   age_h, male_h, race_eth4_h, education_h, marital_h,
#   diabetes_h, depression_h, mi_h, angina_h, stroke_h, hypertension_h
#
# Follows: Hayes-Larson et al. (2022) KHANDLE-weighting approach
#
# INPUTS:
#   - processed_data_dir/04_imputed_stacked.rds
#
# OUTPUTS:
#   - processed_data_dir/05_imputed_with_weights.rds
#   - results_dir/05_covbal_unweighted.rds
#   - results_dir/05_covbal_weighted.rds
#
# ============================================================================

source("config.R")
library(survey)


# ============================================================================
# 1. LOAD DATA
# ============================================================================

imputed <- readRDS(file.path(processed_data_dir, "04_imputed_stacked.rds"))


# ============================================================================
# 2. CUSTOM COVARIATE BALANCE FUNCTION
# ============================================================================
# Computes standardized mean differences (SMD) stratified by race/ethnicity.
# SMD = (mean_ADRC - mean_BRFSS) / SD_BRFSS
# where ADRC means use weight_var (if supplied) and BRFSS means use brfss_sampwt_h.
# Target: |SMD| < 0.25.
#
# For factor covariates, SMDs are computed for each level (dummy-coded).

Checkcovbal <- function(data, weight_var = NULL,
                        covariates, race_var = "race_eth4_h") {

  races <- levels(factor(data[[race_var]]))
  bal_results <- list()

  for (r in races) {
    subset_r <- data %>% filter(.data[[race_var]] == r)
    if (nrow(subset_r) < 20) next

    adrc_r  <- subset_r %>% filter(adrc_h == 1L)
    brfss_r <- subset_r %>% filter(adrc_h == 0L)

    adrc_w  <- if (!is.null(weight_var)) adrc_r[[weight_var]] else rep(1, nrow(adrc_r))
    brfss_w <- brfss_r$brfss_sampwt_h

    smd_list <- list()
    for (v in covariates) {
      if (is.factor(subset_r[[v]])) {
        # Compute SMD for each factor level (dummy coding)
        for (lev in levels(subset_r[[v]])) {
          dummy_a <- as.numeric(adrc_r[[v]] == lev)
          dummy_b <- as.numeric(brfss_r[[v]] == lev)
          mean_a  <- weighted.mean(dummy_a, adrc_w, na.rm = TRUE)
          mean_b  <- weighted.mean(dummy_b, brfss_w, na.rm = TRUE)
          sd_b    <- sqrt(weighted.mean((dummy_b - mean_b)^2, brfss_w, na.rm = TRUE))
          smd_val <- if (sd_b > 0) (mean_a - mean_b) / sd_b else NA_real_
          smd_list[[paste0(v, ":", lev)]] <- tibble(
            race = r, variable = paste0(v, ": ", lev),
            mean_adrc = mean_a, mean_brfss = mean_b, smd = smd_val
          )
        }
      } else {
        vals_a <- adrc_r[[v]]
        vals_b <- brfss_r[[v]]
        mean_a <- weighted.mean(vals_a, adrc_w, na.rm = TRUE)
        mean_b <- weighted.mean(vals_b, brfss_w, na.rm = TRUE)
        sd_b   <- sqrt(weighted.mean((vals_b - mean_b)^2, brfss_w, na.rm = TRUE))
        smd_val <- if (sd_b > 0) (mean_a - mean_b) / sd_b else NA_real_
        smd_list[[v]] <- tibble(
          race = r, variable = v,
          mean_adrc = mean_a, mean_brfss = mean_b, smd = smd_val
        )
      }
    }

    bal_results[[r]] <- bind_rows(smd_list)
  }

  bind_rows(bal_results)
}


# ============================================================================
# 3. ITERATIVE PROPENSITY SCORE MODELS
# ============================================================================
# Three-stage iterative approach (adapted from KHANDLE):
#   Model 1: Core demographics (race, sex, age, education)
#   Model 2: + marital status and medical comorbidities
#   Model 3: + depression (FINAL)
#
# DECISION D06: MI and angina removed from participation model.
# Post-weighting balance diagnostics showed MI and angina flipping from
# small negative SMDs (-0.07, -0.13) to problematic positive SMDs (+0.32,
# +0.20) due to a few heavily-weighted individuals with these conditions.
# Both have low prevalence (~2-4%) and cause weight instability.
# MI and angina are retained in the OUTCOME model (script 06) where they
# contribute to predicting cognitive outcomes.
#
# Race/ethnicity interactions are included at each stage.

covariates_m1 <- c("race_eth4_h", "male_h", "age_h", "education_h")

covariates_m2 <- c(covariates_m1,
                   "marital_h", "diabetes_h", "hypertension_h", "stroke_h")

covariates_m3 <- c(covariates_m2,
                   "depression_h")
# Note: mi_h and angina_h deliberately excluded — see decision log D06

# Process each imputation
imputations <- sort(unique(imputed$imp_h))
weighted_list <- list()

for (imp in imputations) {

  imp_data <- imputed %>% filter(imp_h == imp)

  # ---- Create combined survey design ----
  # BRFSS respondents carry their survey design weights (_LLCPWT).
  # ADRC participants are treated as a census of themselves (weight = 1).
  # svyglm() properly handles these as probability/sampling weights,
  # unlike glm(weights=) which treats them as frequency weights.
  combined_design <- svydesign(
    ids     = ~1,
    weights = ~brfss_sampwt_h,
    data    = imp_data
  )

  # ---- Model 1: Core demographics ----
  m1 <- svyglm(adrc_h ~ race_eth4_h + male_h + age_h + education_h +
                  race_eth4_h:male_h + race_eth4_h:age_h +
                  race_eth4_h:education_h,
                family = quasibinomial,
                design = combined_design)

  imp_data$p1 <- as.numeric(predict(m1, newdata = imp_data, type = "response"))

  # ---- Model 2: + marital & comorbidities ----
  m2 <- svyglm(adrc_h ~ race_eth4_h + male_h + age_h + education_h +
                  marital_h + diabetes_h + hypertension_h + stroke_h +
                  race_eth4_h:male_h + race_eth4_h:age_h +
                  race_eth4_h:education_h +
                  race_eth4_h:diabetes_h + race_eth4_h:hypertension_h,
                family = quasibinomial,
                design = combined_design)

  imp_data$p2 <- as.numeric(predict(m2, newdata = imp_data, type = "response"))

  # ---- Model 3: + depression (FINAL) ----
  # Note: MI and angina excluded from participation model (decision D06)
  m3 <- svyglm(adrc_h ~ race_eth4_h + male_h + age_h + education_h +
                  marital_h + diabetes_h + hypertension_h + stroke_h +
                  depression_h +
                  race_eth4_h:male_h + race_eth4_h:age_h +
                  race_eth4_h:education_h +
                  race_eth4_h:diabetes_h + race_eth4_h:hypertension_h +
                  race_eth4_h:depression_h,
                family = quasibinomial,
                design = combined_design)

  imp_data$p3 <- as.numeric(predict(m3, newdata = imp_data, type = "response"))

  # ---- Compute stabilized IOSW from Model 3 ----
  # Use the SURVEY-WEIGHTED marginal P(ADRC) for stabilization.
  # This is the population-level proportion, not the unweighted sample proportion.
  p_adrc_pop <- sum(imp_data$brfss_sampwt_h[imp_data$adrc_h == 1L]) /
                sum(imp_data$brfss_sampwt_h)

  imp_data <- imp_data %>%
    mutate(
      # Inverse odds weight: P(BRFSS|X) / P(ADRC|X)
      w3  = (1 - p3) / p3,
      # Stabilized: multiply by marginal odds of ADRC participation
      sw3 = w3 * (p_adrc_pop / (1 - p_adrc_pop))
    )

  # Trim at 99th percentile of ADRC weights (sensitivity analysis)
  adrc_sw3 <- imp_data$sw3[imp_data$adrc_h == 1L]
  trim_99  <- quantile(adrc_sw3, 0.99, na.rm = TRUE)

  imp_data <- imp_data %>%
    mutate(
      sw3_trimmed = pmin(sw3, trim_99),
      # Primary analytic weight (untrimmed)
      sw3_final = case_when(
        adrc_h == 1L ~ sw3,
        adrc_h == 0L ~ 1
      ),
      # Trimmed weight (sensitivity analysis)
      sw3_trimmed_final = case_when(
        adrc_h == 1L ~ sw3_trimmed,
        adrc_h == 0L ~ 1
      )
    )

  weighted_list[[imp]] <- imp_data
}

imputed_weighted <- bind_rows(weighted_list)


# ============================================================================
# 4. ITERATIVE COVARIATE BALANCE ASSESSMENT
# ============================================================================
# Following Hayes-Larson et al. (2022): fit model, compute weights, assess
# balance, then expand model if balance is insufficient.
# Assessed on imputation 1 only for diagnostics.

imp1 <- imputed_weighted %>% filter(imp_h == 1)

# Helper: compute stabilized IOSW from a propensity score column
make_sw <- function(data, pcol) {
  p <- data[[pcol]]
  p_adrc <- sum(data$brfss_sampwt_h[data$adrc_h == 1L]) /
            sum(data$brfss_sampwt_h)
  w <- (1 - p) / p * (p_adrc / (1 - p_adrc))
  ifelse(data$adrc_h == 1L, w, 1)
}

imp1$sw1_final <- make_sw(imp1, "p1")
imp1$sw2_final <- make_sw(imp1, "p2")

# Unweighted balance
covbal_unweighted <- Checkcovbal(
  data = imp1, weight_var = NULL, covariates = covariates_m3
)

# Balance after Model 1 (demographics only)
covbal_m1 <- Checkcovbal(
  data = imp1, weight_var = "sw1_final", covariates = covariates_m1
)

# Balance after Model 2 (+ marital, comorbidities)
covbal_m2 <- Checkcovbal(
  data = imp1, weight_var = "sw2_final", covariates = covariates_m2
)

# Balance after Model 3 (FINAL: + depression)
covbal_weighted <- Checkcovbal(
  data = imp1, weight_var = "sw3_final", covariates = covariates_m3
)

# Summarize max |SMD| at each stage
summarize_balance <- function(cb, label) {
  max_smd <- cb %>% summarise(max_abs_smd = max(abs(smd), na.rm = TRUE))
  n_over_25 <- sum(abs(cb$smd) > 0.25, na.rm = TRUE)
  message(sprintf("  %-20s  max|SMD| = %.3f  covariates with |SMD|>0.25: %d",
                  label, max_smd$max_abs_smd, n_over_25))
}

message("\nIterative balance assessment (imputation 1, overall):")
summarize_balance(covbal_unweighted, "Unweighted")
summarize_balance(covbal_m1, "Model 1 (demog)")
summarize_balance(covbal_m2, "Model 2 (+comorbid)")
summarize_balance(covbal_weighted, "Model 3 (FINAL)")


# ============================================================================
# 5. SAVE
# ============================================================================

saveRDS(imputed_weighted, file.path(processed_data_dir, "05_imputed_with_weights.rds"))
saveRDS(covbal_unweighted, file.path(results_dir, "05_covbal_unweighted.rds"))
saveRDS(covbal_m1,         file.path(results_dir, "05_covbal_m1.rds"))
saveRDS(covbal_m2,         file.path(results_dir, "05_covbal_m2.rds"))
saveRDS(covbal_weighted,   file.path(results_dir, "05_covbal_weighted.rds"))

wt_summary <- imputed_weighted %>%
  filter(adrc_h == 1, imp_h == 1) %>%
  summarise(
    n        = n(),
    min_sw3  = min(sw3, na.rm = TRUE),
    p25_sw3  = quantile(sw3, 0.25, na.rm = TRUE),
    med_sw3  = median(sw3, na.rm = TRUE),
    mean_sw3 = mean(sw3, na.rm = TRUE),
    p75_sw3  = quantile(sw3, 0.75, na.rm = TRUE),
    max_sw3  = max(sw3, na.rm = TRUE)
  )

message("Weight distribution (sw3, imputation 1, ADRC only):")
message(paste(capture.output(wt_summary), collapse = "\n"))
message("\nSaved: 05_imputed_with_weights.rds")
