# 06_outcome_models.R
# Purpose: Fit outcome models in ADRC sample for each of the four outcomes.
#          These are the OUTCOME MODEL component of the doubly robust estimator.
#
# Five outcomes with different modeling approaches:
#   1. MCI diagnosis  (binary)        -> logistic regression
#   2. AD diagnosis   (binary)        -> logistic regression
#   3. ADRD diagnosis (binary)        -> logistic regression
#   4. MoCA score     (continuous)     -> linear regression
#   5. CDR-SB score   (zero-inflated) -> two-part model
#      Part 1: P(CDR-SB > 0)         -> logistic regression
#      Part 2: E[CDR-SB | CDR-SB > 0] -> Gamma GLM (log link)
#
# DECISION D18: CDR-SB Part 2 uses Gamma GLM instead of linear regression.
# CDR-SB positive values are right-skewed (skewness = 1.39, mean = 4.6,
# median = 2.5). Gamma GLM with log link had substantially better fit
# (AIC = 1703 vs. 2045 for linear; BIC = 1829 vs. 2171).
#
# Each model uses the same harmonized covariates as the participation model
# so predictions can be generated for the full combined sample.
#
# INPUTS:
#   - processed_data_dir/05_imputed_with_weights.rds
#
# OUTPUTS:
#   - processed_data_dir/06_outcome_predictions.rds
#
# ============================================================================

source("config.R")


# ============================================================================
# 1. LOAD DATA
# ============================================================================

imputed_weighted <- readRDS(file.path(processed_data_dir, "05_imputed_with_weights.rds"))


# ============================================================================
# 2. DEFINE OUTCOME MODEL COVARIATES
# ============================================================================
# Same covariates as participation model (script 05, Model 3)

outcome_formula_base <- ~ race_eth4_h + male_h + age_h + education_h +
  marital_h + diabetes_h + hypertension_h + stroke_h +
  depression_h + mi_h + angina_h +
  race_eth4_h:age_h + race_eth4_h:male_h + race_eth4_h:education_h


# ============================================================================
# 3. FIT OUTCOME MODELS AND GENERATE PREDICTIONS
# ============================================================================

imputations <- sort(unique(imputed_weighted$imp_h))
pred_list <- list()

for (imp in imputations) {

  imp_data <- imputed_weighted %>% filter(imp_h == imp)
  adrc     <- imp_data %>% filter(adrc_h == 1L)

  # ---- 3a. MCI (binary -> logistic) ----
  fit_mci <- glm(update(outcome_formula_base, mci_h ~ .),
                 family = binomial, data = adrc)
  imp_data$pred_mci_h <- predict(fit_mci, newdata = imp_data, type = "response")

  # ---- 3b. AD (binary -> logistic) ----
  fit_ad <- glm(update(outcome_formula_base, ad_h ~ .),
                family = binomial, data = adrc)
  imp_data$pred_ad_h <- predict(fit_ad, newdata = imp_data, type = "response")

  # ---- 3c. ADRD (binary -> logistic) ----
  fit_adrd <- glm(update(outcome_formula_base, adrd_h ~ .),
                  family = binomial, data = adrc)
  imp_data$pred_adrd_h <- predict(fit_adrd, newdata = imp_data, type = "response")

  # ---- 3d. MoCA (continuous -> linear) ----
  fit_moca <- lm(update(outcome_formula_base, moca_h ~ .), data = adrc)
  imp_data$pred_moca_h <- predict(fit_moca, newdata = imp_data)

  # ---- 3e. CDR-SB (zero-inflated -> two-part model) ----
  adrc$cdrsb_pos <- as.integer(adrc$cdrsb_h > 0)

  # Part 1: P(CDR-SB > 0) -- logistic
  fit_cdrsb_p1 <- glm(update(outcome_formula_base, cdrsb_pos ~ .),
                       family = binomial, data = adrc)
  pred_p1 <- predict(fit_cdrsb_p1, newdata = imp_data, type = "response")

  # Part 2: E[CDR-SB | CDR-SB > 0] -- Gamma GLM with log link (D18)
  # Gamma is preferred over linear for right-skewed positive outcomes
  # (AIC 1703 vs 2045). Predictions are on the response scale (always > 0).
  adrc_pos <- adrc %>% filter(cdrsb_h > 0)
  if (nrow(adrc_pos) >= 10) {
    fit_cdrsb_p2 <- glm(update(outcome_formula_base, cdrsb_h ~ .),
                         family = Gamma(link = "log"), data = adrc_pos)
    # Cap predictions at valid CDR-SB range (0-18) to prevent exp() explosion
    # when extrapolating to covariate combinations not in the ADRC training data
    pred_p2 <- pmin(predict(fit_cdrsb_p2, newdata = imp_data, type = "response"), 18)
  } else {
    pred_p2 <- mean(adrc_pos$cdrsb_h, na.rm = TRUE)
    warning("Imputation ", imp, ": <10 CDR-SB > 0 cases; using unconditional mean")
  }

  # Combined: E[CDR-SB] = P(>0) * E[CDR-SB | >0]
  imp_data$pred_cdrsb_h <- pred_p1 * pred_p2

  pred_list[[imp]] <- imp_data
}

outcome_data <- bind_rows(pred_list)


# ============================================================================
# 4. SAVE
# ============================================================================

saveRDS(outcome_data, file.path(processed_data_dir, "06_outcome_predictions.rds"))

message("Saved: 06_outcome_predictions.rds")
message("Prediction columns added: pred_mci_h, pred_ad_h, pred_adrd_h, pred_moca_h, pred_cdrsb_h")
