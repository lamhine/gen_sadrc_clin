# 07_aipsw_analysis.R
# Purpose: Primary doubly robust analysis using Augmented Inverse Probability
#          of Selection Weighting (AIPSW) to generalize ADRC outcomes to the
#          3-county BRFSS target population.
#
# Compares three estimators:
#   1. AIPSW (doubly robust) -- PRIMARY
#   2. IOSW-only (weighting only, as in KHANDLE analyses)
#   3. Outcome-model-only (g-computation)
#
# Bootstrap CIs stratified by study membership (1,000 reps x 40 imputations).
#
# INPUTS:
#   - processed_data_dir/06_outcome_predictions.rds
#
# OUTPUTS:
#   - results_dir/07_aipsw_results.rds
#   - results_dir/07_aipsw_bootests.rds
#
# ============================================================================

source("config.R")
library(boot)
library(survey)


# ============================================================================
# 1. LOAD DATA
# ============================================================================

outcome_data <- readRDS(file.path(processed_data_dir, "06_outcome_predictions.rds"))


# ============================================================================
# 2. DEFINE ESTIMATOR FUNCTIONS
# ============================================================================

# Participation model formula (same as script 05, Model 3)
# Note: MI and angina excluded from participation model (decision D06)
participation_formula <- adrc_h ~ race_eth4_h + male_h + age_h + education_h +
  marital_h + diabetes_h + hypertension_h + stroke_h +
  depression_h +
  race_eth4_h:male_h + race_eth4_h:age_h +
  race_eth4_h:education_h +
  race_eth4_h:diabetes_h + race_eth4_h:hypertension_h +
  race_eth4_h:depression_h

# Outcome model formula (same as script 06)
outcome_formula_base <- ~ race_eth4_h + male_h + age_h + education_h +
  marital_h + diabetes_h + hypertension_h + stroke_h +
  depression_h + mi_h + angina_h +
  race_eth4_h:age_h + race_eth4_h:male_h + race_eth4_h:education_h


#' Compute all three estimators for a given outcome
compute_estimates <- function(data, outcome, family = "binomial",
                              two_part = FALSE) {

  adrc  <- data %>% filter(adrc_h == 1L)
  brfss <- data %>% filter(adrc_h == 0L)
  Y     <- adrc[[outcome]]

  # ---- Participation model (refit on bootstrap sample) ----
  # Use glm() with normalized BRFSS weights inside the bootstrap for speed.
  # Script 05 uses svyglm() for the primary point estimates. Inside bootstrap
  # resamples, the bootstrap itself handles variance estimation, so normalizing
  # BRFSS weights to sum to n_brfss (preserving relative weighting but keeping
  # the ADRC:BRFSS ratio at sample-level) produces equivalent point estimates
  # without the svydesign overhead (~5x faster).
  n_brfss <- sum(data$adrc_h == 0L)
  brfss_wt_sum <- sum(data$brfss_sampwt_h[data$adrc_h == 0L])
  data$boot_wt <- ifelse(
    data$adrc_h == 1L,
    1,
    data$brfss_sampwt_h / brfss_wt_sum * n_brfss
  )

  p_mod <- glm(participation_formula, family = binomial,
               data = data, weights = boot_wt)
  p_hat <- predict(p_mod, newdata = data, type = "response")

  # Stabilization: use sample-level marginal P(ADRC) matching the normalized weights
  p_adrc_sample <- sum(data$adrc_h == 1L) / nrow(data)
  w <- (1 - p_hat) / p_hat
  sw <- w * (p_adrc_sample / (1 - p_adrc_sample))
  adrc_sw <- sw[data$adrc_h == 1L]

  # ---- Outcome model (refit on bootstrap sample) ----
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
      # Cap at valid CDR-SB range to prevent exp() explosion during extrapolation
      pred_p2 <- pmin(predict(p2_mod, newdata = data, type = "response"), 18)
    } else {
      pred_p2 <- mean(adrc_pos$cdrsb_h, na.rm = TRUE)
    }
    mu_hat <- pred_p1 * pred_p2
  }

  mu_adrc  <- mu_hat[data$adrc_h == 1L]
  mu_brfss <- mu_hat[data$adrc_h == 0L]

  # ---- Estimator 1: AIPSW (doubly robust) ----
  # Following Dahabreh et al. (2019) for generalizing from a trial/study
  # to a target population represented by survey data.
  #
  # The AIPSW estimator:
  #   θ_DR = g-comp in target + correction from source
  #
  # g-comp component: survey-weighted mean of μ(X) in BRFSS
  # Correction: IOSW-weighted residuals (Y - μ(X)) in ADRC,
  #             normalized to the ADRC sample (not the target population).
  #
  # The correction term uses IOSW normalized so that they sum to n_ADRC
  # (i.e., the correction is a weighted average of residuals).

  brfss_wt <- brfss$brfss_sampwt_h
  N_target <- sum(brfss_wt)
  n_adrc   <- length(Y)

  # g-computation: survey-weighted mean of predicted outcome in target
  est_gcomp <- sum(mu_brfss * brfss_wt) / N_target

  # IOSW-weighted mean of residuals in ADRC
  # Normalize IOSW so they sum to n_adrc (weighted average of residuals)
  iosw_normalized <- adrc_sw / sum(adrc_sw) * n_adrc
  residuals_adrc  <- Y - mu_adrc
  correction      <- weighted.mean(residuals_adrc, adrc_sw)

  est_aipsw <- est_gcomp + correction

  # ---- Estimator 2: IOSW only ----
  est_iosw <- weighted.mean(Y, adrc_sw)

  # ---- Estimator 3: Outcome model only (g-computation) ----
  # Already computed above as est_gcomp

  # ---- Unweighted ADRC ----
  est_unwtd <- mean(Y, na.rm = TRUE)

  # ---- Race/ethnicity stratified (AIPSW) ----
  races <- c("White NH", "Hispanic", "Asian NH", "NH Other")
  race_ests <- setNames(numeric(length(races)),
                        paste0("aipsw_", gsub(" ", "_", tolower(races))))

  # Race indexing within the ADRC and BRFSS subsets
  adrc_race  <- data$race_eth4_h[data$adrc_h == 1L]
  brfss_race <- data$race_eth4_h[data$adrc_h == 0L]

  for (r in races) {
    r_in_adrc  <- which(adrc_race == r)
    r_in_brfss <- which(brfss_race == r)

    if (length(r_in_adrc) < 5 || length(r_in_brfss) < 5) {
      race_ests[paste0("aipsw_", gsub(" ", "_", tolower(r)))] <- NA_real_
      next
    }

    # g-comp in target for this race
    wt_r    <- brfss_wt[r_in_brfss]
    gcomp_r <- sum(mu_brfss[r_in_brfss] * wt_r) / sum(wt_r)

    # Correction from ADRC for this race
    corr_r <- weighted.mean(
      residuals_adrc[r_in_adrc],
      adrc_sw[r_in_adrc]
    )

    race_ests[paste0("aipsw_", gsub(" ", "_", tolower(r)))] <- gcomp_r + corr_r
  }

  c(aipsw = est_aipsw, iosw = est_iosw, gcomp = est_gcomp,
    unweighted = est_unwtd, race_ests)
}


# ============================================================================
# 3. BOOTSTRAP INFERENCE
# ============================================================================

outcomes <- list(
  list(name = "mci_h",   family = "binomial", two_part = FALSE),
  list(name = "ad_h",    family = "binomial", two_part = FALSE),
  list(name = "adrd_h",  family = "binomial", two_part = FALSE),
  list(name = "moca_h",  family = "gaussian", two_part = FALSE),
  list(name = "cdrsb_h", family = "gaussian", two_part = TRUE)
)

N_BOOT <- 1000  # final analysis (D11)
imputations <- sort(unique(outcome_data$imp_h))

all_bootests <- list()
all_results  <- list()

for (out in outcomes) {
  message("Processing outcome: ", out$name)
  boot_matrix <- list()

  for (imp in imputations) {
    message("  Imputation ", imp, " of ", length(imputations))

    imp_data <- outcome_data %>%
      filter(imp_h == imp) %>%
      filter(!is.na(.data[[out$name]]) | adrc_h == 0L)

    boot_fn <- function(data, indices) {
      d <- data[indices, ]
      tryCatch(
        compute_estimates(d, outcome = out$name,
                          family = out$family,
                          two_part = out$two_part),
        error = function(e) rep(NA_real_, 8)
      )
    }

    boot_result <- boot(
      data      = imp_data,
      statistic = boot_fn,
      R         = N_BOOT,
      strata    = imp_data$adrc_h
    )

    boot_matrix[[imp]] <- boot_result$t
  }

  # ---- Combine across imputations using Rubin's rules (D19) ----
  # For each imputation m:
  #   theta_m = mean of bootstrap estimates (point estimate)
  #   V_m     = variance of bootstrap estimates (within-imputation variance)
  # Then:
  #   theta_bar = mean(theta_m) across M imputations
  #   W_bar     = mean(V_m)  -- average within-imputation variance
  #   B         = var(theta_m) -- between-imputation variance
  #   T         = W_bar + (1 + 1/M) * B  -- total variance
  #   CI: theta_bar +/- t_nu * sqrt(T)
  #   nu = (M-1) * (1 + W_bar / ((1 + 1/M) * B))^2

  est_names <- names(compute_estimates(
    outcome_data %>% filter(imp_h == 1),
    out$name, out$family, out$two_part
  ))

  M <- length(imputations)
  n_est <- length(est_names)

  # Within each imputation: point estimate and variance
  theta_m <- matrix(NA, nrow = M, ncol = n_est)
  V_m     <- matrix(NA, nrow = M, ncol = n_est)

  for (i in seq_len(M)) {
    bmat <- boot_matrix[[imputations[i]]]
    colnames(bmat) <- est_names
    theta_m[i, ] <- colMeans(bmat, na.rm = TRUE)
    V_m[i, ]     <- apply(bmat, 2, var, na.rm = TRUE)
  }

  # Rubin's combining rules
  theta_bar <- colMeans(theta_m, na.rm = TRUE)
  W_bar     <- colMeans(V_m, na.rm = TRUE)
  B         <- apply(theta_m, 2, var, na.rm = TRUE)
  Total_var <- W_bar + (1 + 1/M) * B

  # Degrees of freedom (Barnard & Rubin, 1999)
  # Avoid division by zero when B is very small
  nu <- ifelse(B > 1e-12,
               (M - 1) * (1 + W_bar / ((1 + 1/M) * B))^2,
               Inf)
  # Use normal approximation when df > 500
  t_crit <- ifelse(nu > 500, 1.96, qt(0.975, df = pmax(nu, 1)))

  ci_lower <- theta_bar - t_crit * sqrt(Total_var)
  ci_upper <- theta_bar + t_crit * sqrt(Total_var)

  names(theta_bar) <- est_names
  names(ci_lower)  <- est_names
  names(ci_upper)  <- est_names

  all_results[[out$name]] <- tibble(
    outcome  = out$name,
    estimand = est_names,
    estimate = theta_bar,
    lci      = ci_lower,
    uci      = ci_upper
  )

  # Also store the full bootstrap matrix (stacked) for prevalence ratios
  all_boots <- do.call(rbind, boot_matrix)
  colnames(all_boots) <- est_names
  all_bootests[[out$name]] <- all_boots
}

results_final <- bind_rows(all_results)


# ============================================================================
# 4. SAVE
# ============================================================================

saveRDS(results_final,  file.path(results_dir, "07_aipsw_results.rds"))
saveRDS(all_bootests,   file.path(results_dir, "07_aipsw_bootests.rds"))

message("\nSaved: 07_aipsw_results.rds")
message("Saved: 07_aipsw_bootests.rds")
message("\nResults summary:")
message(paste(capture.output(results_final), collapse = "\n"))
