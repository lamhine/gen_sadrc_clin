# 04_multiple_imputation.R
# Purpose: Perform multiple imputation (MICE) on harmonized ADRC + BRFSS data
#          to handle missing covariates and outcomes.
#
# Follows: Hayes-Larson et al. (2022) KHANDLE-weighting approach
#   - Predictive mean matching (PMM) for continuous/mixed variables
#   - m = 40 imputations, maxit = 10
#   - ADRC and BRFSS imputed separately, then combined
#
# IMPORTANT: Hypertension has structural missingness in BRFSS even years
# (2020, 2022) — the question was not asked. The imputation model must include
# survey_year as a predictor and allow HTN to be imputed from correlated
# covariates in the odd-year respondents.
#
# INPUTS:
#   - processed_data_dir/03_adrc_brfss_harmonized.rds
#
# OUTPUTS:
#   - processed_data_dir/04_imputed_stacked.rds  (all 40 imputations stacked)
#
# ============================================================================

source("config.R")
library(mice)

set.seed(12345)
M     <- 40     # number of imputations
MAXIT <- 10     # max iterations per imputation


# ============================================================================
# 1. LOAD HARMONIZED DATA
# ============================================================================

harmonized <- readRDS(file.path(processed_data_dir, "03_adrc_brfss_harmonized.rds"))

adrc_data  <- harmonized %>% filter(adrc_h == 1L)
brfss_data <- harmonized %>% filter(adrc_h == 0L)


# ============================================================================
# 2. DEFINE IMPUTATION VARIABLES
# ============================================================================

# Harmonized covariates present in both ADRC and BRFSS
shared_covariates <- c(
  "age_h", "male_h", "race_eth_h", "education_h", "marital_h",
  "diabetes_h", "depression_h", "mi_h", "angina_h",
  "stroke_h", "hypertension_h"
)

# Variables to impute in ADRC (covariates + outcomes)
adrc_impute_vars <- c(
  shared_covariates,
  "mci_h", "ad_h", "adrd_h", "moca_h", "cdrsb_h"
)

# Variables to impute in BRFSS (covariates only)
# Include survey_year as auxiliary predictor (not imputed) to help
# predict hypertension for even-year respondents
brfss_impute_vars <- shared_covariates


# ============================================================================
# 3. ASSESS MISSINGNESS PATTERNS
# ============================================================================

adrc_miss  <- colMeans(is.na(adrc_data[, adrc_impute_vars]))
brfss_miss <- colMeans(is.na(brfss_data[, brfss_impute_vars]))

message("ADRC missingness rates:")
message(paste(capture.output(round(adrc_miss[adrc_miss > 0], 3)), collapse = "\n"))
message("\nBRFSS missingness rates:")
message(paste(capture.output(round(brfss_miss[brfss_miss > 0], 3)), collapse = "\n"))

# Report HTN structural missingness separately
n_htn_struct <- sum(brfss_data$htn_structural_missing, na.rm = TRUE)
n_htn_total  <- sum(is.na(brfss_data$hypertension_h))
message("\nHTN: ", n_htn_struct, " of ", n_htn_total,
        " missing values are structural (even-year respondents)")


# ============================================================================
# 4. IMPUTE ADRC DATA
# ============================================================================

adrc_for_mice <- adrc_data %>%
  select(all_of(adrc_impute_vars))

# Set up MICE methods:
#   - PMM for continuous (age_h, moca_h, cdrsb_h)
#   - logreg for binary (male_h, diabetes_h, depression_h, mi_h, angina_h,
#     stroke_h, hypertension_h, mci_h, ad_h)
#   - polyreg for unordered factors (race_eth_h, education_h, marital_h)
# mice() auto-detects method based on variable type when method = ""

adrc_mice <- mice(adrc_for_mice, m = M, maxit = MAXIT,
                  printFlag = FALSE)


# ============================================================================
# 5. IMPUTE BRFSS DATA
# ============================================================================

# Include survey_year as auxiliary predictor for HTN imputation
brfss_for_mice <- brfss_data %>%
  select(all_of(brfss_impute_vars), survey_year)

# Customize predictor matrix: survey_year predicts others but is not imputed
ini <- mice(brfss_for_mice, maxit = 0)
pred <- ini$predictorMatrix
meth <- ini$method

# survey_year should not be imputed (it has no missing values)
meth["survey_year"] <- ""

brfss_mice <- mice(brfss_for_mice, m = M, maxit = MAXIT,
                   method = meth, predictorMatrix = pred,
                   printFlag = FALSE)


# ============================================================================
# 6. EXTRACT AND STACK COMPLETED DATASETS
# ============================================================================

imputed_list <- list()

for (i in 1:M) {
  adrc_imp_i <- complete(adrc_mice, i) %>%
    mutate(
      id_h           = adrc_data$id_h,
      adrc_h         = 1L,
      brfss_sampwt_h = adrc_data$brfss_sampwt_h,
      imp_h          = i,
      survey_year    = NA_integer_,
      htn_structural_missing = FALSE
    )

  brfss_imp_i <- complete(brfss_mice, i) %>%
    mutate(
      id_h           = brfss_data$id_h,
      adrc_h         = 0L,
      brfss_sampwt_h = brfss_data$brfss_sampwt_h,
      imp_h          = i,
      htn_structural_missing = brfss_data$htn_structural_missing,
      # Outcomes remain NA for BRFSS
      mci_h   = NA_integer_,
      ad_h    = NA_integer_,
      adrd_h  = NA_integer_,
      moca_h  = NA_real_,
      cdrsb_h = NA_real_
    )

  imputed_list[[i]] <- bind_rows(adrc_imp_i, brfss_imp_i)
}

imputed_stacked <- bind_rows(imputed_list)

# Derive 4-category race variable for Aims 2-3 (deterministic from imputed race_eth_h)
# Collapse Black NH + Other/Multi NH -> NH Other
imputed_stacked <- imputed_stacked %>%
  mutate(
    race_eth4_h = fct_collapse(race_eth_h,
      "NH Other" = c("Black NH", "Other/Multi NH")
    ) %>%
      fct_relevel("NH Other", after = Inf)
  )


# ============================================================================
# 7. SAVE
# ============================================================================

saveRDS(imputed_stacked, file.path(processed_data_dir, "04_imputed_stacked.rds"))

message("Saved: 04_imputed_stacked.rds (",
        nrow(imputed_stacked), " rows = ",
        M, " imputations x ",
        nrow(imputed_stacked) / M, " observations)")
