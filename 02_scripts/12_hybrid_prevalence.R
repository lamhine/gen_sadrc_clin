# 12_hybrid_prevalence.R
# Purpose: Estimate local AD prevalence using a hybrid model that combines
#          published population-based regression coefficients from Dhana et al.
#          (2023, CHAP) with ADRC-derived race/ethnicity coefficients, applied
#          to individual-level BRFSS microdata.
#
# APPROACH (Decision D24):
#   - Dhana provides: intercept, age, sex, education, Black coefficients
#     (from CHAP, a population-based study)
#   - ADRC provides: Hispanic and Asian NH race coefficients (local Bay Area,
#     fills Dhana's Asian gap, updates Hispanic for local population)
#   - NH Other: uses Dhana's Black coefficient as proxy (ADRC N=18 too small)
#   - Applied to BRFSS individual-level data with survey weights
#
# Also computes:
#   - Dhana-only estimates (for comparison / replication)
#   - Matthews indirect standardization for ADRD (separate from hybrid AD model)
#   - Enrichment diagnostic: ADRC ORs vs Dhana published coefficients
#
# INPUTS:
#   - processed_data_dir/05_imputed_with_weights.rds  (for ADRC race ORs)
#   - processed_data_dir/03_adrc_brfss_harmonized.rds (for BRFSS application)
#
# OUTPUTS:
#   - results_dir/12_hybrid_prevalence.rds
#   - results_dir/12_or_comparison.rds
#
# Reference: Dhana KA et al. (2023). Prevalence of Alzheimer's disease
#   dementia in the 50 US states and 3142 counties. Alzheimers Dement.
#   doi: 10.1002/alz.13081
#
# ============================================================================

source("config.R")
library(survey)


# ============================================================================
# 1. DHANA ET AL. (2023) PUBLISHED COEFFICIENTS
# ============================================================================
# Table 1 from Dhana et al. 2023
# Outcome: AD prevalence (clinical diagnosis from CHAP)
# Model: generalized additive quasibinomial regression
# Population: CHAP participants aged 65+ (Chicago, IL)

dhana_coefs <- list(
  intercept = -3.455,
  # Age groups (ref: 65-69)
  age_70_74 = 0.577,
  age_75_79 = 1.126,
  age_80_84 = 1.800,
  age_85plus = 2.693,
  # Sex (ref: male)
  female = 0.123,
  # Race/ethnicity (ref: White)
  black = 0.915,
  hispanic = 0.548,
  # Note: Asian NOT in Dhana model (assumed = White, i.e., beta = 0)
  # Education (per SD increase in years of schooling)
  # Dhana reports per SD; CHAP education SD is approximately 3.5 years
  education_per_sd = -0.398
)

# CHAP education SD (approximate, from CHAP descriptive statistics)
chap_educ_sd <- 3.5


# ============================================================================
# 2. MATTHEWS ET AL. (2019) ADRD PREVALENCE RATES
# ============================================================================
# For indirect standardization of ADRD (separate from AD hybrid model)

matthews_race <- tibble(
  race_eth4_h = c("White NH", "Hispanic", "Asian NH", "NH Other"),
  adrd_prev_matthews = c(0.103, 0.122, 0.084, 0.138)
)

matthews_age <- tibble(
  age_group = c("65-74", "75-84", "85+"),
  adrd_prev_matthews = c(0.036, 0.136, 0.346)
)

matthews_sex <- tibble(
  sex = c("Female", "Male"),
  adrd_prev_matthews = c(0.122, 0.086)
)


# ============================================================================
# 3. LOAD DATA
# ============================================================================

imputed_weighted <- readRDS(file.path(processed_data_dir,
                                       "05_imputed_with_weights.rds"))
harmonized <- readRDS(file.path(processed_data_dir,
                                 "03_adrc_brfss_harmonized.rds"))


# ============================================================================
# 4. ESTIMATE ADRC RACE/ETHNICITY COEFFICIENTS
# ============================================================================
# Fit AD model in ADRC to extract race ORs (case-control principle: ORs
# preserved under approximately constant enrichment across race groups).
# We match Dhana's model specification: age group + sex + race + education.

imputations <- sort(unique(imputed_weighted$imp_h))

adrc_or_list <- list()

for (imp in imputations) {
  imp_data <- imputed_weighted %>% filter(imp_h == imp)
  adrc <- imp_data %>%
    filter(adrc_h == 1L) %>%
    mutate(
      age_group = cut(age_h, breaks = c(65, 70, 75, 80, 85, Inf),
                      labels = c("65-69", "70-74", "75-79", "80-84", "85+"),
                      right = FALSE),
      female = 1L - male_h,
      # Education in years (approximate from categories for comparison)
      # education_h is a factor: 1=<HS, 2=HS, 3=Some college, 4=College+
      educ_years_approx = case_when(
        education_h == "Less than HS"     ~ 9,
        education_h == "HS diploma/GED"   ~ 12,
        education_h == "Some college"     ~ 14,
        education_h == "College graduate+" ~ 17,
        TRUE ~ NA_real_
      )
    )

  # Fit AD model matching Dhana's specification
  # Note: using ad_h (AD specifically) to match Dhana's outcome
  fit_ad <- glm(ad_h ~ age_group + female + race_eth4_h + educ_years_approx,
                family = binomial, data = adrc)

  coefs <- coef(fit_ad)
  se <- sqrt(diag(vcov(fit_ad)))

  adrc_or_list[[imp]] <- tibble(
    imp_h = imp,
    term = names(coefs),
    beta = coefs,
    se = se,
    or = exp(coefs)
  )
}

adrc_ors <- bind_rows(adrc_or_list)

# Pool across imputations (Rubin's rules for regression coefficients)
adrc_ors_pooled <- adrc_ors %>%
  group_by(term) %>%
  summarise(
    beta_bar = mean(beta),
    W_bar = mean(se^2),       # within-imputation variance
    B = var(beta),             # between-imputation variance
    M = n(),
    T_var = W_bar + (1 + 1/M) * B,
    se_pooled = sqrt(T_var),
    or_pooled = exp(beta_bar),
    or_lci = exp(beta_bar - 1.96 * se_pooled),
    or_uci = exp(beta_bar + 1.96 * se_pooled),
    .groups = "drop"
  )

message("=== ADRC race/ethnicity ORs (pooled across imputations) ===")
message(paste(capture.output(
  print(adrc_ors_pooled %>%
          select(term, beta_bar, se_pooled, or_pooled, or_lci, or_uci) %>%
          mutate(across(where(is.numeric), ~ round(., 3))),
        n = Inf)
), collapse = "\n"))


# ============================================================================
# 5. COMPARE ADRC ORs TO DHANA COEFFICIENTS (Enrichment Diagnostic)
# ============================================================================

# Extract comparable coefficients
or_comparison <- tibble(
  variable = c("Intercept", "Age 70-74", "Age 75-79", "Age 80-84", "Age 85+",
               "Female", "Hispanic", "Education (per year)"),
  dhana_beta = c(
    dhana_coefs$intercept,
    dhana_coefs$age_70_74, dhana_coefs$age_75_79,
    dhana_coefs$age_80_84, dhana_coefs$age_85plus,
    dhana_coefs$female,
    dhana_coefs$hispanic,
    dhana_coefs$education_per_sd / chap_educ_sd  # convert per-SD to per-year
  ),
  dhana_or = exp(c(
    dhana_coefs$intercept,
    dhana_coefs$age_70_74, dhana_coefs$age_75_79,
    dhana_coefs$age_80_84, dhana_coefs$age_85plus,
    dhana_coefs$female,
    dhana_coefs$hispanic,
    dhana_coefs$education_per_sd / chap_educ_sd
  ))
)

# Map ADRC terms to comparison names
adrc_map <- c(
  "(Intercept)" = "Intercept",
  "age_group70-74" = "Age 70-74",
  "age_group75-79" = "Age 75-79",
  "age_group80-84" = "Age 80-84",
  "age_group85+" = "Age 85+",
  "female" = "Female",
  "race_eth4_hHispanic" = "Hispanic",
  "educ_years_approx" = "Education (per year)"
)

adrc_compare <- adrc_ors_pooled %>%
  filter(term %in% names(adrc_map)) %>%
  mutate(variable = adrc_map[term]) %>%
  select(variable, adrc_beta = beta_bar, adrc_or = or_pooled,
         adrc_or_lci = or_lci, adrc_or_uci = or_uci)

or_comparison <- or_comparison %>%
  left_join(adrc_compare, by = "variable") %>%
  mutate(
    or_ratio = adrc_or / dhana_or,
    note = case_when(
      variable == "Intercept" ~ "Expected to differ (enrichment shifts intercept)",
      abs(log(or_ratio)) < 0.3 ~ "Similar (ratio within 0.74-1.35)",
      TRUE ~ "Divergent — enrichment may distort this association"
    )
  )

message("\n=== OR Comparison: ADRC vs Dhana (population-based) ===")
message(paste(capture.output(
  print(or_comparison %>%
          mutate(across(where(is.numeric), ~ round(., 3))) %>%
          select(variable, dhana_or, adrc_or, adrc_or_lci, adrc_or_uci, or_ratio, note),
        n = Inf, width = 120)
), collapse = "\n"))

# Also show Asian NH OR (not in Dhana — this is the ADRC's unique contribution)
asian_or <- adrc_ors_pooled %>%
  filter(term == "race_eth4_hAsian NH")
message("\nADRC Asian NH OR (vs White NH): ",
        round(asian_or$or_pooled, 3),
        " (", round(asian_or$or_lci, 3), ", ", round(asian_or$or_uci, 3), ")")
message("Dhana assumes Asian = White (OR = 1.0)")


# ============================================================================
# 6. BUILD HYBRID MODEL AND APPLY TO BRFSS
# ============================================================================
# Hybrid: Dhana base + ADRC race coefficients for Hispanic and Asian NH

# Extract ADRC race betas (pooled)
adrc_hispanic_beta <- adrc_ors_pooled %>%
  filter(term == "race_eth4_hHispanic") %>% pull(beta_bar)
adrc_asian_beta <- adrc_ors_pooled %>%
  filter(term == "race_eth4_hAsian NH") %>% pull(beta_bar)

message("\nHybrid model race coefficients:")
message("  Hispanic: Dhana = ", dhana_coefs$hispanic,
        " -> ADRC = ", round(adrc_hispanic_beta, 3))
message("  Asian NH: Dhana = 0 (assumed White) -> ADRC = ",
        round(adrc_asian_beta, 3))
message("  NH Other: using Dhana Black = ", dhana_coefs$black,
        " (ADRC N too small)")

# Prepare BRFSS data for prediction
brfss <- harmonized %>%
  filter(adrc_h == 0L) %>%
  mutate(
    age_group = cut(age_h, breaks = c(65, 70, 75, 80, 85, Inf),
                    labels = c("65-69", "70-74", "75-79", "80-84", "85+"),
                    right = FALSE),
    female = 1L - male_h,
    # Map education categories to approximate years
    educ_years_approx = case_when(
      education_h == "Less than HS"      ~ 9,
      education_h == "HS diploma/GED"    ~ 12,
      education_h == "Some college"      ~ 14,
      education_h == "College graduate+" ~ 17,
      TRUE ~ NA_real_
    ),
    # Standardize education relative to CHAP SD
    educ_z = (educ_years_approx - 12) / chap_educ_sd  # center at 12 years
  )

# Drop rows with missing key variables
brfss_complete <- brfss %>%
  filter(!is.na(age_group), !is.na(female), !is.na(race_eth4_h),
         !is.na(educ_years_approx))

message("\nBRFSS complete cases for hybrid model: ", nrow(brfss_complete),
        " / ", nrow(brfss))

# ---- Dhana-only model (replication) ----
brfss_complete <- brfss_complete %>%
  mutate(
    eta_dhana = dhana_coefs$intercept +
      ifelse(age_group == "70-74", dhana_coefs$age_70_74, 0) +
      ifelse(age_group == "75-79", dhana_coefs$age_75_79, 0) +
      ifelse(age_group == "80-84", dhana_coefs$age_80_84, 0) +
      ifelse(age_group == "85+",   dhana_coefs$age_85plus, 0) +
      female * dhana_coefs$female +
      ifelse(race_eth4_h == "NH Other", dhana_coefs$black, 0) +
      ifelse(race_eth4_h == "Hispanic", dhana_coefs$hispanic, 0) +
      # Asian = 0 in Dhana (assumed same as White)
      educ_z * dhana_coefs$education_per_sd,
    p_dhana = plogis(eta_dhana)
  )

# ---- Hybrid model (Dhana base + ADRC race) ----
brfss_complete <- brfss_complete %>%
  mutate(
    eta_hybrid = dhana_coefs$intercept +
      ifelse(age_group == "70-74", dhana_coefs$age_70_74, 0) +
      ifelse(age_group == "75-79", dhana_coefs$age_75_79, 0) +
      ifelse(age_group == "80-84", dhana_coefs$age_80_84, 0) +
      ifelse(age_group == "85+",   dhana_coefs$age_85plus, 0) +
      female * dhana_coefs$female +
      # ADRC race coefficients for Hispanic and Asian
      ifelse(race_eth4_h == "Hispanic", adrc_hispanic_beta, 0) +
      ifelse(race_eth4_h == "Asian NH", adrc_asian_beta, 0) +
      # Dhana Black coefficient for NH Other
      ifelse(race_eth4_h == "NH Other", dhana_coefs$black, 0) +
      educ_z * dhana_coefs$education_per_sd,
    p_hybrid = plogis(eta_hybrid)
  )

# ---- Compute survey-weighted prevalence estimates ----
# Overall
wt <- brfss_complete$brfss_sampwt_h

est_dhana_overall  <- weighted.mean(brfss_complete$p_dhana, wt)
est_hybrid_overall <- weighted.mean(brfss_complete$p_hybrid, wt)

message("\n======================================")
message("LOCAL AD PREVALENCE ESTIMATES (65+)")
message("======================================")
message(sprintf("Dhana-only (replication):  %.1f%%", est_dhana_overall * 100))
message(sprintf("Hybrid (Dhana + ADRC race): %.1f%%", est_hybrid_overall * 100))

# By race/ethnicity
race_estimates <- brfss_complete %>%
  group_by(race_eth4_h) %>%
  summarise(
    n_brfss = n(),
    dhana_prev  = weighted.mean(p_dhana, brfss_sampwt_h) * 100,
    hybrid_prev = weighted.mean(p_hybrid, brfss_sampwt_h) * 100,
    .groups = "drop"
  )

message("\nBy race/ethnicity:")
message(paste(capture.output(
  print(race_estimates %>% mutate(across(where(is.numeric) & !matches("n_"),
                                          ~ round(., 1))),
        n = Inf)
), collapse = "\n"))

# By age group
age_estimates <- brfss_complete %>%
  group_by(age_group) %>%
  summarise(
    n_brfss = n(),
    dhana_prev  = weighted.mean(p_dhana, brfss_sampwt_h) * 100,
    hybrid_prev = weighted.mean(p_hybrid, brfss_sampwt_h) * 100,
    .groups = "drop"
  )

message("\nBy age group:")
message(paste(capture.output(
  print(age_estimates %>% mutate(across(where(is.numeric) & !matches("n_"),
                                         ~ round(., 1))),
        n = Inf)
), collapse = "\n"))

# By sex
sex_estimates <- brfss_complete %>%
  mutate(sex_label = ifelse(female == 1, "Female", "Male")) %>%
  group_by(sex_label) %>%
  summarise(
    n_brfss = n(),
    dhana_prev  = weighted.mean(p_dhana, brfss_sampwt_h) * 100,
    hybrid_prev = weighted.mean(p_hybrid, brfss_sampwt_h) * 100,
    .groups = "drop"
  )

message("\nBy sex:")
message(paste(capture.output(
  print(sex_estimates %>% mutate(across(where(is.numeric) & !matches("n_"),
                                         ~ round(., 1))),
        n = Inf)
), collapse = "\n"))


# ============================================================================
# 7. MATTHEWS INDIRECT STANDARDIZATION FOR ADRD
# ============================================================================
# Apply Matthews race x age x sex rates to local population structure

brfss_adrd <- brfss_complete %>%
  mutate(
    age_group_matthews = case_when(
      age_h >= 65 & age_h < 75 ~ "65-74",
      age_h >= 75 & age_h < 85 ~ "75-84",
      age_h >= 85              ~ "85+",
      TRUE ~ NA_character_
    ),
    sex_label = ifelse(female == 1, "Female", "Male")
  ) %>%
  left_join(matthews_race, by = "race_eth4_h") %>%
  left_join(matthews_age, by = c("age_group_matthews" = "age_group")) %>%
  left_join(matthews_sex, by = c("sex_label" = "sex"))

# Simple indirect standardization using race-specific rates
adrd_race_std <- weighted.mean(brfss_adrd$adrd_prev_matthews.x,
                                brfss_adrd$brfss_sampwt_h, na.rm = TRUE)

message("\n======================================")
message("LOCAL ADRD PREVALENCE (Matthews indirect standardization)")
message("======================================")
message(sprintf("Overall (race-standardized): %.1f%%", adrd_race_std * 100))

adrd_by_race <- brfss_adrd %>%
  group_by(race_eth4_h) %>%
  summarise(
    n = n(),
    matthews_rate = first(adrd_prev_matthews.x) * 100,
    .groups = "drop"
  )
message("\nMatthews ADRD rates applied to local population:")
message(paste(capture.output(print(adrd_by_race, n = Inf)), collapse = "\n"))


# ============================================================================
# 8. SAVE
# ============================================================================

saveRDS(list(
  # Hybrid AD prevalence
  overall = tibble(
    model = c("Dhana only", "Hybrid (Dhana + ADRC race)"),
    ad_prev_pct = c(est_dhana_overall * 100, est_hybrid_overall * 100)
  ),
  by_race = race_estimates,
  by_age = age_estimates,
  by_sex = sex_estimates,
  # ADRD indirect standardization
  adrd_overall_pct = adrd_race_std * 100,
  adrd_by_race = adrd_by_race,
  # Model coefficients
  dhana_coefs = dhana_coefs,
  adrc_hispanic_beta = adrc_hispanic_beta,
  adrc_asian_beta = adrc_asian_beta
), file.path(results_dir, "12_hybrid_prevalence.rds"))

saveRDS(list(
  or_comparison = or_comparison,
  adrc_ors_pooled = adrc_ors_pooled
), file.path(results_dir, "12_or_comparison.rds"))

message("\nSaved: 12_hybrid_prevalence.rds")
message("Saved: 12_or_comparison.rds")
