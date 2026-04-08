# 03_harmonize.R
# Purpose: Harmonize variable names and coding across ADRC Clinical Core
#          (NACC UDS v3 format) and California BRFSS (2019-2023, 3-county
#          catchment) to create a combined analysis dataset.
#
# AUTHORITATIVE REFERENCE: User-provided harmonization map (March 2026).
# NACC form references: UDS v3 Initial Visit Packet (forms A1, A5, D1, C2, B4).
# BRFSS references: brfss_vars_summary.xlsx (by_year and by_var sheets).
#
# All harmonized variables use the "_h" suffix.
# Study membership indicator: adrc_h (1 = ADRC, 0 = BRFSS)
#
# INPUTS:
#   - processed_data_dir/02_adrc_clean.rds
#   - processed_data_dir/02_brfss_catchment.rds
#
# OUTPUTS:
#   - processed_data_dir/03_adrc_brfss_harmonized.rds
#   - processed_data_dir/03_adrc_zip_for_acs.rds
#
# ============================================================================
# CONFIRMED HARMONIZABLE VARIABLES (for participation model)
# ============================================================================
#
# +--------------------+---------------------------------+-----------------------------------+--------+
# | Harmonized name    | NACC source                     | BRFSS source                      | Status |
# +--------------------+---------------------------------+-----------------------------------+--------+
# | age_h              | a1_birthyr + a1_visitdate       | AGE (continuous)                  | Full   |
# | male_h             | a1_sex (A1 Q4: 1=M, 2=F)       | sex_unified (SEX2/BIRTHSEX/_SEX)  | Full   |
# | race_eth_h (5-cat) | a1_hispanic (Q7) + a1_race (Q9) | HISPANC3 + race_unified           | Full   |
# | education_h (4-cat)| a1_educ (Q11: years -> 4 cats)  | _EDUCAG (stable 2019-2023)        | Full   |
# | marital_h (6-cat)  | a1_maristat (Q12, recoded)      | MARITAL (stable 2019-2023)        | Full   |
# | diabetes_h         | a5_diabetes (A5 Q5a: 0/1/2/9)   | diabetes_unified (DIABCOR3/DIABETE4)| Full |
# | depression_h       | a5_dep2yrs (A5 Q7d: 0/1/2/9)   | depression_unified (DEPRESS1/ADDEPEV3)| Full|
# | mi_h               | a5_cvhatt (A5 Q2a: 0/1/2/9)    | mi_unified (HEART2/CVDINFR4)      | Full   |
# | angina_h           | a5_cvangina (A5 Q2g: 0/1/2/9)  | angina_unified (ANGINA/CVDCRHD4)  | Full   |
# | stroke_h           | a5_cbstroke (A5 Q3a: 0/1/2/9)   | stroke_unified (STROKE2/CVDSTRK3) | Full   |
# | hypertension_h     | a5_hyperten (A5 Q5b: 0/1/2/9)  | htn_unified (BPHIGH2/3/6)         | Full*  |
# +--------------------+---------------------------------+-----------------------------------+--------+
# * HTN: NOT asked in BRFSS even years (2020, 2022) -> structural missingness.
#   Must be handled in imputation model (script 04) with year indicator.
#
# ============================================================================
# CONFIRMED GAPS (omit from participation model, document as limitations)
# ============================================================================
#
# Income: NACC has only a1_zip (first 3 digits), no individual income.
#         BRFSS has household income (INCOM02/INCOM03/INCOME3).
#         Cannot harmonize. Omit. For Aim 1 only, may link ZIP to ACS area-level income.
#
# Self-rated health (GENHLTH): NOT available in NACC UDS v3.
#
# Exercise: NOT available in NACC UDS v3.
#
# Vision/blindness: NOT in NACC A1 or A5.
#
# Military/veteran: NOT reliably harmonizable (a1_veteran has high missingness
#                   and different coding from BRFSS). Omit from participation model.
#
# ADL walking/dressing: BRFSS has DIFFWALK and DIFDRES2/DIFFDRES.
#                       NACC has no direct equivalent. Do not proxy. Omit.
#
# ============================================================================
# OUTCOMES (ADRC only)
# ============================================================================
#   mci_h   = MCI diagnosis (binary): d1_mciamem | d1_mciaplus | d1_mcinon1 | d1_mcinon2
#   ad_h    = Alzheimer's disease etiology (binary): d1_alzdis == 1
#   adrd_h  = AD + related dementias (binary): d1_alzdis | d1_lbdis | d1_cvd | d1_ftldnos
#   moca_h  = MoCA total score (continuous, 0-30): c2_mocatots
#   cdrsb_h = CDR Sum of Boxes (continuous, 0-18, zero-inflated): b4_cdrsum
#
# ============================================================================

source("config.R")


# ============================================================================
# 1. LOAD CLEANED DATA
# ============================================================================

adrc_clean <- readRDS(file.path(processed_data_dir, "02_adrc_clean.rds"))
brfss_clean <- readRDS(file.path(processed_data_dir, "02_brfss_catchment.rds"))


# ============================================================================
# 2. HARMONIZE ADRC (NACC) VARIABLES
# ============================================================================

adrc_h <- adrc_clean %>%
  mutate(
    # ---- Identifiers & study membership ----
    id_h   = as.character(record_id),
    adrc_h = 1L,
    brfss_sampwt_h = 1,   # unit weight for ADRC participants

    # ---- Age ----
    # Computed from A1 birth year/month to visit date
    # Visit date is in m/d/yy format (e.g., "9/10/15"); use lubridate::mdy()
    visit_date_parsed = lubridate::mdy(a1_visitdate),
    birth_date_approx = as.Date(paste0(a1_birthyr, "-",
                                        sprintf("%02d", a1_birthmo), "-15")),
    age_h = as.numeric(
      difftime(visit_date_parsed, birth_date_approx, units = "days")
    ) / 365.25,

    # ---- Sex ----
    # NACC Form A1 Q4: 1=Male, 2=Female
    # Same coding direction as BRFSS. Harmonize to binary male indicator.
    male_h = case_when(
      a1_sex == 1L ~ 1L,   # Male
      a1_sex == 2L ~ 0L,   # Female
      TRUE ~ NA_integer_
    ),

    # ---- Race / Ethnicity (5 categories) ----
    # NACC Form A1: a1_hispanic (Q7: 0=No, 1=Yes Hispanic, 9=Unknown)
    #               a1_race (Q9: 1=White, 2=Black, 3=AIAN, 4=NHPI,
    #                            5=Asian, 50=Other, 99=Unknown)
    # Hispanic ethnicity takes priority over race (any race).
    # Target categories: White NH, Black NH, Hispanic, Asian NH, Other/Multiracial NH
    race_eth_h = case_when(
      a1_hispanic == 1L                                     ~ "Hispanic",
      a1_hispanic == 0L & a1_race == 1L                     ~ "White NH",
      a1_hispanic == 0L & a1_race == 2L                     ~ "Black NH",
      a1_hispanic == 0L & a1_race == 5L                     ~ "Asian NH",
      # AIAN, NHPI, Other, Multiracial, secondary/tertiary race -> Other/Multi NH
      a1_hispanic == 0L & a1_race %in% c(3L, 4L, 50L)      ~ "Other/Multi NH",
      a1_hispanic == 0L & !is.na(a1_racesec) &
        a1_racesec != 99L                                    ~ "Other/Multi NH",
      a1_hispanic == 9L                                     ~ NA_character_,
      TRUE ~ NA_character_
    ),

    # ---- Education (4 categories) ----
    # NACC Form A1 Q11: a1_educ = years of education (0-36; 99=unknown)
    # Recode to 4 categories matching BRFSS _EDUCAG
    education_h = case_when(
      a1_educ >= 0  & a1_educ <= 11 ~ 1L,   # Did not graduate high school
      a1_educ == 12                  ~ 2L,   # Graduated high school (HS/GED)
      a1_educ >= 13 & a1_educ <= 15 ~ 3L,   # Attended college or technical school
      a1_educ >= 16 & a1_educ <= 36 ~ 4L,   # Graduated from college or technical school
      a1_educ == 99                  ~ NA_integer_,
      TRUE ~ NA_integer_
    ),

    # ---- Marital status (6 categories, recoded to match BRFSS) ----
    # NACC Form A1 Q12: a1_maristat
    #   NACC coding:  1=Married, 2=Widowed, 3=Divorced, 4=Separated,
    #                 5=Never married, 6=Living as married/domestic partner, 9=Unknown
    #   BRFSS coding: 1=Married, 2=Divorced, 3=Widowed, 4=Separated,
    #                 5=Never married, 6=Unmarried couple, 9=Refused
    # CRITICAL: NACC and BRFSS swap codes for Widowed (NACC=2/BRFSS=3) and
    #           Divorced (NACC=3/BRFSS=2). Recode NACC to match BRFSS.
    marital_h = case_when(
      a1_maristat == 1L ~ 1L,   # Married -> 1 (same in both)
      a1_maristat == 2L ~ 3L,   # NACC Widowed=2 -> BRFSS Widowed=3
      a1_maristat == 3L ~ 2L,   # NACC Divorced=3 -> BRFSS Divorced=2
      a1_maristat == 4L ~ 4L,   # Separated -> 4 (same in both)
      a1_maristat == 5L ~ 5L,   # Never married -> 5 (same in both)
      a1_maristat == 6L ~ 6L,   # Living as married -> 6 (unmarried couple)
      a1_maristat == 9L ~ NA_integer_,   # Unknown -> NA
      TRUE ~ NA_integer_
    ),

    # ---- Diabetes (binary) ----
    # NACC Form A5 Q5a: a5_diabetes
    #   0=absent, 1=recent/active, 2=remote/inactive, 9=unknown
    #   Recode: 1 or 2 -> yes (ever diagnosed), 0 -> no, 9 -> NA
    diabetes_h = case_when(
      a5_diabetes %in% c(1L, 2L) ~ 1L,   # Ever diagnosed (active or remote)
      a5_diabetes == 0L           ~ 0L,   # Absent
      a5_diabetes == 9L           ~ NA_integer_,
      TRUE ~ NA_integer_
    ),

    # ---- Depression (binary) ----
    # NACC Form A5 Q7d: a5_dep2yrs (depression lasting 2+ years)
    #   0=absent, 1=recent/active, 2=remote/inactive, 9=unknown
    depression_h = case_when(
      a5_dep2yrs %in% c(1L, 2L) ~ 1L,
      a5_dep2yrs == 0L           ~ 0L,
      a5_dep2yrs == 9L           ~ NA_integer_,
      TRUE ~ NA_integer_
    ),

    # ---- Heart attack / MI (binary) ----
    # NACC Form A5 Q2a: a5_cvhatt (heart attack/cardiac arrest)
    #   0=absent, 1=recent/active, 2=remote/inactive, 9=unknown
    mi_h = case_when(
      a5_cvhatt %in% c(1L, 2L) ~ 1L,
      a5_cvhatt == 0L           ~ 0L,
      a5_cvhatt == 9L           ~ NA_integer_,
      TRUE ~ NA_integer_
    ),

    # ---- Angina / CVD (binary) ----
    # NACC Form A5 Q2g: a5_cvangina
    #   0=absent, 1=recent/active, 2=remote/inactive, 9=unknown
    angina_h = case_when(
      a5_cvangina %in% c(1L, 2L) ~ 1L,
      a5_cvangina == 0L           ~ 0L,
      a5_cvangina == 9L           ~ NA_integer_,
      TRUE ~ NA_integer_
    ),

    # ---- Stroke (binary) ----
    # NACC Form A5 Q3a: a5_cbstroke
    #   0=absent, 1=recent/active, 2=remote/inactive, 9=unknown
    stroke_h = case_when(
      a5_cbstroke %in% c(1L, 2L) ~ 1L,
      a5_cbstroke == 0L           ~ 0L,
      a5_cbstroke == 9L           ~ NA_integer_,
      TRUE ~ NA_integer_
    ),

    # ---- Hypertension (binary) ----
    # NACC Form A5 Q5b: a5_hyperten
    #   0=absent, 1=recent/active, 2=remote/inactive, 9=unknown
    hypertension_h = case_when(
      a5_hyperten %in% c(1L, 2L) ~ 1L,
      a5_hyperten == 0L           ~ 0L,
      a5_hyperten == 9L           ~ NA_integer_,
      TRUE ~ NA_integer_
    ),

    # ---- Structural missingness flag for HTN ----
    # ADRC always has HTN from A5 -> no structural missingness
    htn_structural_missing = FALSE,

    # ---- Outcomes (ADRC only) ----

    # MCI diagnosis (binary): any MCI subtype from Form D1
    # D1: d1_mciamem (amnestic single), d1_mciaplus (amnestic multi),
    #     d1_mcinon1 (non-amnestic single), d1_mcinon2 (non-amnestic multi)
    mci_h = case_when(
      as.integer(d1_mciamem) == 1L | as.integer(d1_mciaplus) == 1L |
        as.integer(d1_mcinon1) == 1L | as.integer(d1_mcinon2) == 1L ~ 1L,
      as.integer(d1_normcog) == 1L                                   ~ 0L,
      as.integer(d1_normcog) == 0L                                   ~ 0L,
      TRUE ~ NA_integer_
    ),

    # Alzheimer's disease etiology (binary): Form D1
    # d1_alzdis: 0 = No AD, 1 = AD is contributing/primary etiology
    # d1_alzdis is only coded for participants with dementia etiology assessment.
    # For cognitively normal participants (d1_normcog==1), AD is 0 by definition.
    # For impaired participants without d1_alzdis coded, treat as 0 (non-AD).
    ad_h = case_when(
      as.integer(d1_alzdis) == 1L              ~ 1L,   # AD etiology present
      as.integer(d1_alzdis) == 0L              ~ 0L,   # Assessed, no AD
      as.integer(d1_normcog) == 1L             ~ 0L,   # Normal cognition -> no AD
      !is.na(oci) & is.na(as.integer(d1_alzdis)) ~ 0L, # Impaired, etiology not AD
      TRUE ~ NA_integer_
    ),

    # ADRD — Alzheimer's disease and related dementias (binary): Form D1
    # Includes AD (d1_alzdis), Lewy body (d1_lbdis), cerebrovascular (d1_cvd),
    # and frontotemporal (d1_ftldnos). Maps to Medicare claims-based ADRD
    # definition used in Matthews et al. (2019) for population prevalence.
    adrd_h = case_when(
      as.integer(d1_alzdis) == 1L | as.integer(d1_lbdis) == 1L |
        as.integer(d1_cvd) == 1L | as.integer(d1_ftldnos) == 1L ~ 1L,
      as.integer(d1_normcog) == 1L                               ~ 0L, # Normal -> no ADRD
      !is.na(oci)                                                 ~ 0L, # Impaired, no ADRD etiology
      TRUE ~ NA_integer_
    ),

    # MoCA total score (continuous, 0-30): Form C2
    moca_h = as.numeric(c2_mocatots),

    # CDR Sum of Boxes (continuous, 0-18, zero-inflated): Form B4
    cdrsb_h = as.numeric(b4_cdrsum)
  ) %>%
  select(
    id_h, adrc_h, brfss_sampwt_h,
    age_h, male_h, race_eth_h, education_h, marital_h,
    diabetes_h, depression_h, mi_h, angina_h, stroke_h, hypertension_h,
    htn_structural_missing,
    mci_h, ad_h, adrd_h, moca_h, cdrsb_h,
    a1_zip   # keep for ACS linkage
  )


# ============================================================================
# 3. HARMONIZE BRFSS VARIABLES
# ============================================================================
# Uses unified columns created in 02_clean_data.R:
#   sex_unified, stroke_unified, htn_unified, diabetes_unified,
#   depression_unified, mi_unified, angina_unified, race_unified

brfss_h <- brfss_clean %>%
  mutate(
    # ---- Identifiers & study membership ----
    id_h   = paste0("brfss_", row_number()),
    adrc_h = 0L,

    # Survey weight: _LLCPWT is the combined landline + cell phone final weight.
    # Per CDC guidance on pooling multiple BRFSS years: divide by the number of
    # pooled years (N_BRFSS_YEARS) so that weights sum to the target population
    # once rather than N times. See:
    #   CDC. Complex Sampling Weights and Preparing Module Data for Analysis. 2022.
    # Decision log: FLAG-F / D13.
    brfss_sampwt_h = as.numeric(`_LLCPWT`) / N_BRFSS_YEARS,

    # ---- Age ----
    # BRFSS AGE: continuous (reported age in years)
    # Missing codes: 7, 9 are NOT NA for AGE per by_var sheet ("7,9=NA; 77,99=valid")
    # Wait — the by_var sheet says "7,9=NA; 77,99=valid" which means single-digit
    # 7 and 9 are NA. Ages 77 and 99 are valid ages.
    age_h = {
      a <- as.numeric(AGE)
      ifelse(a %in% c(7, 9), NA_real_, a)
    },

    # ---- Sex ----
    # sex_unified created in 02_clean_data.R: coalesce(_SEX, SEX2, BIRTHSEX)
    # All coded 1=Male, 2=Female. Same coding as NACC.
    male_h = case_when(
      sex_unified == 1L ~ 1L,   # Male
      sex_unified == 2L ~ 0L,   # Female
      TRUE ~ NA_integer_
    ),

    # ---- Race / Ethnicity (5 categories) ----
    # BRFSS: HISPANC3 (2019-2023, stable)
    #   1-4 = Hispanic subgroups (1=Mexican, 2=Puerto Rican, 3=Cuban, 4=Other Hispanic)
    #   5 = No (not Hispanic)
    #   7 = Don't know, 9 = Refused
    #
    # BRFSS: race_unified (MRACASC1 for 2019-2021/2023, MRACASC2 for 2022)
    #   10=White, 20=Black, 30=AIAN, 40=Asian, 41-47=Asian subgroups,
    #   50=NHPI, 51-54=NHPI subgroups, 60=Other, 77=DK, 88=None, 99=Refused
    #   Multi-race: concatenated codes (e.g., 1020 = White+Black)
    #
    # Hispanic takes priority. Then map race to 5 target categories.
    hisp_val = as.integer(HISPANC3),
    race_val = race_unified,

    race_eth_h = case_when(
      # Hispanic: any subgroup (1-4)
      hisp_val %in% 1:4                              ~ "Hispanic",
      # Not Hispanic (hisp_val == 5): classify by race
      hisp_val == 5L & race_val == 10L               ~ "White NH",
      hisp_val == 5L & race_val == 20L               ~ "Black NH",
      hisp_val == 5L & race_val %in% c(40:47)        ~ "Asian NH",
      # AIAN, NHPI, Other, or multi-race codes -> Other/Multi NH
      hisp_val == 5L & race_val %in% c(30L, 50:54, 60L) ~ "Other/Multi NH",
      hisp_val == 5L & race_val > 99L                ~ "Other/Multi NH",   # multi-race codes
      # DK/Refused for either -> NA
      TRUE ~ NA_character_
    ),

    # ---- Education (4 categories) ----
    # BRFSS _EDUCAG: stable across 2019-2023, already 4 categories
    #   1=Did not graduate high school
    #   2=Graduated high school
    #   3=Attended college or technical school
    #   4=Graduated from college or technical school
    #   9=Don't know/not sure/missing -> NA
    edu_val = as.integer(`_EDUCAG`),
    education_h = case_when(
      edu_val %in% 1:4 ~ edu_val,
      TRUE ~ NA_integer_
    ),

    # ---- Marital status (6 categories) ----
    # BRFSS MARITAL (stable 2019-2023):
    #   1=Married, 2=Divorced, 3=Widowed, 4=Separated,
    #   5=Never married, 6=Unmarried couple, 9=Refused
    # NOTE: BRFSS coding is the TARGET coding. NACC was recoded to match this.
    mar_val = as.integer(MARITAL),
    marital_h = case_when(
      mar_val %in% 1:6 ~ mar_val,
      mar_val == 9L     ~ NA_integer_,
      TRUE ~ NA_integer_
    ),

    # ---- Diabetes (binary) ----
    # diabetes_unified from 02_clean_data.R: coalesce(DIABCOR3, DIABETE4)
    # BRFSS coding: 1=Yes, 2=Yes but only during pregnancy (-> yes),
    #               3=No, 4=No/borderline/pre-diabetes (-> no),
    #               7=DK, 9=Refused -> NA
    diabetes_h = case_when(
      diabetes_unified %in% c(1L, 2L) ~ 1L,   # Yes (including pregnancy-only)
      diabetes_unified %in% c(3L, 4L) ~ 0L,   # No (including borderline)
      TRUE ~ NA_integer_                        # 7=DK, 9=Refused, NA
    ),

    # ---- Depression (binary) ----
    # depression_unified from 02_clean_data.R: coalesce(DEPRESS1, ADDEPEV3)
    # BRFSS YESNO: 1=Yes, 2=No, 7/77=DK, 9/99=Refused
    depression_h = case_when(
      depression_unified == 1L ~ 1L,
      depression_unified == 2L ~ 0L,
      TRUE ~ NA_integer_
    ),

    # ---- Heart attack / MI (binary) ----
    # mi_unified from 02_clean_data.R: coalesce(HEART2, CVDINFR4)
    # BRFSS YESNO: 1=Yes, 2=No, 7/77=DK, 9/99=Refused
    mi_h = case_when(
      mi_unified == 1L ~ 1L,
      mi_unified == 2L ~ 0L,
      TRUE ~ NA_integer_
    ),

    # ---- Angina / CVD (binary) ----
    # angina_unified from 02_clean_data.R: coalesce(ANGINA, CVDCRHD4)
    # BRFSS YESNO: 1=Yes, 2=No, 7/77=DK, 9/99=Refused
    angina_h = case_when(
      angina_unified == 1L ~ 1L,
      angina_unified == 2L ~ 0L,
      TRUE ~ NA_integer_
    ),

    # ---- Stroke (binary) ----
    # stroke_unified from 02_clean_data.R: coalesce(STROKE2, CVDSTRK3)
    # BRFSS YESNO: 1=Yes, 2=No, 7/77=DK, 9/99=Refused
    stroke_h = case_when(
      stroke_unified == 1L ~ 1L,
      stroke_unified == 2L ~ 0L,
      TRUE ~ NA_integer_
    ),

    # ---- Hypertension (binary) ----
    # htn_unified from 02_clean_data.R: coalesce(BPHIGH2, BPHIGH3, BPHIGH6)
    # BRFSS BPHIGHB coding: 1=Yes, 2=Yes only during pregnancy (-> yes),
    #                       3=No, 4=Borderline (-> no),
    #                       7/77=DK, 9/99=Refused, 96=Other -> NA
    # STRUCTURAL MISSINGNESS: Not asked in even years (2020, 2022).
    # htn_structural_missing flag carried from 02_clean_data.R.
    hypertension_h = case_when(
      htn_unified %in% c(1L, 2L) ~ 1L,   # Yes (including pregnancy-only)
      htn_unified %in% c(3L, 4L) ~ 0L,   # No (including borderline)
      htn_unified %in% c(7L, 9L, 77L, 96L, 99L) ~ NA_integer_,
      TRUE ~ NA_integer_                   # Includes structural missing from even years
    ),

    # ---- Outcomes: NA for BRFSS ----
    mci_h   = NA_integer_,
    ad_h    = NA_integer_,
    adrd_h  = NA_integer_,
    moca_h  = NA_real_,
    cdrsb_h = NA_real_
  ) %>%
  select(
    id_h, adrc_h, brfss_sampwt_h,
    age_h, male_h, race_eth_h, education_h, marital_h,
    diabetes_h, depression_h, mi_h, angina_h, stroke_h, hypertension_h,
    htn_structural_missing,
    mci_h, ad_h, adrd_h, moca_h, cdrsb_h,
    survey_year   # keep for imputation model (HTN structural missingness)
  )


# ============================================================================
# 4. COMBINE ADRC + BRFSS INTO SINGLE HARMONIZED DATASET
# ============================================================================

# Separate a1_zip before binding (ADRC only, for ACS linkage)
adrc_zip <- adrc_h %>% select(id_h, a1_zip)
adrc_h   <- adrc_h %>% select(-a1_zip)

# Add survey_year column to ADRC (NA — not applicable)
adrc_h <- adrc_h %>% mutate(survey_year = NA_integer_)

harmonized_all <- bind_rows(adrc_h, brfss_h) %>%
  mutate(
    # Convert race_eth_h to factor with consistent ordering (5 categories for Aim 1)
    race_eth_h = factor(race_eth_h,
                        levels = c("White NH", "Black NH", "Hispanic",
                                   "Asian NH", "Other/Multi NH")),
    # Convert education_h to labeled factor
    education_h = factor(education_h,
                         levels = 1:4,
                         labels = c("Less than HS", "HS diploma/GED",
                                    "Some college", "College graduate+")),
    # Convert marital_h to labeled factor
    marital_h = factor(marital_h,
                       levels = 1:6,
                       labels = c("Married", "Divorced", "Widowed",
                                  "Separated", "Never married",
                                  "Unmarried couple"))
  )


# ============================================================================
# 5. RESTRICT TO TARGET POPULATION AGE (D20)
# ============================================================================
# Primary analysis: adults 65+ (standard in dementia epidemiology, aligns with
#   published prevalence estimates from Matthews et al. 2019 and Alzheimer's
#   Association Facts & Figures, and enables case-control calibration analysis)
# Supplementary: adults 50+ (broader, retains more ADRC participants)
#
# The all-ages dataset is saved separately for Aim 1 descriptive comparison.

n_before <- nrow(harmonized_all)

# 50+ supplementary dataset
harmonized_50plus <- harmonized_all %>%
  filter(age_h >= 50 | is.na(age_h))
message("Age >= 50 filter: ", n_before, " -> ", nrow(harmonized_50plus),
        " (dropped ", n_before - nrow(harmonized_50plus), " rows)")
message("  ADRC: ", sum(harmonized_50plus$adrc_h == 1),
        "  BRFSS: ", sum(harmonized_50plus$adrc_h == 0))

# 65+ primary dataset
harmonized <- harmonized_all %>%
  filter(age_h >= 65 | is.na(age_h))   # keep NAs for imputation
n_after <- nrow(harmonized)
message("Age >= 65 filter: ", n_before, " -> ", n_after,
        " (dropped ", n_before - n_after, " rows)")
message("  ADRC: ", sum(harmonized$adrc_h == 1),
        "  BRFSS: ", sum(harmonized$adrc_h == 0))


# ============================================================================
# 6. CREATE 4-CATEGORY RACE VARIABLE FOR AIMS 2-3
# ============================================================================
# Aim 1 uses the full 5-category race_eth_h to show representativeness.
# Aims 2-3 (generalizability models) need stable race-stratified estimates,
# but ADRC has only 13 Black NH and 7 Other/Multi NH participants.
# Collapse to 4 categories: Hispanic (any race), NH Asian, NH White, NH Other.

harmonized <- harmonized %>%
  mutate(
    race_eth4_h = fct_collapse(race_eth_h,
      "NH Other" = c("Black NH", "Other/Multi NH")
    ) %>%
      fct_relevel("NH Other", after = Inf)  # put NH Other last
  )

message("\nrace_eth4_h (for Aims 2-3):")
message(paste(capture.output(
  table(harmonized$race_eth4_h, harmonized$adrc_h, useNA = "ifany",
        dnn = c("race_eth4_h", "adrc_h"))
), collapse = "\n"))


# ============================================================================
# 7. SAVE
# ============================================================================

# Full dataset (no age filter) for Aim 1 descriptive representativeness
saveRDS(harmonized_all, file.path(processed_data_dir,
                                   "03_adrc_brfss_harmonized_all_ages.rds"))

# Age 65+ primary analytic dataset for Aims 2-3 (imputation, weighting, AIPSW)
saveRDS(harmonized, file.path(processed_data_dir,
                               "03_adrc_brfss_harmonized.rds"))

# Age 50+ supplementary analytic dataset
saveRDS(harmonized_50plus, file.path(processed_data_dir,
                                      "03_adrc_brfss_harmonized_50plus.rds"))

# ADRC zip codes for ACS linkage
saveRDS(adrc_zip, file.path(processed_data_dir, "03_adrc_zip_for_acs.rds"))

message("\nSaved: 03_adrc_brfss_harmonized_all_ages.rds (",
        nrow(harmonized_all), " rows, no age filter, for Aim 1)")
message("Saved: 03_adrc_brfss_harmonized.rds (",
        sum(harmonized$adrc_h == 1), " ADRC + ",
        sum(harmonized$adrc_h == 0), " BRFSS = ",
        nrow(harmonized), " rows, age 65+, for Aims 2-3)")
message("Saved: 03_adrc_brfss_harmonized_50plus.rds (",
        nrow(harmonized_50plus), " rows, age 50+, supplementary)")
message("Saved: 03_adrc_zip_for_acs.rds")
message("\nHarmonized variables:")
message(paste(names(harmonized), collapse = ", "))
message("\nHTN structural missingness (BRFSS, age 65+):")
message(paste(capture.output(
  table(harmonized$htn_structural_missing[harmonized$adrc_h == 0],
        useNA = "ifany")
), collapse = "\n"))
