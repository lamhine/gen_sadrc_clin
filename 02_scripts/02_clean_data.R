# 02_clean_data.R
# Purpose: Clean ADRC Clinical Core (NACC UDS v3) and California BRFSS data
#          for generalizability analysis.
#
# NACC variable references: UDS v3 Initial Visit Packet
#   Form A1 = Subject Demographics
#   Form A5 = Subject Health History
#   Form D1 = Clinician Diagnosis
#   Form C2 = Neuropsychological Battery (MoCA)
#   Form B4 = CDR Plus NACC FTLD
#
# BRFSS variable references: brfss_vars_summary.xlsx (by_year and by_var sheets)
#   Variable names change across survey years — see year-specific notes below.

# ---------------------- #
# LOAD CONFIGURATION
# ---------------------- #
source("config.R")

# ---------------------- #
# LOAD ADRC CLINICAL CORE AND CA BRFSS DATA
# ---------------------- #
adrc <- read_csv(file.path(adrc_data_dir, "z1x_uds3_12092025.csv"))
ca_bound <- readRDS(file = file.path(processed_data_dir, "01_ca_bound.rds"))


# ============================================================================
# ADRC CLINICAL CORE: CLEAN AND EXTRACT
# ============================================================================

# Create objective cognitive impairment (OCI) at the row/event level
adrc_event <- adrc %>%
  mutate(
    visit_num = as.integer(str_match(redcap_event_name, "^visit_(\\d+)_")[,2]),
    d1_normcog  = as.integer(d1_normcog),
    d1_demented = as.integer(d1_demented),
    # Form D1: Normal cognition vs. any cognitive impairment
    oci = case_when(
      d1_normcog == 1 ~ 0L,
      d1_normcog == 0 ~ 1L,
      TRUE ~ NA_integer_
    ),
    dementia = case_when(
      d1_demented == 1 ~ 1L,
      d1_demented == 0 ~ 0L,
      TRUE ~ NA_integer_
    )
  )

# Choose most recent visit with D1 observed
adrc_latest <- adrc_event %>%
  filter(!is.na(oci)) %>%
  group_by(record_id) %>%
  slice_max(order_by = visit_num, n = 1, with_ties = FALSE) %>%
  ungroup()

# Pull baseline covariates from A1 (demographics) and A5 (health history)
adrc_baseline <- adrc_event %>%
  filter(redcap_event_name == "visit_1_arm_1") %>%
  transmute(
    record_id,

    # ---- Form A1: Subject Demographics ----
    a1_visitdate = a1_visitdate,
    a1_birthyr   = as.integer(a1_birthyr),     # A1 Q3: Year of birth
    a1_birthmo   = as.integer(a1_birthmo),     # A1 Q3: Month of birth
    a1_sex       = as.integer(a1_sex),         # A1 Q4: 1=Male, 2=Female
    a1_educ      = as.integer(a1_educ),        # A1 Q11: Years of education (0-36, 99=unknown)
    a1_hispanic  = as.integer(a1_hispanic),    # A1 Q7: 0=No, 1=Yes Hispanic, 9=Unknown
    a1_hispor    = as.integer(a1_hispor),      # A1 Q8: Hispanic subgroup (1=Mex, 2=PR, etc.)
    a1_race      = as.integer(a1_race),        # A1 Q9: 1=White, 2=Black, 3=AIAN, 4=NHPI,
                                                #         5=Asian, 50=Other, 99=Unknown
    a1_racesec   = as.integer(a1_racesec),     # A1 Q9: Secondary race
    a1_raceter   = as.integer(a1_raceter),     # A1 Q9: Tertiary race
    a1_primlang  = as.integer(a1_primlang),    # A1 Q10: Primary language
    a1_maristat  = as.integer(a1_maristat),    # A1 Q12: 1=Married, 2=Widowed, 3=Divorced,
                                                #          4=Separated, 5=Never married,
                                                #          6=Living as married, 9=Unknown
    a1_zip       = a1_zip,                     # A1 Q17: First 3 digits of ZIP

    # ---- Form A5: Subject Health History ----
    # Tobacco (A5 Q1a-d)
    a5_tobac30   = as.integer(a5_tobac30),     # A5 Q1a: Tobacco in last 30 days (0=No, 1=Yes)
    a5_tobac100  = as.integer(a5_tobac100),    # A5 Q1b: Smoked 100+ cigarettes (0=No, 1=Yes)
    a5_smokyrs   = as.integer(a5_smokyrs),     # A5 Q1c: Years smoked
    a5_quitsmok  = as.integer(a5_quitsmok),    # A5 Q1d: Year quit smoking

    # Cardiovascular (A5 Q2a, Q2g, Q3a)
    # All coded: 0=absent, 1=recent/active, 2=remote/inactive, 9=unknown
    a5_cvhatt    = as.integer(a5_cvhatt),      # A5 Q2a: Heart attack/cardiac arrest
    a5_cvangina  = as.integer(a5_cvangina),    # A5 Q2g: Angina/CVD
    a5_cbstroke  = as.integer(a5_cbstroke),     # A5 Q3a: Stroke

    # Other conditions (A5 Q5a, Q5b, Q7d)
    # All coded: 0=absent, 1=recent/active, 2=remote/inactive, 9=unknown
    a5_diabetes  = as.integer(a5_diabetes),    # A5 Q5a: Diabetes
    a5_hyperten  = as.integer(a5_hyperten),    # A5 Q5b: Hypertension
    a5_hypercho  = as.integer(a5_hypercho),    # A5 Q5c: Hypercholesterolemia (auxiliary)
    a5_dep2yrs   = as.integer(a5_dep2yrs)      # A5 Q7d: Depression (2+ years)
  )

# Merge baseline covariates (visit 1) to outcomes from latest visit with D1
adrc_person <- adrc_baseline %>%
  inner_join(
    adrc_latest %>% select(
      record_id, visit_num, oci, dementia,
      # Form D1: MCI subtypes for composite MCI indicator
      d1_normcog, d1_mciamem, d1_mciaplus, d1_mcinon1, d1_mcinon2,
      # Form D1: Dementia etiologies (for ADRD composite)
      d1_alzdis, d1_alzdisif, d1_lbdis, d1_cvd, d1_ftldnos,
      # Form C2: MoCA total score (0-30)
      c2_mocatots,
      # Form B4: CDR Sum of Boxes (0-18) and global CDR
      b4_cdrsum, b4_cdrglob
    ),
    by = "record_id"
  )


# ============================================================================
# BRFSS: FILTER TO CATCHMENT AND TARGET YEARS
# ============================================================================

# Create unified survey year
# IYEAR is populated for 2023+ data; YEAR is populated for 2014-2022 data
ca_bound <- ca_bound %>%
  mutate(
    survey_year = coalesce(as.integer(IYEAR), as.integer(YEAR))
  )

# Filter to catchment area (Alameda=001, San Mateo=081, Santa Clara=085)
catchment <- ca_bound %>%
  mutate(
    COUNTY1_clean  = ifelse(COUNTY1  %in% c(777, 888, 999, 116), NA_integer_, as.integer(COUNTY1)),
    CTYCODE2_clean = ifelse(CTYCODE2 %in% c(777, 888, 999, 116), NA_integer_, as.integer(CTYCODE2)),
    COUNTY_clean   = ifelse(COUNTY   %in% c(777, 999),           NA_integer_, as.integer(COUNTY)),

    county_fips = coalesce(
      COUNTY1_clean,
      CTYCODE2_clean,
      case_when(
        COUNTY_clean == 1  ~ 1L,    # Alameda -> 001
        COUNTY_clean == 41 ~ 81L,   # San Mateo -> 081
        COUNTY_clean == 43 ~ 85L,   # Santa Clara -> 085
        TRUE ~ NA_integer_
      )
    ),

    county_name = case_when(
      county_fips == 1  ~ "Alameda",
      county_fips == 81 ~ "San Mateo",
      county_fips == 85 ~ "Santa Clara",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(county_fips %in% c(1L, 81L, 85L))

# Filter BRFSS to 2015-2024 (aligned with ADRC enrollment period)
brfss_catchment <- catchment %>%
  filter(survey_year %in% 2015:2024)


# ============================================================================
# BRFSS: UNIFY YEAR-SPECIFIC VARIABLE NAMES
# ============================================================================
# BRFSS renames variables across years. Create unified columns using coalesce
# so downstream scripts can use a single name per construct.
#
# Variable name changes documented in brfss_vars_summary.xlsx (by_year sheet):
#   Sex:        SEX (2014-2015), SEX1 (2016-2017), SEX2 (2018-2021),
#               BIRTHSEX (2022), _SEX (2020-2024 imputed)
#   Stroke:     STROKE2 (2015-2021) -> CVDSTRK3 (2022-2024)
#   HTN:        BPHIGH2 (2015,2017,2019), BPHIGH3 (2021), BPHIGH6 (2023,2024)
#               NOT ASKED in 2016, 2018, 2020, 2022 -> structural missingness
#   Diabetes:   DIABCOR3 (2015-2021) -> DIABETE4 (2022-2024)
#   Depression: DEPRESS1 (2015-2021) -> ADDEPEV3 (2022-2024)
#   MI:         HEART2 (2015-2021) -> CVDINFR4 (2022-2024)
#   Angina/CVD: ANGINA (2015-2021) -> CVDCRHD4 (2022-2024)
#   Race:       MRACASC1 (2015-2021, 2023-2024) -> MRACASC2 (2022 only)
#   Hispanic:   HISPANC3 (all years, stable)
#   Education:  _EDUCAG (all years, stable)

brfss_catchment <- brfss_catchment %>%
  mutate(
    # Sex: coalesce across year-specific names; all use 1=Male, 2=Female
    # _SEX (imputed) available 2020+; SEX2 for 2018-2021; SEX1 for 2016-2017; SEX for 2014-2015
    sex_unified = coalesce(
      as.integer(`_SEX`),
      as.integer(SEX2),
      as.integer(SEX1),
      as.integer(SEX),
      as.integer(BIRTHSEX)
    ),

    # Stroke: STROKE2 (2015-2021) / CVDSTRK3 (2022-2024)
    # BRFSS YESNO: 1=Yes, 2=No, 7/77=DK, 9/99=Refused, BLANK=NA
    stroke_unified = coalesce(as.integer(STROKE2), as.integer(CVDSTRK3)),

    # Hypertension: BPHIGH2 (2015,2017,2019) / BPHIGH3 (2021) / BPHIGH6 (2023,2024)
    # NOT asked in 2016, 2018, 2020, 2022 -> structural missingness
    # Note: 2024 restored HTN (BPHIGH6), breaking the even-year pattern
    # BPHIGHB coding: 1=Yes, 2=Yes during pregnancy, 3=No, 4=Borderline,
    #                 7/77=DK, 9/99=Refused, 96=Other, BLANK=NA
    htn_unified = coalesce(
      as.integer(BPHIGH2),
      as.integer(BPHIGH3),
      as.integer(BPHIGH6)
    ),

    # Diabetes: DIABCOR3 (2015-2021) / DIABETE4 (2022-2024)
    # 1=Yes, 2=Yes pregnancy, 3=No, 4=Borderline/pre-DM, 7=DK, 9=Refused
    diabetes_unified = coalesce(as.integer(DIABCOR3), as.integer(DIABETE4)),

    # Depression: DEPRESS1 (2015-2021) / ADDEPEV3 (2022-2024)
    # YESNO: 1=Yes, 2=No, 7/77=DK, 9/99=Refused
    depression_unified = coalesce(as.integer(DEPRESS1), as.integer(ADDEPEV3)),

    # Heart attack/MI: HEART2 (2015-2021) / CVDINFR4 (2022-2024)
    mi_unified = coalesce(as.integer(HEART2), as.integer(CVDINFR4)),

    # Angina/CVD: ANGINA (2015-2021) / CVDCRHD4 (2022-2024)
    angina_unified = coalesce(as.integer(ANGINA), as.integer(CVDCRHD4)),

    # Race: MRACASC1 (2015-2021, 2023-2024) / MRACASC2 (2022 only)
    # Both use same coding: 10=White, 20=Black, 30=AIAN, 40=Asian, etc.
    race_unified = coalesce(as.integer(MRACASC1), as.integer(MRACASC2)),

    # Flag structural missingness for hypertension (not asked in 2016, 2018, 2020, 2022)
    # Note: 2024 DOES have HTN (BPHIGH6), breaking the even-year pattern
    htn_structural_missing = survey_year %in% c(2016, 2018, 2020, 2022)
  )


# ============================================================================
# SAVE CLEANED DATA
# ============================================================================
saveRDS(adrc_person,     file.path(processed_data_dir, "02_adrc_clean.rds"))
saveRDS(brfss_catchment, file.path(processed_data_dir, "02_brfss_catchment.rds"))

message("Saved: 02_adrc_clean.rds (", nrow(adrc_person), " participants)")
message("Saved: 02_brfss_catchment.rds (", nrow(brfss_catchment),
        " BRFSS records, years ",
        paste(sort(unique(brfss_catchment$survey_year)), collapse = ", "), ")")
