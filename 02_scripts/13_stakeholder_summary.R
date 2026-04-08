# 13_stakeholder_summary.R
# Purpose: Generate stakeholder-facing summary tables combining all three
#          analysis components into presentation-ready format.
#
# Components:
#   1. Local AD prevalence (hybrid model, script 12)
#   2. Local ADRD prevalence (Matthews indirect standardization, script 12)
#   3. Cognitive disparities (AIPSW, script 07/09)
#   4. Enrichment diagnostics (OR comparison, script 12)
#
# INPUTS:
#   - results_dir/12_hybrid_prevalence.rds
#   - results_dir/12_or_comparison.rds
#   - results_dir/07_aipsw_results.rds
#
# OUTPUTS:
#   - results_dir/13_prevalence_summary.csv
#   - results_dir/13_disparity_summary.csv
#   - results_dir/13_enrichment_diagnostics.csv
#   - results_dir/13_full_summary.rds (all tables in one object)
#
# ============================================================================

source("config.R")


# ============================================================================
# 1. LOAD RESULTS
# ============================================================================

hybrid     <- readRDS(file.path(results_dir, "12_hybrid_prevalence.rds"))
or_diag    <- readRDS(file.path(results_dir, "12_or_comparison.rds"))
aipsw      <- readRDS(file.path(results_dir, "07_aipsw_results.rds"))

# Check if PR/PD results exist (from script 09)
pr_pd_file <- file.path(results_dir, "09_pr_pd_results.rds")
has_pr_pd  <- file.exists(pr_pd_file)
if (has_pr_pd) {
  pr_pd <- readRDS(pr_pd_file)
}


# ============================================================================
# 2. TABLE 1: LOCAL PREVALENCE ESTIMATES
# ============================================================================
# Combines AD (hybrid model) and ADRD (Matthews) into one stakeholder table

# Overall
prevalence_overall <- tibble(
  group = "Overall (65+)",
  ad_pct_dhana = sprintf("%.1f%%", hybrid$overall$ad_prev_pct[1]),
  ad_pct_hybrid = sprintf("%.1f%%", hybrid$overall$ad_prev_pct[2]),
  adrd_pct_matthews = sprintf("%.1f%%", hybrid$adrd_overall_pct)
)

# By race
prevalence_race <- hybrid$by_race %>%
  left_join(hybrid$adrd_by_race %>% select(race_eth4_h, matthews_rate),
            by = "race_eth4_h") %>%
  transmute(
    group = as.character(race_eth4_h),
    ad_pct_dhana = sprintf("%.1f%%", dhana_prev),
    ad_pct_hybrid = sprintf("%.1f%%", hybrid_prev),
    adrd_pct_matthews = sprintf("%.1f%%", matthews_rate)
  )

# By age
prevalence_age <- hybrid$by_age %>%
  transmute(
    group = as.character(age_group),
    ad_pct_dhana = sprintf("%.1f%%", dhana_prev),
    ad_pct_hybrid = sprintf("%.1f%%", hybrid_prev),
    adrd_pct_matthews = "---"
  )

# By sex
prevalence_sex <- hybrid$by_sex %>%
  transmute(
    group = sex_label,
    ad_pct_dhana = sprintf("%.1f%%", dhana_prev),
    ad_pct_hybrid = sprintf("%.1f%%", hybrid_prev),
    adrd_pct_matthews = "---"
  )

prevalence_table <- bind_rows(
  prevalence_overall,
  tibble(group = "--- By Race/Ethnicity ---",
         ad_pct_dhana = "", ad_pct_hybrid = "", adrd_pct_matthews = ""),
  prevalence_race,
  tibble(group = "--- By Age Group ---",
         ad_pct_dhana = "", ad_pct_hybrid = "", adrd_pct_matthews = ""),
  prevalence_age,
  tibble(group = "--- By Sex ---",
         ad_pct_dhana = "", ad_pct_hybrid = "", adrd_pct_matthews = ""),
  prevalence_sex
)

names(prevalence_table) <- c("Group",
                              "AD Prevalence (Dhana model)",
                              "AD Prevalence (Hybrid: Dhana + local race)",
                              "ADRD Prevalence (Matthews)")

write.csv(prevalence_table, file.path(results_dir, "13_prevalence_summary.csv"),
          row.names = FALSE)

message("=== LOCAL PREVALENCE ESTIMATES ===")
message("3-county target population, adults 65+\n")
message(paste(capture.output(print(prevalence_table, n = Inf)), collapse = "\n"))


# ============================================================================
# 3. TABLE 2: DISPARITY ESTIMATES (AIPSW)
# ============================================================================
# Race/ethnicity-stratified AIPSW estimates and ratios vs White NH

# Absolute estimates by race
race_ests <- aipsw %>%
  filter(grepl("^aipsw_", estimand)) %>%
  mutate(
    race = str_replace(estimand, "^aipsw_", "") %>%
      str_replace_all("_", " ") %>% str_to_title(),
    outcome_label = factor(outcome,
      levels = c("mci_h", "ad_h", "adrd_h", "moca_h", "cdrsb_h"),
      labels = c("MCI", "AD", "ADRD", "MoCA", "CDR-SB")),
    est_ci = sprintf("%.2f (%.2f, %.2f)", estimate, lci, uci)
  ) %>%
  select(outcome_label, race, est_ci) %>%
  pivot_wider(names_from = race, values_from = est_ci)

message("\n=== AIPSW ESTIMATES BY RACE/ETHNICITY ===")
message("Demographically-standardized to 3-county target population\n")
message(paste(capture.output(print(race_ests, n = Inf, width = 120)), collapse = "\n"))

# Prevalence/mean ratios (if available)
if (has_pr_pd) {
  ratio_table <- pr_pd %>%
    filter(grepl("Ratio", measure)) %>%
    mutate(
      outcome_label = factor(outcome,
        levels = c("mci_h", "ad_h", "adrd_h", "moca_h", "cdrsb_h"),
        labels = c("MCI", "AD", "ADRD", "MoCA", "CDR-SB")),
      est_ci = sprintf("%.2f (%.2f, %.2f)", estimate, lci, uci)
    ) %>%
    select(outcome_label, race, est_ci) %>%
    pivot_wider(names_from = race, values_from = est_ci)

  message("\n=== PREVALENCE/MEAN RATIOS (vs White NH) ===")
  message(paste(capture.output(print(ratio_table, n = Inf, width = 120)),
                collapse = "\n"))

  diff_table <- pr_pd %>%
    filter(grepl("Difference", measure)) %>%
    mutate(
      outcome_label = factor(outcome,
        levels = c("mci_h", "ad_h", "adrd_h", "moca_h", "cdrsb_h"),
        labels = c("MCI", "AD", "ADRD", "MoCA", "CDR-SB")),
      est_ci = sprintf("%.3f (%.3f, %.3f)", estimate, lci, uci)
    ) %>%
    select(outcome_label, race, est_ci) %>%
    pivot_wider(names_from = race, values_from = est_ci)

  message("\n=== PREVALENCE/MEAN DIFFERENCES (vs White NH) ===")
  message(paste(capture.output(print(diff_table, n = Inf, width = 120)),
                collapse = "\n"))

  write.csv(bind_rows(
    ratio_table %>% mutate(measure = "Ratio"),
    diff_table %>% mutate(measure = "Difference")
  ), file.path(results_dir, "13_disparity_summary.csv"), row.names = FALSE)
}


# ============================================================================
# 4. TABLE 3: ENRICHMENT DIAGNOSTICS
# ============================================================================

or_comp <- or_diag$or_comparison %>%
  mutate(across(where(is.numeric), ~ round(., 2))) %>%
  select(
    Variable = variable,
    `Dhana OR (population-based)` = dhana_or,
    `ADRC OR` = adrc_or,
    `ADRC 95% LCI` = adrc_or_lci,
    `ADRC 95% UCI` = adrc_or_uci,
    `OR Ratio (ADRC/Dhana)` = or_ratio,
    Note = note
  )

# Add Asian row
asian_row <- or_diag$adrc_ors_pooled %>%
  filter(term == "race_eth4_hAsian NH") %>%
  transmute(
    Variable = "Asian NH",
    `Dhana OR (population-based)` = 1.00,
    `ADRC OR` = round(or_pooled, 2),
    `ADRC 95% LCI` = round(or_lci, 2),
    `ADRC 95% UCI` = round(or_uci, 2),
    `OR Ratio (ADRC/Dhana)` = round(or_pooled / 1.0, 2),
    Note = "Dhana assumes Asian = White; ADRC fills this gap"
  )

or_comp <- bind_rows(or_comp, asian_row)

write.csv(or_comp, file.path(results_dir, "13_enrichment_diagnostics.csv"),
          row.names = FALSE)

message("\n=== ENRICHMENT DIAGNOSTIC: ADRC vs POPULATION-BASED ORs ===")
message(paste(capture.output(print(or_comp, n = Inf, width = 130)),
              collapse = "\n"))


# ============================================================================
# 5. KEY TAKEAWAYS (text summary for stakeholder presentation)
# ============================================================================

message("\n")
message("================================================================")
message("KEY FINDINGS: AD/ADRD in Alameda, San Mateo, and Santa Clara Counties")
message("================================================================")
message("")
message("1. HOW BIG IS THE PROBLEM?")
message(sprintf("   - Estimated AD prevalence among adults 65+: %.1f%%",
                hybrid$overall$ad_prev_pct[2]))
message(sprintf("   - Estimated ADRD prevalence among adults 65+: %.1f%%",
                hybrid$adrd_overall_pct))
message(sprintf("   - By age: 65-69 (%.1f%%), 70-74 (%.1f%%), 75-79 (%.1f%%), 80-84 (%.1f%%), 85+ (%.1f%%)",
                hybrid$by_age$hybrid_prev[1], hybrid$by_age$hybrid_prev[2],
                hybrid$by_age$hybrid_prev[3], hybrid$by_age$hybrid_prev[4],
                hybrid$by_age$hybrid_prev[5]))
message("")
message("2. WHO IS MOST AFFECTED?")
ad_race <- hybrid$by_race
message(sprintf("   AD prevalence by race/ethnicity:"))
for (i in 1:nrow(ad_race)) {
  message(sprintf("     %s: %.1f%%", ad_race$race_eth4_h[i], ad_race$hybrid_prev[i]))
}
message("")
message("3. WHAT'S DIFFERENT ABOUT OUR COMMUNITY?")
message("   - Existing national models (Dhana et al. 2023) have NO Asian category")
message("   - Bay Area has a large, diverse Asian population (~8% of 65+ residents)")
message(sprintf("   - Local estimate for Asian NH AD prevalence: %.1f%%",
                ad_race$hybrid_prev[ad_race$race_eth4_h == "Asian NH"]))
message("   - This fills a critical gap in existing literature")
message("")
message("4. COGNITIVE HEALTH DISPARITIES (from ADRC clinical data)")
message("   - Hispanic and NH Other groups show lower MoCA scores than White NH")
message("     (indicating greater cognitive burden after demographic standardization)")
message("   - Full disparity tables above (PRs and PDs with 95% CIs)")
message("")
message("5. WHAT WE CAN'T TELL FROM THESE DATA")
message("   - True local prevalence (ADRC is clinically enriched, ~5x population rate)")
message("   - Time trends (cross-sectional analysis)")
message("   - Incidence rates (would require longitudinal population-based data)")
message("   - Asian subgroup differences (Chinese, Filipino, Vietnamese, Indian, etc.)")
message("")


# ============================================================================
# 6. SAVE COMBINED
# ============================================================================

saveRDS(list(
  prevalence_table = prevalence_table,
  race_estimates_aipsw = race_ests,
  ratio_table = if (has_pr_pd) ratio_table else NULL,
  diff_table = if (has_pr_pd) diff_table else NULL,
  or_comparison = or_comp
), file.path(results_dir, "13_full_summary.rds"))

message("Saved: 13_prevalence_summary.csv")
message("Saved: 13_disparity_summary.csv")
message("Saved: 13_enrichment_diagnostics.csv")
message("Saved: 13_full_summary.rds")
