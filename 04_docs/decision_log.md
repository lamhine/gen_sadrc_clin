# Decision Log: Stanford ADRC Generalizability Analysis

Tracks analytic decisions for manuscript methods section. Each entry records
what was decided, why, and what alternatives were considered.

---

## D01: Target population age restriction (age 50+)

**Date:** 2026-03-24
**Decision:** Restrict both ADRC and BRFSS samples to adults aged 50+.
**Rationale:** Without restriction, the BRFSS sample includes all adults 18+
(mean age ~50, SD ~18) while the ADRC sample has mean age ~72. This creates
severe positivity violations — half the BRFSS sample has no ADRC counterpart.
Age 50+ retains all but 12 ADRC participants and defines a target population
relevant to ADRC's mission (cognitive aging/dementia).
**Alternatives considered:** Age 55+ or 65+ would improve overlap but reduce
BRFSS sample size for race-stratified estimates. Age 50+ balances overlap
and statistical power.
**Impact:** ADRC: 665 -> 653 participants. BRFSS: 4,182 -> 2,112 respondents.

---

## D02: Race/ethnicity variable — two versions for different aims

**Date:** 2026-03-24
**Decision:** Use 5-category race/ethnicity (White NH, Black NH, Hispanic,
Asian NH, Other/Multi NH) for Aim 1 representativeness, and 4-category
(White NH, Hispanic, Asian NH, NH Other) for Aims 2-3 generalizability models.
**Rationale:** Aim 1 should show the full representativeness picture, including
the very small Black NH group (n=13 in ADRC). For Aims 2-3, Black NH (n=13)
and Other/Multi NH (n=7) are too small for stable race-stratified estimates
or interactions in the propensity model. Collapsing to "NH Other" (n=20)
provides a single group for model stability, though estimates for this group
will have wide CIs and should be interpreted cautiously.
**Alternatives considered:** Dropping small groups entirely; using a 3-category
variable (White NH, Asian NH, Other). The 4-category approach preserves
Hispanic as a distinct group, which is important given the heavily Hispanic
catchment area.

---

## D03: Survey weights in propensity model — svyglm() with design-based weights

**Date:** 2026-03-24
**Decision:** Use `svyglm()` with `quasibinomial` family and a `svydesign()`
object for the primary propensity score estimates (script 05). ADRC
participants enter the combined design with weight = 1 (census of themselves);
BRFSS respondents carry their `_LLCPWT` survey weights.
**Rationale:** The BRFSS weights are probability/sampling weights representing
the complex survey design. Using `glm(weights = brfss_sampwt_h)` treats them
as frequency weights, which produced propensity scores near zero for all ADRC
participants (P(ADRC) ~ 0.0003) and absurdly large IOSW (median ~1,100).
`svyglm()` properly handles design-based weights, and combined with
population-level stabilization (using survey-weighted marginal P(ADRC) rather
than unweighted sample proportion), produces stabilized IOSW with mean ~1
and reasonable range (0.03 to 27).
**Alternatives considered:** (a) Raw survey weights in `glm()` — produced
broken weights. (b) Normalizing BRFSS weights to sum to n_BRFSS before
`glm()` — produces similar point estimates and is used inside the bootstrap
for computational speed (see D03a). `svyglm()` is the more principled
approach for the primary analysis.
**Impact:** Stabilized IOSW: min=0.03, median=0.33, mean=0.95, max=27.
The KHANDLE analyses (Repos 1-2) used `glm(weights=)` with raw survey
weights; this may warrant a note if comparing approaches directly.

### D03a: Bootstrap uses glm() with normalized weights for speed

**Date:** 2026-03-24
**Decision:** Inside the bootstrap (script 07), use `glm()` with BRFSS
weights normalized to sum to n_BRFSS, rather than `svyglm()`.
**Rationale:** `svyglm()` + `svydesign()` per bootstrap replicate is ~5x
slower than `glm()` with normalized weights. With 50 bootstrap reps x 40
imputations x 4 outcomes = 8,000 model fits, this is the difference between
~10 minutes and ~50+ minutes (and laptop sleep kills long runs). Normalizing
BRFSS weights preserves relative weighting across respondents while keeping
the ADRC:BRFSS ratio at the sample level. Point estimates are nearly
identical; the bootstrap itself handles variance estimation.
**Impact:** Runtime reduced from ~50 min to ~10 min. Final bootstrap CIs
will use 1,000 reps (currently 50 for development).

---

## D04: Stabilization of IOSW — population-level marginal probability

**Date:** 2026-03-24
**Decision:** Stabilize the inverse odds weights using the survey-weighted
marginal P(ADRC), not the unweighted sample proportion.
**Rationale:** The stabilized IOSW formula is:
  sw_i = [(1 - p_i) / p_i] * [P(S=1) / P(S=0)]
When the propensity model uses survey weights, p_i reflects the
population-level probability. The stabilization constant must be on the same
scale: P(S=1) = sum(ADRC weights) / sum(all weights) ~ 0.00009, not the
unweighted 653/2765 ~ 0.24. Using the unweighted proportion produced weights
in the thousands; the population-level proportion produces weights with
mean ~ 1.
**Impact:** This is a direct consequence of D03. The two must be consistent.

---

## D05: Hypertension structural missingness in BRFSS even years

**Date:** 2026-03-24
**Decision:** Include hypertension in the participation model despite
structural missingness in BRFSS 2020 and 2022 (the question was not asked).
Handle via multiple imputation with survey year as an auxiliary predictor.
**Rationale:** Hypertension is an important health covariate for the
participation model (prevalence differs between ADRC and catchment).
1,595 of 2,112 BRFSS respondents (76% of even-year respondents) have
structural missing HTN. The imputation model includes `survey_year` as a
predictor so that HTN is imputed using odd-year data and correlated
covariates (age, diabetes, stroke, etc.).
**Alternatives considered:** (a) Drop hypertension entirely — loses an
important health indicator. (b) Restrict BRFSS to odd years only — cuts
sample size substantially. (c) Current approach retains full sample and
leverages MI to fill the gap.
**Note for methods:** This is missingness by survey design, not individual
nonresponse. Should be described as structural missingness, not MCAR/MAR.

---

## D06: Covariates in participation model — iterative specification

**Date:** 2026-03-24
**Decision (initial):** Three-stage iterative propensity model:
  M1: race, sex, age, education
  M2: + marital status, diabetes, hypertension, stroke
  M3: + depression, MI, angina (FINAL)
All models include race x covariate interactions.

**Update (pending):** Remove MI and angina from participation model based
on covariate balance diagnostics. Retain them in outcome model only.

**Rationale for update:** Post-weighting SMD diagnostics showed MI and angina
flipping from small negative SMDs (-0.07, -0.13) to problematic positive
SMDs (+0.32, +0.20). Both have low prevalence (~2-4% in ADRC) and the few
heavily-weighted individuals with MI/angina (e.g., one White NH participant
with weight 26.9, MI=1, angina=1) disproportionately affect the weighted
distribution. Removing them from the participation model eliminates this
instability; they remain in the outcome model (script 06) where they
contribute to predicting cognitive outcomes.

**Alternatives considered:** (a) Weight trimming at 95th percentile — would
help but treats the symptom, not the cause. (b) Keeping MI/angina with
additional interactions — increases model complexity with sparse data.
(c) Current approach: remove from participation model, retain in outcome model.

**Confirmed covariates in participation model (revised):**
  age_h, male_h, race_eth4_h, education_h, marital_h,
  diabetes_h, hypertension_h, stroke_h, depression_h
  + race interactions with: male, age, education, diabetes, hypertension,
    depression

**Covariates in outcome model only (not in participation model):**
  mi_h, angina_h (plus all participation model covariates)

---

## D07: Gap variables — omitted from participation model

**Date:** 2026-03-24
**Decision:** The following KHANDLE covariates are NOT available in NACC UDS
v3 and are omitted from the participation model:
  - Individual income (NACC has only 3-digit ZIP; no individual income)
  - Self-rated health (GENHLTH — not in NACC)
  - Exercise (EXERANY — not in NACC)
  - ADL walking difficulty (DIFFWALK — no NACC equivalent; do not proxy)
  - ADL dressing difficulty (DIFFDRES — no NACC equivalent; do not proxy)
  - Vision/blindness (BLIND — not in NACC)
  - Military/veteran status (unreliable in NACC; high missingness)

**Rationale:** These variables cannot be harmonized between NACC and BRFSS.
Proxying (e.g., using CDR functional domains for ADLs) would introduce
measurement error and is not recommended. The available NACC medical history
variables (diabetes, hypertension, stroke, depression) partially compensate
for the self-rated health gap by providing objective health status information.

**Impact on interpretation:** The participation model may be misspecified if
these omitted variables independently predict both ADRC participation and
cognitive outcomes. This is a key limitation. The doubly robust estimator
provides some protection: if the outcome model is correctly specified,
consistent estimation is achieved even if the participation model is not.

**For Aim 1 only:** ZIP-level ACS median income linkage (via a1_zip) may be
used as a supplementary descriptive comparison, not as a model covariate.

---

## D08: Doubly robust AIPSW as primary estimator

**Date:** 2026-03-24
**Decision:** Use augmented inverse probability of selection weighting (AIPSW)
as the primary estimator, with IOSW-only and g-computation-only as
sensitivity/comparison analyses.
**Rationale:** The AIPSW estimator is consistent if either the participation
model or the outcome model is correctly specified. Given the gap variables
(D07), the participation model is likely somewhat misspecified. The outcome
model provides a second chance at consistent estimation. The three estimators
are reported side-by-side to assess robustness.
**AIPSW formula:**
  theta_DR = [survey-weighted mean of mu(X) in BRFSS]
             + [IOSW-weighted mean of (Y - mu(X)) in ADRC]
where mu(X) = E[Y | X, S=1] is the outcome model fitted in ADRC, and IOSW
are the stabilized inverse odds of selection weights.

---

## D08a: Participation and outcome models use different covariate sets

**Date:** 2026-03-24
**Decision:** The participation model and outcome model intentionally use
different covariate sets. The participation model includes 9 covariates
(age, sex, race, education, marital status, diabetes, hypertension, stroke,
depression). The outcome model includes these 9 plus MI and angina (11 total).
**Rationale:** This is a feature of the doubly robust framework, not a
limitation. The participation model is fit on the combined ADRC + BRFSS
sample, so it can only use harmonizable variables, and must produce stable
weights. The outcome model is fit only in ADRC, so it can include any
variable that improves prediction of cognitive outcomes — MI and angina are
clinically relevant predictors of dementia/cognitive impairment. Dahabreh
et al. (2019) and Buchanan et al. (2018) both note that the outcome model
can (and often should) include predictors beyond those in the participation
model. The DR consistency property holds as long as either model is correctly
specified, regardless of whether they share the same covariates.

---

## D09: Four outcomes with different modeling approaches

**Date:** 2026-03-24
**Decision:**
  1. MCI diagnosis (binary) -> logistic outcome model
  2. AD diagnosis (binary) -> logistic outcome model
  3. MoCA score (continuous, 0-30) -> linear outcome model
  4. CDR-SB (continuous, 0-18, zero-inflated) -> two-part model
     Part 1: P(CDR-SB > 0) via logistic regression
     Part 2: E[CDR-SB | CDR-SB > 0] via linear regression
**Rationale:** CDR-SB has a large mass at zero (cognitively normal
participants) with a right-skewed positive tail. A single linear model
would be misspecified. The two-part model separately estimates the
probability of any impairment and the severity given impairment.

---

## D10: Multiple imputation — 40 imputations, separate ADRC/BRFSS

**Date:** 2026-03-24
**Decision:** m = 40 imputations, maxit = 10, ADRC and BRFSS imputed
separately then combined. MICE with auto-detected methods (PMM for continuous,
logreg for binary, polyreg for unordered factors).
**Rationale:** Separate imputation avoids imposing the ADRC outcome
distribution on BRFSS (outcomes are structurally absent in BRFSS) and vice
versa for gap variables. m = 40 provides stable pooled estimates.
**Note:** Survey year is included as an auxiliary predictor in BRFSS
imputation to help predict hypertension for even-year respondents (see D05).

---

## D11: Bootstrap variance estimation — stratified by study membership

**Date:** 2026-03-24
**Decision:** Use stratified bootstrap (stratified by adrc_h) for confidence
intervals. Development runs: 50 reps. Final analysis: 1,000 reps.
**Rationale:** Stratification maintains the ADRC:BRFSS ratio in each
bootstrap sample. Without stratification, some replicates could
under-sample ADRC participants, causing model fitting failures.
Standard errors are pooled across 40 imputations using Rubin's rules
(implicitly via the bootstrap-then-pool approach: bootstrap within each
imputation, stack all bootstrap draws, take percentiles).

---

## D12: NACC "unknown" (9) and BRFSS "don't know/refused" (7/9) treated as NA

**Date:** 2026-03-25
**Decision:** All NACC unknown codes (9 for categorical variables, 99 for
education) and all BRFSS missing codes (7/77 = don't know, 9/99 = refused,
BLANK = not asked) are recoded to `NA` during harmonization (script 03) and
subsequently imputed by MICE (script 04).
**Rationale:** In this analysis context, "unknown" carries no meaningful
analytic information distinct from "missing." For the participation model,
keeping "unknown" as a separate factor level would create an ADRC-only
category (since BRFSS has no equivalent "unknown" level), distorting
propensity scores. For imputation, coding 9 as a valid value would cause MICE
to treat it as a real response rather than a value to predict from other
covariates.
**Exception:** `a1_educ = 9` in NACC means "9 years of education" (9th grade),
not "unknown." The NACC codebook uses 99 (not 9) for unknown education. This
value is correctly retained and mapped to education category 1 (less than HS).
**Verified:** All harmonized binary variables contain only 0, 1, and NA.
No residual 9s, 99s, 7s, or 77s remain in the harmonized dataset.
**Affected variables (NACC):** a1_hispanic (2 unknowns), a1_maristat (1),
a5_diabetes (3), a5_dep2yrs (5), a5_cvhatt (2), a5_cbstroke (4),
a5_hyperten (3). All small counts, handled appropriately by MICE.

---

## D13: BRFSS survey weights divided by number of pooled years

**Date:** 2026-03-25
**Decision:** Divide each BRFSS respondent's `_LLCPWT` weight by 10 (the number
of pooled survey years: 2015-2024) when constructing `brfss_sampwt_h`.
**Rationale:** Each year of BRFSS weights sums to approximately the full
population of the target area (~1.4M for the 3 counties). Pooling 10 years
without adjustment inflates the effective target population by 10x,
which distorts the population-level P(ADRC) used in weight stabilization
(D04) and the g-computation denominator. Dividing by N_years is standard
CDC guidance for combining multiple BRFSS years (CDC, "Complex Sampling
Weights and Preparing Module Data for Analysis," 2022).
**Implementation:** `brfss_sampwt_h = _LLCPWT / N_BRFSS_YEARS` in
script 03_harmonize.R. `N_BRFSS_YEARS = 10` is defined in config.R.
**Note:** The CDC also recommends making strata and PSU identifiers unique
across years (concatenating year with _STSTR and _PSU). We do not currently
use stratified survey design in the bootstrap (D03a), but this should be
addressed if we switch to a full survey design for variance estimation.
**Update (2026-03-25):** Changed from 5 years (2019-2023) to 10 years
(2015-2024) per D20.

---

## D14: County FIPS verification across BRFSS years

**Date:** 2026-03-25
**Decision:** Verified that county identification uses COUNTY1 (FIPS format,
available 2019-2022) as primary source and CTYCODE2 (FIPS format, available
2022-2023) as fallback. The COUNTY variable uses sequential CA county codes
in some years and FIPS codes in others (e.g., 2021 switched mid-series),
so it is used only as a tertiary fallback with an explicit crosswalk.
**Verification results:**
  2019: 1,114 respondents (Alameda=428, San Mateo=195, Santa Clara=491)
  2020: 563 respondents (196, 113, 254)
  2021: 688 respondents (293, 121, 274)
  2022: 1,032 respondents (444, 168, 420)
  2023: 785 respondents (310, 142, 333)
  Total: 4,182 respondents across 5 years
Proportions across counties are stable year-to-year and consistent with
relative population sizes. No action needed.

---

## D15: AD prevalence (33%) reflects ADRC enriched enrollment

**Date:** 2026-03-25
**Observation:** AD prevalence in the ADRC sample is 33% (218/665), which
appears high but is expected for a clinical ADRC that enrolls both healthy
controls and cognitively impaired participants. The sample is approximately
48% normal cognition (n=320) and 52% impaired (n=345). Among the impaired,
63% have AD etiology (218/345). The remaining 127 impaired participants have
d1_alzdis=NA (not d1_alzdis=0), indicating the NACC form only records AD
when present. These are coded as ad_h=0 (no AD).
**Note for methods:** AD prevalence reflects the clinical enrichment of the
ADRC sample, not community prevalence. The generalizability analysis
reweights this to estimate what AD prevalence would be in the target
population. The unweighted 33% should not be interpreted as a population
estimate.

---

## D16: HTN sensitivity analysis — results robust to HTN exclusion

**Date:** 2026-03-25
**Decision:** Retain hypertension in the primary analysis models but report
sensitivity analysis showing results are robust to its exclusion.
**Rationale:** Hypertension has ~76% structural missingness in BRFSS because
the question (RFHYPE6/RFHYPE5) is only asked in odd survey years (2019, 2021,
2023). Even after multiple imputation with survey_year as an auxiliary predictor,
the quality of imputed HTN values was uncertain. A sensitivity analysis re-fit
both the participation model and outcome models without hypertension_h (and its
race interaction) and compared AIPSW point estimates:
  - MCI:    0.294 (main) vs. 0.300 (no HTN) — 2.3% difference, within main CI
  - AD:     0.368 (main) vs. 0.360 (no HTN) — 2.1% difference, within main CI
  - MoCA:  20.869 (main) vs. 20.778 (no HTN) — 0.4% difference, within main CI
  - CDR-SB: 2.906 (main) vs. 3.083 (no HTN) — 6.1% difference, within main CI
All no-HTN estimates fall within the main analysis 95% CIs, indicating that
the structural missingness of HTN does not meaningfully affect conclusions.
**Script:** `02_scripts/10_sensitivity_htn.R`

---

## D17: Prevalence ratios and differences (vs. White NH reference)

**Date:** 2026-03-25
**Decision:** Compute prevalence/mean ratios (PR) and differences (PD) for each
non-White-NH racial/ethnic group relative to White NH, using the bootstrap
distributions from the main AIPSW analysis.
**Method:** For each bootstrap replicate b:
  PR_b = est_race_b / est_white_nh_b
  PD_b = est_race_b - est_white_nh_b
Point estimates are means across all bootstrap × imputation draws; 95% CIs
are 2.5th/97.5th percentiles. This preserves the correlation structure between
race-stratified estimates within each bootstrap replicate.
**Script:** `02_scripts/09_results_summary.R` (Section 4)

---

## D18: CDR-SB two-part model — Gamma GLM for positive part

**Date:** 2026-03-25
**Decision:** Use Gamma GLM with log link (instead of linear regression) for the
second part of the CDR-SB two-part model (E[CDR-SB | CDR-SB > 0]).
**Rationale:** CDR-SB positive values are right-skewed (skewness = 1.39,
mean = 4.6, median = 2.5, range 0.5-18). Model comparison on imputation 1:
  - Linear regression: AIC = 2044.9, BIC = 2170.8
  - Gamma GLM (log link): AIC = 1703.4, BIC = 1829.3
Gamma is preferred by ~340 AIC points — a decisive difference. The Gamma GLM
also naturally constrains predictions to be positive, avoiding the need for
`pmax(pred, 0)` flooring.
**Scripts:** `02_scripts/06_outcome_models.R`, `02_scripts/07_aipsw_analysis.R`

---

## D19: Rubin's rules for combining MI + bootstrap

**Date:** 2026-03-25
**Decision:** Use Rubin's rules (Rubin 1987; Barnard & Rubin 1999) to combine
point estimates and variance across 40 imputed datasets, replacing the
stack-and-percentile approach.
**Method:**
  1. Within each imputation m: theta_m = mean of bootstrap draws (point est),
     V_m = variance of bootstrap draws (within-imputation variance)
  2. Combined point estimate: theta_bar = mean(theta_m)
  3. Within-imputation variance: W_bar = mean(V_m)
  4. Between-imputation variance: B = var(theta_m)
  5. Total variance: T = W_bar + (1 + 1/M) * B
  6. df = (M-1)(1 + W_bar / ((1+1/M)*B))^2
  7. CI: theta_bar +/- t_df * sqrt(T)
**Rationale:** Properly decomposes uncertainty into sampling variability
(within-imputation) and missing data uncertainty (between-imputation), with
appropriate degrees-of-freedom adjustment. The previous stack-and-percentile
approach is simpler and used by Hayes-Larson et al. (2022), but Rubin's rules
is the methodologically preferred standard.
**Note:** Bootstrap draws are still stored stacked for prevalence ratio
computation (D17), which requires the joint distribution within replicates.
**Script:** `02_scripts/07_aipsw_analysis.R`

---

## D20: BRFSS expanded to 2015-2024 (10 years) and analysis restricted to age 65+

**Date:** 2026-03-25
**Decision:** (a) Expand BRFSS target population years from 2019-2023 to
2015-2024 (10 years). (b) Restrict the primary analytic sample to adults
aged 65+ (supplementary analysis retains 50+).
**Rationale for year expansion:** Restricting to 65+ reduces the BRFSS
catchment sample from ~2,100 (50+, 5 years) to ~1,045 (65+, 5 years).
Expanding to 10 years recovers sample size (~2,000 at 65+) and aligns the
BRFSS window with the ADRC enrollment period (Sept 2015 – Mar 2025).
Demographic composition in the 3-county catchment is stable across the
decade (White NH: 72% → 70%, education: ~61% college+, mean age: 73.5 → 74.5).
survey_year is included as a covariate in the propensity model to adjust for
any temporal shifts.
**Rationale for 65+ restriction:** (a) Standard age cutoff in dementia
epidemiology; directly comparable to published prevalence estimates (Matthews
et al. 2019, Alzheimer's Association Facts & Figures). (b) The ADRC recruits
healthy volunteers at age 70+, so 50-64 participants enter almost exclusively
through clinical referral — creating a group with near-100% clinical selection
that cannot be corrected by demographic reweighting. (c) Enables the planned
case-control calibration analysis (secondary), which requires external
prevalence estimates available only for 65+.
**Variable availability for 2015-2024:** All harmonized variables have
coverage across the full range via coalesce of year-specific names:
  - Sex: SEX (2015) / SEX1 (2016-17) / SEX2 (2018-21) / _SEX (2020-24)
  - Stroke: STROKE2 (2015-21) / CVDSTRK3 (2022-24)
  - Diabetes: DIABCOR3 (2015-21) / DIABETE4 (2022-24)
  - Depression: DEPRESS1 (2015-21) / ADDEPEV3 (2022-24)
  - MI: HEART2 (2015-21) / CVDINFR4 (2022-24)
  - Angina: ANGINA (2015-21) / CVDCRHD4 (2022-24)
  - HTN: BPHIGH2 (2015,17,19) / BPHIGH3 (2021) / BPHIGH6 (2023-24)
    Structural missingness in even years (2016,2018,2020,2022,2024) = 50%
**Impact:** N_BRFSS_YEARS updated from 5 to 10 in config.R. Year filter in
02_clean_data.R updated from 2019:2023 to 2015:2024. HTN structural
missingness flag extended to all even years. Sex coalesce chain extended
with SEX1 and SEX.
**Scripts modified:** config.R, 02_scripts/02_clean_data.R

---

## D21: ADRD added as fifth outcome

**Date:** 2026-04-05
**Decision:** Add all-cause Alzheimer's disease and related dementias (ADRD)
as a fifth outcome, defined from NACC Form D1 as: AD (d1_alzdis) OR Lewy body
disease (d1_lbdis) OR cerebrovascular disease (d1_cvd) OR frontotemporal
dementia (d1_ftldnos).
**Rationale:** ADRD maps directly to the Medicare claims-based ADRD definition
used in Matthews et al. (2019), enabling comparison with published national
prevalence estimates. The ADRC has AD-only and etiology-specific diagnoses, but
ADRD is the standard construct in dementia epidemiology and the most relevant
for stakeholder communication.
**ADRC composition (65+):** 274 ADRD cases (51.8%): AD=160, Lewy body=123,
cerebrovascular=25, FTD=1. Note: Lewy body is unusually high (45% of ADRD
cases), reflecting Stanford ADRC's research focus areas.
**Impact:** Five outcomes: MCI, AD, ADRD, MoCA, CDR-SB. All five run through
the AIPSW pipeline (scripts 06-07). AD is retained alongside ADRD because it
represents a clinically distinct etiology, and many stakeholders specifically
ask about AD.
**Scripts modified:** 02_clean_data.R (added d1_lbdis, d1_cvd, d1_ftldnos to
ADRC select), 03_harmonize.R (adrd_h outcome), 04_multiple_imputation.R
(added to imputation variables), 06_outcome_models.R (logistic model),
07_aipsw_analysis.R (added to outcomes list), 09_results_summary.R (labels).

---

## D22: Analysis reframed — disparity analysis primary, prevalence via external models

**Date:** 2026-04-06
**Decision:** Reframe the analysis from "ADRC-based population prevalence
estimation" to a three-component analysis:

  1. **Local AD/ADRD prevalence via external models** (Dhana et al. 2023
     coefficients applied to BRFSS microdata, augmented with ADRC-derived
     race coefficients): answers "how big is the problem here?"

  2. **AIPSW disparity analysis** (demographically-standardized prevalence
     ratios, odds ratios, and cognitive score gaps across race/ethnicity):
     answers "who is most affected?"

  3. **Clinical enrichment characterization** (enrichment ratios, referral
     source analysis, OR comparisons with published population-based
     estimates): transparency diagnostic that honestly characterizes what
     the ADRC can and cannot tell us.

**Rationale:** The ADRC's clinical enrichment (~52% impaired vs ~10% in the
general 65+ population) means that AIPSW with demographic-only reweighting
cannot recover population prevalence — selection depends on the outcome
(Y → S), violating the conditional exchangeability assumption. This was
explored extensively:

  - Intercept-shift calibration using Matthews et al. (2019) prevalence was
    implemented but abandoned: it was partly circular (overall estimate
    matched Matthews by construction) and relied on strong, unverifiable
    assumptions (conditional independence of demographic and clinical
    selection, transportability of all ADRC model coefficients).
  - Restricting to self-referred "volunteers" was considered but rejected:
    self-referrals are still 40% ADRD (4x population rate), and education
    variation within the ADRC is too narrow (~16-17 years mean across all
    subgroups) to reliably estimate education effects.
  - Multiplying IOSW by clinical calibration weights was considered but
    has no established methodological framework.

However, the ADRC's race ORs ARE defensible under the case-control principle:
if clinical enrichment is approximately constant across race groups (tested
empirically: White NH 4.8x, Hispanic 3.8x, Asian NH 4.5x, NH Other 4.0x),
odds ratios comparing race groups are preserved. This supports both the
disparity analysis and the hybrid prevalence model.

**Key references for reframing:**
  - Dhana et al. (2023): Published regression coefficients for AD prevalence
    from CHAP, applied to county-level Census data. We apply to individual-
    level BRFSS microdata.
  - Matthews et al. (2019, PMC6333531): Medicare FFS claims-based ADRD
    prevalence by race, age, sex. Used for enrichment ratio calculations
    and as external benchmark.
  - Mayeda et al. (2016, 2017): KPNC Northern California dementia incidence
    by race/ethnicity including Asian subgroups. Local validation reference.
**Alternatives considered:** (a) Full Bayesian evidence synthesis across
data sources — methodologically appealing but beyond pilot scope. (b) Two-
phase sampling framework (Breslow & Chatterjee 1999) — requires known Phase 2
selection probabilities, which the ADRC doesn't have. (c) Reporting AIPSW
absolute estimates as prevalence with caveats — rejected as misleading given
the known ~5x enrichment.

---

## D23: Referral source characterization of clinical enrichment

**Date:** 2026-04-06
**Decision:** Use NACC Form A1 question 2a (a1_refersc: principal referral
source) to characterize the clinical enrichment mechanism. This is a
descriptive/diagnostic analysis, not a correction.
**NACC referral source codes:**
  1 = Self-referral
  2 = Non-professional contact (spouse/partner, relative, friend, etc.)
  3 = ADC participant referral
  4 = ADC clinician, staff, or investigator referral
  5 = Nurse, doctor, or other health care provider
  6 = Other research study clinician/staff/investigator (non-ADC)
  8 = Other
  9 = Unknown
**Key findings (65+, N=531 with non-missing referral source):**
  - ADC clinician/staff referral (code 4): N=216 (41%), 69.4% impaired,
    67.6% ADRD — primary source of clinical enrichment
  - Self-referral (code 1): N=148 (28%), 37.8% impaired, 39.9% ADRD —
    moderately enriched (still 4x population rate)
  - ADC participant referral (code 3): N=26 (5%), 30.8% impaired — least
    enriched
  - Outside physician (code 5): N=26 (5%), 69.2% impaired
  - Registry/other research (codes 5-6): N=81 (15%), ~35% impaired
**Education by referral pathway:** Uniformly high across all pathways
(mean 16.4-17.2 years regardless of referral source or cognitive status).
This confirms that education-based self-referral bias exists but is in the
X → S pathway (handled by IOSW), not differentially interacting with Y → S.
However, the narrow education range means the ADRC cannot reliably estimate
education-outcome associations across the full population range.
**Implication:** Clinical enrichment is pervasive — even self-referrals are
4x enriched. Dropping any single referral pathway does not solve the problem.
The referral source analysis characterizes the mechanism transparently but
does not enable correction.
**Script:** Ad hoc analysis (referral_explore.R), findings reported in
methods/discussion.

---

## D24: Hybrid prevalence model — Dhana base + ADRC race coefficients

**Date:** 2026-04-06
**Decision:** Estimate local AD prevalence using a hybrid model that combines
published regression coefficients from Dhana et al. (2023, CHAP-based) with
ADRC-derived race/ethnicity coefficients, applied to individual-level BRFSS
microdata.

**Model specification:**
  logit P(AD | X) = Dhana intercept (-3.455)
                  + Dhana age coefficients (70-74: 0.577, 75-79: 1.126,
                    80-84: 1.800, 85+: 2.693)
                  + Dhana sex coefficient (female: 0.123)
                  + Dhana education coefficient (per SD: -0.398)
                  + ADRC race coefficients (Hispanic, Asian NH vs White NH)

**Rationale for hybrid:**
  - Dhana's intercept, age, sex, and education coefficients come from CHAP,
    a population-based cohort with good variation on all covariates. These
    are well-estimated and transportable.
  - Dhana has NO Asian category (assumes Asian = White). Matthews shows
    Asian/PI prevalence (8.4%) is lower than White (10.3%). The Bay Area
    has a large, diverse Asian population — this gap matters.
  - Dhana's Black coefficient is from a Chicago cohort. Our collapsed
    "NH Other" category (N=18, mostly Black) is too small for stable
    estimation. Use Dhana's Black coefficient for this group.
  - ADRC race ORs for Hispanic (N=58) and Asian NH (N=40) are defensible
    under the case-control principle: enrichment ratios are approximately
    constant across race groups (3.8x-4.8x), so ORs are approximately
    preserved. The ADRC provides LOCAL race coefficients for a Bay Area
    population that may differ from Chicago (Dhana) or national Medicare
    (Matthews).
**Limitations:**
  - ADRC race ORs estimated from modest samples (Hispanic N=58, Asian N=40)
  - Case-control OR preservation assumes enrichment is constant across race
    groups — tested empirically (3.8x-4.8x) but with wide uncertainty given
    sample sizes
  - Education coefficient from Dhana (Chicago) may not transfer perfectly to
    Bay Area, but the ADRC cannot provide a better estimate (insufficient
    variation)
  - Model is for AD only (Dhana's outcome). ADRD prevalence uses Matthews
    rates via indirect standardization as the ADRC's ADRD composition is
    unusual (48% Lewy body, reflecting Stanford research focus)
**Script:** 12_hybrid_prevalence.R

---

## D25: Trimmed weight sensitivity analysis (FLAG-P)

**Date:** 2026-04-06
**Decision:** 99th percentile trimming of stabilized IOSW has minimal impact
on overall AIPSW estimates (<5% change for all outcomes). Race-specific
estimates for Hispanic are sensitive to trimming (MCI: 22%, CDR-SB: 18%
change), consistent with small sample size (N=58). White NH, Asian NH, and
NH Other estimates are stable (<4% change).
**Rationale:** ~6 of 529 ADRC observations per imputation are trimmed, with
max weight dropping from ~19 to ~7.7. The overall estimates barely move,
confirming that extreme weights are not distorting aggregate results. The
Hispanic sensitivity is expected given small N and is already reflected in
the wide bootstrap confidence intervals.
**Conclusion:** Untrimmed weights are appropriate for the primary analysis.
Hispanic-specific estimates should be interpreted with caution (noted in
limitations). No need to use trimmed weights as a primary specification.
**Script:** 10_sensitivity_trimmed_weights.R
