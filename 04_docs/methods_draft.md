# Methods — DRAFT

Working draft for manuscript methods section. Remaining flags for decisions/issues
marked with **[FLAG]**. References to decision log entries in parentheses (e.g., D01).

---

## Study design

This study combined multiple data sources to estimate the local burden and
demographic disparities in Alzheimer's disease (AD), AD and related dementias
(ADRD), and cognitive outcomes among adults aged 65 and older residing in
Alameda, San Mateo, and Santa Clara counties, California. The analysis had
three components:

1. **Local AD prevalence estimation** using published regression coefficients
   from the Chicago Health and Aging Project (CHAP; Dhana et al. 2023),
   augmented with locally-derived race/ethnicity coefficients from the Stanford
   ADRC, applied to individual-level population-representative data from the
   California BRFSS (D22, D24).

2. **Cognitive outcome disparity analysis** using a doubly-robust
   generalizability framework (AIPSW) to estimate demographically-standardized
   prevalence ratios, odds ratios, and mean differences in cognitive outcomes
   across racial/ethnic groups, comparing the Stanford ADRC source sample to
   the three-county target population represented by BRFSS (D08, D22).

3. **Clinical enrichment characterization** using referral source data and
   external prevalence benchmarks to quantify and describe the ADRC's
   outcome-dependent selection, as a transparency diagnostic (D22, D23).

The three-county target population was defined based on the Stanford ADRC's
primary catchment area. The target population is non-institutionalized adults
aged 65 and older, consistent with the BRFSS sampling frame.

## Data sources

### Stanford ADRC Clinical Core

The Stanford ADRC Clinical Core follows a standardized protocol for enrollment
and assessment of participants with normal cognition, mild cognitive impairment,
and dementia. Data are collected using the National Alzheimer's Coordinating
Center (NACC) Uniform Data Set (UDS) version 3 Initial Visit Packet. We
included all participants with a completed baseline visit.

The data extract was pulled in December 2024 (z1x_uds3_12092024.csv),
covering participants enrolled between September 2015 and December 2024.

The ADRC enrolls participants through multiple pathways: ADC clinician/staff
referrals (41%), self-referral (28%), registry or other research study
recruitment (15%), outside physician referrals (5%), and other sources (11%)
(D23). This recruitment strategy produces a sample that is clinically enriched
for cognitive impairment relative to the general population: approximately 52%
of participants aged 65+ had cognitive impairment at baseline, compared with
an estimated ~10% ADRD prevalence in the general 65+ population (Matthews et
al. 2019). This clinical enrichment is pervasive across all referral pathways:
ADC clinician referrals (69% impaired), outside physician referrals (69%),
self-referrals (38%), and research/registry recruitment (30%). The
implications of this enrichment for interpreting study results are discussed
in the Analysis sections below and in D22.

After restricting to adults aged 65+ (D20), the analytic ADRC sample included
529 participants.

### California BRFSS

The BRFSS is an annual, state-based telephone survey of non-institutionalized
adults conducted by the Centers for Disease Control and Prevention. We used
California BRFSS data from 2015-2024, restricted to respondents residing in
the three-county catchment area (Alameda: FIPS 001, San Mateo: FIPS 081,
Santa Clara: FIPS 085) and aged 65+ (D20). County was identified using the
BRFSS county FIPS codes, with COUNTY1 as the primary identifier and CTYCODE2
as a secondary source; these were harmonized via coalesce across survey years
(D14).

The 10-year BRFSS window (2015-2024) was chosen to align with the ADRC
enrollment period (September 2015 – March 2025) and to recover sample size
after the age 65+ restriction (D20). Demographic composition in the 3-county
catchment was stable across the decade (White NH: 72% → 70%; college
graduate+: ~61%; mean age: 73.5 → 74.5), and survey year was included as a
covariate in the propensity model to adjust for any temporal demographic
shifts.

After county and age restrictions, the BRFSS analytic sample included 2,150
respondents. BRFSS respondents were assigned their survey design weights
(_LLCPWT, the combined landline and cell phone weight), divided by 10 (the
number of pooled survey years) per CDC guidance for combining multiple years
of BRFSS data, to represent the non-institutionalized adult population of the
target area in a single year (D13).

### External prevalence data sources

Two external data sources provided population-based prevalence estimates used
for the hybrid model and enrichment characterization:

- **Dhana et al. (2023)** estimated AD prevalence across all US states and
  counties by applying regression coefficients from the Chicago Health and
  Aging Project (CHAP), a population-based cohort study, to 2020 Census
  population estimates. Their published model coefficients (intercept, age
  group, sex, race/ethnicity, education) are used in our hybrid prevalence
  model. A key limitation is the absence of an Asian category; Asian Americans
  were assumed to have the same risk as White Americans.

- **Matthews et al. (2019)** estimated ADRD prevalence among Medicare
  fee-for-service beneficiaries aged 65+ using the CMS Chronic Conditions
  Warehouse algorithm. Race/ethnicity-specific estimates (White 10.3%, Black
  13.8%, Hispanic 12.2%, Asian/Pacific Islander 8.4%) and age- and
  sex-specific estimates are used for enrichment ratio calculations and as
  an external benchmark.

## Covariate harmonization

We harmonized covariates across the NACC UDS v3 forms and BRFSS to enable
estimation of participation probabilities in the combined sample. The following
covariates were available in both data sources and included in the
participation model:

- **Age** (continuous): NACC — computed from date of birth to visit date.
  BRFSS — reported age (AGE variable, available all years).
- **Sex** (binary: male/female): NACC Form A1 Q4 (1=male, 2=female). BRFSS
  SEX (2014-2015) / SEX1 (2016-2017) / SEX2 (2018-2021) / _SEX (2020-2024)
  / BIRTHSEX (2022+), harmonized via coalesce. Same coding direction.
- **Race/ethnicity** (categorical): Constructed from NACC Hispanic ethnicity
  (A1 Q7) and race (A1 Q9) and from BRFSS Hispanic origin (HISPANC3) and race
  (MRACASC1/MRACASC2). Hispanic ethnicity took priority over race. For Aim 1
  (representativeness), a 5-category variable was used (White non-Hispanic,
  Black non-Hispanic, Hispanic, Asian non-Hispanic, Other/Multi non-Hispanic)
  to display the full distribution of race/ethnicity in the ADRC relative to
  the target population. For generalizability models, a 4-category variable
  was used: White non-Hispanic, Hispanic (any race), Asian non-Hispanic, and
  non-Hispanic Other (D02). The non-Hispanic Other category collapsed Black
  non-Hispanic (n=12 in ADRC 65+), American Indian/Alaska Native, Native
  Hawaiian/Pacific Islander, multiracial, and other races due to small cell
  sizes that precluded stable race-specific model estimation.

The collapsing of Black non-Hispanic participants into "NH Other" was driven
   by the small ADRC sample size (n=12 at age 65+), which precluded stable
   race-specific model estimation. This does not reflect the population
   composition of the catchment area — Alameda County in particular has a
   substantial Black population. Aim 1 representativeness analyses use the
   5-category race/ethnicity variable to make this recruitment disparity
   visible. Increasing Black participant enrollment is a priority for future
   ADRC recruitment efforts and is discussed in the Limitations.

- **Education** (4 categories): NACC years of education (A1 Q11: numeric,
  where 12 = HS/GED, 16 = bachelor's, 18 = master's, 20 = doctorate,
  99 = unknown) recoded to match BRFSS _EDUCAG: <12 years (no HS diploma),
  12 years (HS/GED), 13-15 years (some college), ≥16 years (college
  graduate+).
- **Marital status** (6 categories): NACC Form A1 Q12 recoded to match BRFSS
  MARITAL coding (NACC swaps codes for widowed=2 and divorced=3 relative to
  BRFSS widowed=3 and divorced=2).
- **Diabetes** (binary: ever diagnosed): NACC Form A5 Q5a (0=absent,
  1=recent/active, 2=remote/inactive; 1 or 2 coded as yes). BRFSS DIABCOR3
  (2015-2021) / DIABETE4 (2022-2024).
- **Depression** (binary: ever diagnosed, lasting 2+ years): NACC Form A5 Q7d.
  BRFSS DEPRESS1 (2015-2021) / ADDEPEV3 (2022-2024).
- **Stroke** (binary: ever): NACC Form A5 Q3a. BRFSS STROKE2 (2015-2021) /
  CVDSTRK3 (2022-2024).
- **Hypertension** (binary: ever): NACC Form A5 Q5b. BRFSS BPHIGH2
  (2015, 2017, 2019) / BPHIGH3 (2021) / BPHIGH6 (2023-2024). Hypertension
  was not ascertained in BRFSS even years (2016, 2018, 2020, 2022), resulting
  in structural missingness for approximately 50% of BRFSS respondents; this
  was addressed in the multiple imputation model using survey year as an
  auxiliary predictor (see below) (D05). A sensitivity analysis excluding
  hypertension from all models confirmed that results were robust to this
  structural missingness (see Sensitivity Analyses) (D16).

Unknown or refused responses in both data sources (NACC code 9; BRFSS codes
7/77 for "don't know" and 9/99 for "refused") were recoded to NA prior to
imputation (D12).

The outcome model additionally included myocardial infarction (MI) and angina
as predictors. These variables were available in both NACC (Form A5) and BRFSS
but were excluded from the participation model because they caused weight
instability (post-weighting SMDs flipped from small negative values to
problematic positive values due to low prevalence and heavily-weighted
individuals) (D06). Including different covariates in the participation and
outcome models is a feature of the doubly robust framework in which the AIPSW
estimator remains consistent if either model is correctly specified (Dahabreh
et al. 2019) (D08a).

The following covariates used in prior KHANDLE generalizability analyses
(Hayes-Larson et al. 2022) were not available in the NACC UDS v3 and were
omitted: individual income, self-rated health, exercise, difficulty walking,
difficulty dressing, and vision impairment (D07). Military/veteran status was
also omitted due to unreliable harmonization. The available NACC medical
history variables (diabetes, hypertension, stroke, depression) partially
substitute for self-rated health by providing objective health status
information, but unmeasured confounding from these omitted covariates remains
a limitation.

## Outcomes

Five cognitive outcomes were assessed at the ADRC (cross-sectional, using
baseline covariates with most recent visit cognitive status):

1. **Mild cognitive impairment (MCI)** (binary): Defined as any MCI subtype
   from NACC Form D1 (amnestic single-domain, amnestic multi-domain,
   non-amnestic single-domain, or non-amnestic multi-domain). Participants
   with normal cognition or dementia were coded as MCI = 0, as MCI represents
   an intermediate clinical state between normal cognition and dementia.

2. **Alzheimer's disease (AD)** (binary): Defined as AD etiology present from
   NACC Form D1 (d1_alzdis = 1, indicating AD is a contributing or primary
   etiologic diagnosis). Participants with normal cognition were coded as
   AD = 0. Participants with cognitive impairment but without an AD etiology
   designation were also coded as AD = 0. AD was retained alongside ADRD
   because it represents a clinically distinct etiology, and stakeholders
   frequently ask about AD specifically (D15, D21).

3. **Alzheimer's disease and related dementias (ADRD)** (binary): Defined as
   the presence of any of the following etiologic diagnoses from NACC Form D1:
   AD (d1_alzdis), Lewy body disease (d1_lbdis), cerebrovascular disease
   contributing to cognitive impairment (d1_cvd), or frontotemporal dementia
   (d1_ftldnos) (D21). This composite maps to the broad ADRD construct used in
   Medicare claims-based surveillance (Matthews et al. 2019). Among 274 ADRD
   cases in the ADRC (65+): AD=160, Lewy body=123, cerebrovascular=25, FTD=1.
   The high proportion of Lewy body disease (45% of ADRD cases) reflects
   Stanford ADRC's research focus and should be noted when comparing to
   population-based estimates.

4. **Montreal Cognitive Assessment (MoCA) score** (continuous, 0-30): Total
   score from NACC Form C2. Higher scores indicate better cognitive function.

5. **Clinical Dementia Rating Sum of Boxes (CDR-SB)** (continuous, 0-18):
   From NACC Form B4. Higher scores indicate greater impairment. CDR-SB has
   a substantial point mass at zero (participants with no clinical impairment)
   and was modeled using a two-part approach (see Statistical Analysis) (D09).

## Missing data

Missing covariate values were addressed using multiple imputation by chained
equations (MICE) with 40 imputed datasets and 10 iterations per dataset (D10).
ADRC and BRFSS data were imputed separately to avoid imposing the outcome
distribution (available only in ADRC) on BRFSS respondents and vice versa for
covariates available only in BRFSS (e.g., survey year). Outcome variables
(MCI, AD, ADRD, MoCA, CDR-SB) were included in the ADRC imputation model to
ensure compatibility between the imputation and analysis models.

MICE methods were auto-detected based on variable type: predictive mean
matching for continuous variables, logistic regression for binary variables,
and polytomous logistic regression for unordered categorical variables.

For BRFSS imputation, survey year was included as an auxiliary predictor (but
not imputed) to inform imputation of hypertension for even-year respondents
where hypertension was structurally missing (D05).

Convergence was assessed by examining between-imputation variability in
imputed means across the 40 datasets. Variability was minimal: for BRFSS
(which had the most missingness due to structural hypertension gaps), the
standard deviation of imputed means across imputations was 0.019 for age,
0.0007 for diabetes, and 0.012 for hypertension, indicating stable convergence
by 10 iterations.

## Statistical analysis

### Component 1: Local AD prevalence estimation (hybrid model)

To estimate AD prevalence in the three-county 65+ population, we used a hybrid
approach combining published population-based regression coefficients with
locally-derived race/ethnicity effects (D24).

Dhana et al. (2023) published logistic regression coefficients for AD
prevalence from the Chicago Health and Aging Project (CHAP), a population-based
cohort: intercept = -3.455, with age group (70-74: 0.577, 75-79: 1.126,
80-84: 1.800, 85+: 2.693), sex (female: 0.123), race/ethnicity (Black: 0.915,
Hispanic: 0.548), and education (per SD increase: -0.398) as predictors. A key
limitation of the Dhana model is the absence of an Asian category; Asian
Americans were assumed to have the same risk as White Americans.

We augmented the Dhana model by substituting ADRC-derived race/ethnicity
coefficients for the Hispanic and Asian non-Hispanic groups, while retaining
Dhana's intercept, age, sex, education, and Black coefficients. The ADRC race
coefficients are defensible under a case-control analogy: the ADRC is clinically
enriched for cognitive impairment, but if this enrichment is approximately
constant across racial/ethnic groups, odds ratios comparing race groups are
preserved (Breslow & Day 1980). We tested this assumption empirically by
computing enrichment ratios (ADRC prevalence / Matthews published prevalence)
by race and found approximately constant ratios (White NH: 4.8x, Hispanic:
3.8x, Asian NH: 4.5x, NH Other: 4.0x) (D22).

The hybrid model was applied to individual-level BRFSS data with survey
weights to produce survey-weighted AD prevalence estimates for the three-county
target population, overall and stratified by race/ethnicity, age group, and
sex. This approach improves on Dhana's county-level application by using
individual-level microdata and fills a critical gap for the Bay Area's large
Asian population.

For ADRD prevalence, we used Matthews et al. (2019) race × age × sex-specific
rates applied to the local population structure via indirect standardization,
as the ADRC's ADRD etiological composition (high Lewy body proportion) makes
direct ADRC-to-population inference less appropriate for ADRD than for AD.

### Component 2: Disparity analysis (AIPSW)

#### Participation model (inverse odds of selection weights)

We estimated the probability of ADRC participation (vs. being a BRFSS
respondent in the target population) using logistic regression fit with
`svyglm()` from the R `survey` package, with a `quasibinomial` family, on the
combined ADRC + BRFSS dataset (D03). BRFSS respondents carried their adjusted
survey design weights (_LLCPWT / 10); ADRC participants entered with weight = 1.
The participation model included the following covariates: race/ethnicity (4
categories), sex, age, education, marital status, diabetes, hypertension,
stroke, and depression, along with interactions between race/ethnicity and sex,
age, education, diabetes, hypertension, and depression.

Stabilized inverse odds of selection weights (sIOSW) were computed as:

  sw_i = [(1 - p_i) / p_i] * [P(S=1) / P(S=0)]

where p_i = P(S=1 | X_i) is the predicted participation probability and the
stabilization constant P(S=1) / P(S=0) used the survey-weighted marginal
probability of ADRC participation (D04).

#### Outcome models

Outcome models were fit in the ADRC sample only, predicting each cognitive
outcome from the same covariates as the participation model plus myocardial
infarction and angina (D08a). Model families were chosen based on outcome
distribution (D09, D18):

- **Binary outcomes** (MCI, AD, ADRD): logistic regression
- **MoCA** (continuous, 0-30): linear regression
- **CDR-SB** (continuous, 0-18, zero-inflated): two-part model. Part 1:
  logistic regression for P(CDR-SB > 0). Part 2: Gamma GLM with log link
  for E[CDR-SB | CDR-SB > 0]. The Gamma GLM was preferred over linear
  regression for the positive part based on AIC comparison (1703 vs 2045)
  and the right-skewed distribution of positive CDR-SB values (D18).

#### Doubly robust estimation (AIPSW)

The primary estimates were obtained using the augmented inverse probability of
selection weighting (AIPSW) estimator (Dahabreh et al. 2019; Buchanan et al.
2018) (D08):

  theta_DR = [survey-weighted mean of mu(X) in BRFSS]
             + [IOSW-weighted mean of (Y - mu(X)) in ADRC]

where mu(X) = E[Y | X, S=1] is the outcome model prediction. The first term
(g-computation) transports ADRC-fitted predictions to the target population
covariate distribution. The second term (residual correction) adjusts for any
remaining bias due to outcome model misspecification, weighted by the inverse
odds of selection.

**Interpretation note:** Because the ADRC is clinically enriched for cognitive
impairment (~52% impaired vs ~10% in the general 65+ population), AIPSW
absolute estimates reflect the demographically-standardized ADRC outcome
distribution, not population prevalence. The demographic reweighting adjusts
for differences in age, sex, race, education, and comorbidities between the
ADRC and target population, but cannot correct for outcome-dependent selection
(clinical enrichment). The primary inferential targets from this component are
therefore relative comparisons — prevalence ratios, odds ratios, and mean
differences across racial/ethnic groups — which are preserved under
approximately constant enrichment across groups (D22).

Race/ethnicity-stratified AIPSW estimates were computed by applying the same
doubly robust framework within each racial/ethnic group, using within-group
g-computation and within-group IOSW-weighted residual correction. Prevalence
ratios (PR) and prevalence/mean differences (PD) were computed relative to
White non-Hispanic as the reference group, using the bootstrap distributions
to obtain 95% confidence intervals (D17).

For comparison, we also report: (a) IOSW-only estimates (weighting without
outcome modeling), (b) g-computation-only estimates (outcome modeling without
weighting), and (c) unweighted ADRC estimates.

#### Variance estimation

Confidence intervals were obtained via stratified bootstrap (stratified by
ADRC/BRFSS membership to maintain the source/target ratio) with 1,000
replicates per imputed dataset (D11). Within each bootstrap replicate, the
participation model was refit using `glm()` with BRFSS weights normalized to
sum to n_BRFSS (preserving relative weighting while avoiding `svydesign()`
overhead) and the outcome model was refit on the resampled ADRC data (D03a).

Bootstrap estimates were pooled across 40 imputed datasets using Rubin's rules
(Rubin 1987; Barnard & Rubin 1999) (D19): the combined point estimate is the
mean of imputation-specific estimates, with total variance decomposed into
within-imputation variance (mean of bootstrap variances) and between-imputation
variance (variance of point estimates across imputations), with appropriate
degrees-of-freedom adjustment.

### Component 3: Clinical enrichment characterization

To quantify the ADRC's clinical enrichment and assess the plausibility of key
analytic assumptions, we conducted the following diagnostic analyses:

1. **Enrichment ratios by race/ethnicity:** We computed the ratio of ADRC ADRD
   prevalence (after demographic reweighting) to Matthews et al. (2019)
   published prevalence for each racial/ethnic group. Approximately constant
   enrichment ratios across groups support the assumption that race-specific
   odds ratios are preserved despite clinical enrichment (D22).

2. **Referral source analysis:** Using NACC Form A1 question 2a (a1_refersc),
   we characterized the distribution of referral pathways by cognitive status,
   education, and comorbidities to identify how clinical enrichment operates
   and whether it interacts with measured covariates (D23).

3. **ADRC ORs vs. Dhana published coefficients:** We compared ADRC-derived
   odds ratios for age, sex, education, and race/ethnicity against Dhana et
   al.'s population-based coefficients from CHAP. Under the case-control
   principle, odds ratios should be approximately equal despite the ADRC's
   enrichment (the intercept absorbs the enrichment), while substantial
   divergence would suggest enrichment-covariate interactions that compromise
   OR transportability.

### Covariate balance assessment

Covariate balance between the (reweighted) ADRC sample and the
survey-weighted BRFSS target population was assessed using standardized mean
differences (SMDs). SMDs were computed as
(mean_ADRC - mean_BRFSS) / SD_BRFSS for all covariates including dummy-coded
levels of categorical variables (race/ethnicity, education, marital status).
Balance was assessed both overall and stratified by race/ethnicity. Balance
was considered acceptable at |SMD| < 0.25, following prior generalizability
analyses (Stuart et al. 2011; Hayes-Larson et al. 2022). Before and after
weighting balance was displayed using Love plots.

### Sensitivity analyses

We conducted the following sensitivity analyses to assess robustness of
the AIPSW disparity estimates:

1. **Estimator comparison:** IOSW-only and g-computation-only estimates were
   compared to the primary AIPSW estimates. Agreement across estimators
   provides evidence against severe model misspecification.

2. **Hypertension exclusion (D16):** Because hypertension had structural
   missingness in BRFSS even years (2016, 2018, 2020, 2022), we re-fit both
   the participation and outcome models excluding hypertension and its
   interactions. All outcome estimates from the no-hypertension models fell
   within the 95% confidence intervals of the primary analysis, indicating
   that results are robust to this structural missingness.

3. **Propensity score overlap / positivity assessment:** Overlap of propensity
   score distributions between ADRC and BRFSS was assessed using density plots.
   Positivity diagnostics included the range of propensity scores in the ADRC
   sample and the proportion of extreme weights.

4. **Weight trimming sensitivity:** We compared AIPSW estimates using untrimmed
   stabilized IOSW weights to estimates with 99th percentile trimming
   (approximately 6 of 529 ADRC observations trimmed per imputation, max
   weight reduced from ~19 to ~7.7). Overall estimates changed by <5% for
   all outcomes, confirming that extreme weights do not distort aggregate
   results. Race-specific estimates for Hispanic participants were more
   sensitive to trimming (up to 22% for MCI), consistent with the small
   sample size (N=58) and reflected in the wide bootstrap confidence
   intervals for this group.

## Software

All analyses were conducted in R version 4.5.0 (2025-04-11) using the
following packages: `mice` 3.17.0 (multiple imputation), `survey` 4.4.2
(design-based propensity model), `boot` 1.3.31 (bootstrap inference),
`tidyverse` 2.0.0 (data management), and `ggplot2` 3.5.2 (visualization).

---

# Summary of remaining flags

All flags resolved.

### Flags resolved since previous draft

| Flag | Resolution |
|------|------------|
| D | Resolved — ADRC sampling frame described via referral source analysis (D23). ADRC recruits through multiple pathways: ADC clinician referral (41%), self-referral (28%), registry/research (15%), outside physician (5%), other (11%). Explicitly described as clinically enriched, not population-based. |
| E | County FIPS coding verified across all BRFSS years (D14). |
| F | BRFSS weights divided by 10 per CDC guidance for pooling 10 survey years (D13, D20). |
| I | Sensitivity analysis completed (D16). Excluding hypertension from all models changed estimates by 0.4-6.1%; all no-HTN estimates fell within main analysis 95% CIs. |
| J | Incorporated into text. MI/angina excluded from participation model (D06), retained in outcome model (D08a). Different covariate sets are a feature of the DR framework. |
| K | Clarified in text. MCI = 0 for dementia is clinically correct. Specified that outcomes are cross-sectional/prevalent. |
| L | Verified (D15). AD prevalence of 33% is correct for this enriched ADRC sample. |
| O | Reframed as a single model specification rather than claiming iterative building. |
| Q | Resolved — CDR-SB two-part model uses Gamma GLM with log link for positive part (D18). AIC comparison: Gamma 1703 vs linear 2045. |
| R | Implemented (D17). Prevalence ratios and differences computed from bootstrap distributions. |
| S | Resolved — Rubin's rules implemented for combining MI + bootstrap (D19). Point estimates averaged across imputations, total variance = within-imputation + (1 + 1/M) × between-imputation variance, with df adjustment. |
| T | Citation added (Stuart et al. 2011; Hayes-Larson et al. 2022 for 0.25 SMD threshold). |
| U | Sensitivity analysis section updated with four analyses including weight trimming. |
| C | Resolved — enrollment dates (September 2015 – December 2024) added from data extract filename and BRFSS alignment text. |
| N | Resolved — MICE convergence verified numerically. Between-imputation SD of means across 40 datasets is tiny (age: 0.019, diabetes: 0.0007, HTN: 0.012 for BRFSS). |
| P | Resolved — 99th percentile trimming changes overall AIPSW estimates by <5%. Hispanic race-specific estimates more sensitive (up to 22%) due to small N=58. Untrimmed weights retained for primary analysis (D25). |
| V | Resolved — R 4.5.0, mice 3.17.0, survey 4.4.2, boot 1.3.31, tidyverse 2.0.0, ggplot2 3.5.2. |
| A | Resolved — 3-county catchment confirmed as ADRC's primary recruitment area. Stated in study design. |
| B | Resolved — Target population defined as non-institutionalized adults 65+, consistent with BRFSS sampling frame. |
| G | Resolved — ACS not needed; removed. BRFSS is sufficient for all analyses. |
| H | Resolved — Black NH collapse framed as sample-size-driven, not analytic preference. 5-category variable used in Aim 1 to make disparity visible. Black recruitment flagged as priority in limitations. |
