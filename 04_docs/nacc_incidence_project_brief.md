# Project Brief: Generalizable Dementia Incidence Rates by Race/Ethnicity from the National NACC UDS

*Draft brief for a new analysis project. Paste into a new Claude session to get started.*

---

## Overview

This project will estimate **generalizable incidence rates of MCI and dementia by race/ethnicity** in the US older adult population, using longitudinal data from the National Alzheimer's Coordinating Center (NACC) Uniform Data Set (UDS) reweighted to the national Behavioral Risk Factor Surveillance System (BRFSS) target population. A key contribution is the first reported incidence rates for **multiracial adults** (n=1,340 in NACC, 985 aged 65+).

## Background and motivation

- The NACC UDS contains longitudinal cognitive assessments from ~40,000+ participants across 30+ Alzheimer's Disease Research Centers (ADRCs), but the sample is a convenience/volunteer sample that overrepresents White, highly educated, and healthy individuals
- Chan et al. (2025, doi:10.1002/alz.70657) documented substantial between-center heterogeneity in NACC participant characteristics and proposed (but did not implement) a meta-analytic framework treating each ADRC as its own study
- Hayes-Larson et al. (2022, doi:10.1002/alz.12491) demonstrated generalizability weighting for the single-site KHANDLE study using IOSW and BRFSS as the target population
- Our Stanford ADRC project (gen_sadrc_clin repo) extended this to doubly robust AIPSW estimation for a single site
- **No study has yet implemented a multi-site meta-analytic generalizability framework across NACC ADRCs**
- NACC has 1,340 multiracial participants -- one of the largest samples of multiracial older adults in any cognitive aging study. Most studies collapse multiracial into "Other" due to small cell sizes. This project can report race-specific incidence rates for a 6-category race/ethnicity variable.

## Proposed approach

### Option A: Pooled NACC with site effects
1. Pool all NACC sites with longitudinal follow-up
2. Include ADRC site as a random effect in the participation model
3. Compute IOSW for each NACC participant relative to national BRFSS (adults 65+)
4. Fit IOSW-weighted Cox proportional hazards model for time-to-MCI and time-to-dementia
5. Generalize marginal hazard/survival curves to the US target population
6. Stratify by 6-category race/ethnicity (White NH, Black NH, Hispanic, Asian NH, Multiracial NH, Other NH)

### Option B: Multi-site meta-analysis (preferred, novel)
1. For each ADRC site k: identify the local catchment population using state/county BRFSS
2. Compute site-specific IOSW: P(BRFSS_k | X) / P(ADRC_k | X)
3. Estimate site-specific incidence rates by race/ethnicity using weighted survival models
4. Meta-analyze across sites using random-effects meta-analysis (DerSimonian-Laird or REML)
5. Report pooled incidence rates with heterogeneity statistics (I^2, tau^2)
6. Explore sources of heterogeneity (site recruitment practices, geography, urban/rural)

### Doubly robust extension
- Outcome model: parametric survival model (Weibull or flexible parametric) fit within NACC
- Participation model: logistic regression predicting P(NACC | X) vs P(BRFSS | X)
- AIPSW for survival: combine weighted Kaplan-Meier/Cox with outcome model predictions
- Reference: Dahabreh et al. (2019) Section 4 on survival outcomes

## Data sources

### Source sample: NACC UDS
- Request full national UDS dataset from NACC (all sites, all visits)
- Key forms: A1 (demographics), A5 (health history), D1 (clinical diagnosis), C2 (MoCA), B4 (CDR)
- Inclusion: participants with normal cognition at baseline (D1 normal cognition diagnosis) and >= 1 follow-up visit
- Outcome: time from baseline normal cognition to first MCI or dementia diagnosis
- Censoring: last visit date if no progression; death date if available
- Site identifier: ADRC center ID

### Target population: BRFSS
- National BRFSS (not just California), pooled 2019-2023
- For Option A: national BRFSS adults 65+
- For Option B: state-level or county-level BRFSS matching each ADRC's catchment
- Survey weights: _LLCPWT / N_years_pooled (per CDC guidance)

### Supplementary: ACS PUMS
- For catchment area demographics and validation
- PUMA-to-county crosswalks for defining site catchment areas

## Covariate harmonization (NACC <-> BRFSS)
Same framework as the Stanford ADRC project. Key harmonized covariates:
- Age, sex, race/ethnicity (6 categories), education (4 categories)
- Marital status, diabetes, depression, stroke, hypertension
- Gaps (not in BRFSS): income, self-rated health, exercise, ADL difficulties

## Key methodological considerations

1. **Left truncation:** participants enter at different ages; condition on survival to study entry
2. **Informative censoring:** participants with severe dementia may drop out of follow-up
3. **Competing risks:** death competes with dementia incidence; consider Fine-Gray subdistribution hazards or cause-specific hazards
4. **Prevalent cases at baseline:** must be excluded (restrict to baseline normal cognition)
5. **Site heterogeneity:** random effects or meta-analysis to handle between-site variation
6. **Time scale:** age (preferred for incidence) vs. time-on-study
7. **Multiple imputation:** MICE for missing covariates, stratified by site

## Race/ethnicity — the multiracial innovation

- NACC has 1,340 multiracial participants (985 aged 65+)
- This is one of the largest samples of multiracial older adults in cognitive aging research
- 6-category race/ethnicity: White NH, Black NH, Hispanic, Asian NH, Multiracial NH, Other NH
- Can estimate race-specific incidence rates AND incidence rate ratios (vs. White NH reference)
- This addresses a major gap: virtually no published dementia incidence data for multiracial adults
- Frame as health equity contribution: multiracial population is the fastest-growing demographic group in the US

## Expected outputs

1. Generalized incidence rates (per 1,000 person-years) for MCI and dementia by race/ethnicity
2. Incidence rate ratios relative to White NH
3. Weighted survival curves by race/ethnicity
4. Heterogeneity assessment across ADRC sites (for Option B)
5. Comparison of generalized vs. unweighted incidence rates (quantifying selection bias)

## References

- Chan KCG, Xia F, Kukull WA (2025). NACC data: Who is represented over time and across centers, and implications for generalizability. Alzheimer's & Dementia, 21(9), e70657.
- Hayes-Larson E, et al. (2022). Accounting for lack of representation in dementia research. Alzheimer's & Dementia, 18(11), 2137-2146.
- Dahabreh IJ, et al. (2019). Generalizing causal inferences from randomized trials. Biometrics, 75(2), 685-694.
- Buchanan AL, et al. (2018). Generalizing evidence from randomized trials using inverse probability of sampling weights. JRSS-A, 181(4), 1193-1209.
- Rubin DB (1987). Multiple Imputation for Nonresponse in Surveys. Wiley.

## Existing codebase to adapt

The gen_sadrc_clin repository (/Users/lamhine/Documents/GitHub/gen_sadrc_clin) contains a complete single-site generalizability pipeline:
- 01_load_brfss.R: BRFSS loading (adapt to national)
- 02_clean_data.R: NACC UDS cleaning (adapt to multi-site)
- 03_harmonize.R: NACC <-> BRFSS covariate harmonization (reuse directly)
- 04_multiple_imputation.R: MICE imputation (adapt to include site)
- 05_weight_development.R: IOSW via svyglm (adapt to site-specific or pooled)
- 06_outcome_models.R: replace cross-sectional models with survival models
- 07_aipsw_analysis.R: adapt AIPSW to survival outcomes
- config.R: project configuration

Key adaptations needed:
- Cross-sectional outcome models -> Cox/parametric survival models
- Single-site weights -> multi-site (pooled or meta-analytic)
- 4-category race -> 6-category race
- 3-county BRFSS -> national BRFSS (or state-matched)
- Add competing risks framework
- Add meta-analysis combining step (for Option B)
