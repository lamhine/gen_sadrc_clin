# Generalizability of the Stanford ADRC Clinical Core

## Overview

This repository contains analysis code for a generalizability assessment of the Stanford Alzheimer's Disease Research Center (ADRC) Clinical Core relative to the surrounding three-county population (Alameda, San Mateo, and Santa Clara counties, California). The analysis combines ADRC participant data (NACC UDS v3) with population-representative California BRFSS microdata (2015-2024) to address three aims:

1. **Local AD prevalence estimation** — A hybrid model combining published regression coefficients from the Chicago Health and Aging Project (Dhana et al. 2023) with locally-derived race/ethnicity coefficients from the Stanford ADRC, applied to individual-level BRFSS data.

2. **Cognitive outcome disparity analysis** — Augmented Inverse Probability of Selection Weighting (AIPSW), a doubly-robust generalizability estimator, to estimate demographically-standardized prevalence ratios for cognitive outcomes across racial/ethnic groups.

3. **Clinical enrichment characterization** — Quantification of outcome-dependent selection in the ADRC using referral source data and external prevalence benchmarks, as a transparency diagnostic.

## Repository Structure

```
gen_sadrc_clin/
├── config.R                  # User-specific paths and analysis parameters
├── gen_sadrc_clin.Rproj      # RStudio project file
├── 01_data/                  # Processed data (gitignored)
├── 02_scripts/               # Analysis pipeline (see below)
├── 03_results/               # Output tables, figures, and model objects
└── 04_docs/                  # Methods draft, decision log, poster
```

## Analysis Pipeline

Scripts are numbered and should be run sequentially from the repository root:

| Script | Purpose |
|--------|---------|
| `01_load_brfss.R` | Load and merge California BRFSS microdata (2015-2024) |
| `02_clean_data.R` | Clean ADRC (NACC UDS v3) and BRFSS data |
| `03_harmonize.R` | Harmonize variables across ADRC and BRFSS into a combined dataset |
| `04_multiple_imputation.R` | Multiple imputation via MICE for missing covariates |
| `05_weight_development.R` | Develop stabilized inverse odds of selection weights (sIOSW) |
| `06_outcome_models.R` | Fit outcome models in ADRC sample (outcome component of AIPSW) |
| `07_aipsw_analysis.R` | Primary AIPSW analysis with 1,000-rep bootstrap inference |
| `08_aim1_representativeness.R` | Covariate balance diagnostics (love plots, positivity) |
| `09_results_summary.R` | Publication-ready figures and tables |
| `10_sensitivity_htn.R` | Sensitivity analysis: hypertension in participation model |
| `10_sensitivity_trimmed_weights.R` | Sensitivity analysis: 99th-percentile weight trimming |
| `12_hybrid_prevalence.R` | Hybrid AD prevalence estimation (Aim 1) |
| `13_stakeholder_summary.R` | Combined summary tables across all three aims |

## Data Availability

- **ADRC data**: Stanford ADRC Clinical Core data (NACC UDS v3) are available through the [National Alzheimer's Coordinating Center](https://naccdata.org/) under a data use agreement.
- **BRFSS data**: California BRFSS microdata are publicly available from the [CDC BRFSS website](https://www.cdc.gov/brfss/annual_data/annual_data.htm). This analysis uses the 2015-2024 survey years, restricted to the three-county catchment area (FIPS: 001, 081, 085) and adults aged 65+.

Individual-level data files are not included in this repository (`01_data/` is gitignored).

## Setup

1. Clone the repository and open `gen_sadrc_clin.Rproj` in RStudio.

2. Edit `config.R` to set local paths:
   - `adrc_data_dir` — Path to ADRC data extract (Box or local)
   - `brfss_data_dir` — Path to downloaded BRFSS files
   - `processed_data_dir` — Path for intermediate data files
   - `results_dir` — Path for output tables and figures

3. Run scripts sequentially from the repo root:
   ```r
   source("02_scripts/01_load_brfss.R")
   source("02_scripts/02_clean_data.R")
   # ... etc.
   ```

## Key Dependencies

R 4.5.0 with: `mice`, `survey`, `boot`, `tidyverse`, `ggplot2`, `posterdown` (for poster output)

## Authors

Tracy Lam-Hine, Lisa Goldman-Rosas, Victor W. Henderson  
Stanford University School of Medicine
