# config.R
# Purpose: Define user-specific settings for file paths and directories

# Set data directories
adrc_data_dir <- "/Users/lamhine/Library/CloudStorage/Box-Box/Tracy\ Lam-Hine\'s\ Files/ADRC"
brfss_data_dir <- "/Users/lamhine/Library/CloudStorage/Box-Box/Tracy Lam-Hine's Files/brfss_disagg/CA-BRFSS/2014-2023"
processed_data_dir <- "/Users/lamhine/Documents/GitHub/gen_sadrc_clin/01_data"
results_dir <- "/Users/lamhine/Documents/GitHub/gen_sadrc_clin/03_results"


# Confirm directory setup
message("Using ADRC data directory: ", adrc_data_dir)
message("Using BRFSS data directory: ", brfss_data_dir)
message("Using processed data directory: ", processed_data_dir)
message("Using results directory: ", results_dir)

# Load commonly used packages
library(rstudioapi)
library(tidyverse)
library(summarytools)
library(janitor)
library(haven)      
library(data.table)  

# Suppress scientific notation for cleaner output
options(scipen = 999)

# ---- Analysis parameters ----

# Number of BRFSS years pooled (for weight adjustment per CDC guidance)
# When combining N years, divide each respondent's _LLCPWT by N so that
# weights represent the target population once, not N times.
N_BRFSS_YEARS <- 10  # 2015-2024 (aligned with ADRC enrollment Sept 2015 - Mar 2025)

# Set survey design settings
options(survey.lonely.psu = "adjust")