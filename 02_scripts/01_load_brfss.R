# 01_load_brfss.R 
# Purpose: Load and merge California BRFSS microdata (2014-2023)

# ---------------------- #
# LOAD CONFIGURATION
# ---------------------- #
source("config.R")

# ---------------------- #
# LOAD & MERGE 2014-2023 CA BRFSS DATA
# ---------------------- #

# Get list of all `.sas7bdat` files in `brfss_data_dir`
brfss_files <- list.files(path = brfss_data_dir,
                          pattern = "\\.sas7bdat$",
                          full.names = TRUE)


# Read in BRFSS `.sas7bdat` files into a list
ca_brfss <- lapply(brfss_files, read_sas)

# Convert DATE column to date format for datasets where it is inconsistent
date_fix_indices <- c(1, 2, 5, 6, 7, 8)  # Specify datasets with mismatched DATE formats
for (i in date_fix_indices) {
  if ("DATE" %in% names(ca_brfss[[i]])) {
    ca_brfss[[i]]$DATE <- mdy(ca_brfss[[i]]$DATE)
  }
}

# Bind all list elements into one data frame
ca_bound <- rbindlist(ca_brfss, fill = TRUE) %>% as.data.frame()

# Ensure processed data directory exists before saving
if (!dir.exists(processed_data_dir)) {
  dir.create(processed_data_dir, recursive = TRUE)
  message("Created missing processed data directory: ",
          processed_data_dir)
}

# ---------------------- #
# SAVE FILES TO PROCESSED DATA DIRECTORY
# ---------------------- #

saveRDS(ca_bound, file = file.path(processed_data_dir, "01_ca_bound.rds"))
message("Saved raw BRFSS dataset: ",
        file.path(processed_data_dir, "01_ca_bound.rds"))

# Quick check of loaded data
str(ca_bound)

# End of script
message("01_load_data.R completed successfully.")