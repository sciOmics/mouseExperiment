# Test script for updated functions with composite identifiers
# This script demonstrates how to use the updated functions with the new data structure

# Load the package
library(mouseExperiment)

# Test with combo treatment data
cat("Testing with combo treatment data...\n")

# Read combo treatment data
data_path <- "/Users/david/Library/Mobile Documents/com~apple~CloudDocs/Programming/mouseExperiment/data/combo_treatment_synthetic_data.csv"
combo_data <- read.csv(data_path)

# Calculate tumor volumes
combo_data <- calculate_volume(combo_data)

# Calculate days from start date - the function will auto-detect the MM/DD/YYYY format
combo_data <- calculate_dates(combo_data, start_date = "03/24/2025")

# Use plot_tumor_growth with the new parameters
cat("Testing plot_tumor_growth with Treatment column...\n")
result <- try(plot_tumor_growth(combo_data, 
                             treatment_column = "Treatment", 
                             cage_column = "Cage", 
                             ID_column = "ID"))
if (inherits(result, "try-error")) {
  cat("ERROR: plot_tumor_growth function failed with combo data\n")
} else {
  cat("SUCCESS: plot_tumor_growth function worked with combo data\n")
}

# Use plot_survival with the new parameters
cat("\nTesting plot_survival with Treatment column...\n")
result <- try(plot_survival(combo_data, 
                         treatment_column = "Treatment", 
                         cage_column = "Cage", 
                         id_column = "ID"))
if (inherits(result, "try-error")) {
  cat("ERROR: plot_survival function failed with combo data\n")
} else {
  cat("SUCCESS: plot_survival function worked with combo data\n")
}

# Test with dose levels data
cat("\nTesting with dose levels data...\n")

# Read dose levels data
data_path <- "/Users/david/Library/Mobile Documents/com~apple~CloudDocs/Programming/mouseExperiment/data/dose_levels_synthetic_data.csv"
dose_data <- read.csv(data_path)

# Calculate tumor volumes
dose_data <- calculate_volume(dose_data)

# Calculate days from start date
dose_data <- calculate_dates(dose_data, 
                           start_date = "24-Mar", 
                           date_format = "%d-%b", 
                           year = 2023)

# Use plot_tumor_growth with the new parameters including dose
cat("Testing plot_tumor_growth with Treatment and Dose columns...\n")
result <- try(plot_tumor_growth(dose_data, 
                             treatment_column = "Treatment", 
                             cage_column = "Cage", 
                             ID_column = "ID", 
                             dose_column = "Dose"))
if (inherits(result, "try-error")) {
  cat("ERROR: plot_tumor_growth function failed with dose data\n")
} else {
  cat("SUCCESS: plot_tumor_growth function worked with dose data\n")
}

# Use plot_survival with the new parameters including dose
cat("\nTesting plot_survival with Treatment and Dose columns...\n")
result <- try(plot_survival(dose_data, 
                         treatment_column = "Treatment", 
                         cage_column = "Cage", 
                         id_column = "ID", 
                         dose_column = "Dose"))
if (inherits(result, "try-error")) {
  cat("ERROR: plot_survival function failed with dose data\n")
} else {
  cat("SUCCESS: plot_survival function worked with dose data\n")
}

cat("\nAll tests completed!\n")