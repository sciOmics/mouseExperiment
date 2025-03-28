# Script to analyze why 2 mice from the aPD1 group are missing in AUC results
library(mouseExperiment)

# Load data
data(combo_treatment_synthetic_data)

# Use both aPD1 and Control treatment groups to avoid contrast error
test_data <- subset(combo_treatment_synthetic_data, Treatment %in% c("aPD1", "Control"))

# Check for NaN values in aPD1 data
aPD1_data <- subset(test_data, Treatment == "aPD1")
cat("Total NaN values in aPD1 data:", sum(is.nan(aPD1_data$Volume)), "\n\n")

# Check each mouse ID in aPD1 group
cat("NaN values by Mouse ID in aPD1 group:\n")
for (id in unique(aPD1_data$ID)) {
  mouse_data <- aPD1_data[aPD1_data$ID == id, ]
  nan_count <- sum(is.nan(mouse_data$Volume))
  total_count <- nrow(mouse_data)
  
  cat("ID:", id, "- NaN values:", nan_count, "/", total_count, 
      "(", round(100 * nan_count / total_count, 1), "%)\n")
}

# Run tumor_auc_analysis and check which mice are included
cat("\nRunning tumor_auc_analysis on aPD1 and Control data...\n")
auc_results <- tumor_auc_analysis(
  df = test_data,
  time_column = "Day",
  volume_column = "Volume", 
  treatment_column = "Treatment",
  id_column = "ID",
  cage_column = "Cage"
)

# Check the sample sizes in the results
cat("\nSample sizes in AUC summary:\n")
print(auc_results$auc_summary[, c("Treatment", "N")])

# Check which aPD1 mice are included in the results
aPD1_auc_data <- subset(auc_results$auc_data, Treatment == "aPD1")
cat("\nMice included in aPD1 AUC results:\n")
for (id in unique(aPD1_data$ID)) {
  if (id %in% aPD1_auc_data$Subject) {
    cat("ID:", id, "is included\n")
  } else {
    cat("ID:", id, "is NOT included\n")
  }
}

# Print how many time points are required for AUC calculation
# (Determined by examining the tumor_auc_analysis.R code)
cat("\nTime points needed for AUC calculation:\n")
cat("The tumor_auc_analysis function requires at least 2 valid time points to calculate AUC\n")

# Check how many non-NaN measurements each mouse has
cat("\nNumber of valid measurements by Mouse ID in aPD1 group:\n")
for (id in unique(aPD1_data$ID)) {
  mouse_data <- aPD1_data[aPD1_data$ID == id, ]
  valid_count <- sum(!is.nan(mouse_data$Volume))
  total_count <- nrow(mouse_data)
  
  cat("ID:", id, "- Valid measurements:", valid_count, "/", total_count, "\n")
  
  # Print the first valid measurement
  if (valid_count > 0) {
    first_valid <- mouse_data[!is.nan(mouse_data$Volume), ][1, ]
    cat("  First valid measurement: Day", first_valid$Day, "Volume", first_valid$Volume, "\n")
  }
}

# Print the actual AUC data for aPD1 group
cat("\nAUC data for aPD1 group:\n")
print(aPD1_auc_data) 