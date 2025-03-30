# Test script for post_power_analysis function
# This script tests the updated function to ensure it properly returns
# the three required objects: effect_sizes, post_power_analysis, and sample_size_estimates

library(mouseExperiment)

# Load demo data for a tumor growth experiment
data(combo_treatment_synthetic_data)

# Use synthetic data directly
processed_data <- combo_treatment_synthetic_data

cat("Running post_power_analysis with AUC method...\n")

# Perform post-hoc power analysis with AUC method
power_results_auc <- post_power_analysis(
  data = processed_data,
  alpha = c(0.05, 0.01),
  power = c(0.8, 0.9, 0.95),
  method = "auc"
)

# Check if the required output objects exist
cat("\nChecking if required output objects exist:\n")
cat("effect_sizes object exists: ", !is.null(power_results_auc$effect_sizes), "\n")
cat("post_power_analysis object exists: ", !is.null(power_results_auc$post_power_analysis), "\n")
cat("sample_size_estimates object exists: ", !is.null(power_results_auc$sample_size_estimates), "\n")

# Print the structure of each object
cat("\nStructure of effect_sizes object:\n")
str(power_results_auc$effect_sizes)

cat("\nStructure of post_power_analysis object:\n")
str(power_results_auc$post_power_analysis)

cat("\nStructure of sample_size_estimates object:\n")
str(power_results_auc$sample_size_estimates)

# Check if custom alpha and power levels are correctly used
cat("\nChecking if custom alpha levels are used in post_power_analysis:\n")
unique_alphas <- unique(power_results_auc$post_power_analysis$Alpha)
cat("Custom alpha levels found: ", paste(unique_alphas, collapse = ", "), "\n")
cat("Match expected values (0.05, 0.01): ", all(sort(unique_alphas) == sort(c(0.05, 0.01))), "\n")

cat("\nChecking if custom power levels are used in sample_size_estimates:\n")
unique_powers <- unique(power_results_auc$sample_size_estimates$Target_Power)
cat("Custom power levels found: ", paste(unique_powers, collapse = ", "), "\n")
cat("Match expected values (0.8, 0.9, 0.95): ", all(sort(unique_powers) == sort(c(0.8, 0.9, 0.95))), "\n")

cat("\nTest completed.\n") 