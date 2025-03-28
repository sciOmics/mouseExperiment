# Test script to verify that tumor_growth_statistics correctly reports sample sizes
# with the combo_treatment_synthetic_data
library(mouseExperiment)

# Load the synthetic dataset
data(combo_treatment_synthetic_data)

# Run tumor_growth_statistics with AUC model type
cat("Running tumor_growth_statistics with AUC model...\n")
results_auc <- tumor_growth_statistics(
  df = combo_treatment_synthetic_data,
  time_column = "Day",
  volume_column = "Volume", 
  treatment_column = "Treatment",
  id_column = "ID",
  cage_column = "Cage",
  model_type = "auc"
)

# Run tumor_growth_statistics with lme4 model type
cat("\nRunning tumor_growth_statistics with lme4 model...\n")
results_lme <- tumor_growth_statistics(
  df = combo_treatment_synthetic_data,
  time_column = "Day",
  volume_column = "Volume", 
  treatment_column = "Treatment",
  id_column = "ID",
  cage_column = "Cage",
  model_type = "lme4"
)

# Check the reported sample sizes
cat("\nAUC model data_description:\n")
print(results_auc$summary$data_description)

cat("\nlme4 model data_description:\n")
print(results_lme$summary$data_description)

# Calculate the correct sample size directly from the data
correct_sample_size <- length(unique(paste(
  combo_treatment_synthetic_data$ID, 
  combo_treatment_synthetic_data$Treatment, 
  combo_treatment_synthetic_data$Cage, 
  sep="_"
)))

cat("\nIndependent calculation of total unique mice:", correct_sample_size, "\n")

# Count mice per treatment group
cat("\nMice per treatment group:\n")
for (treatment in unique(combo_treatment_synthetic_data$Treatment)) {
  treatment_data <- subset(combo_treatment_synthetic_data, Treatment == treatment)
  unique_mice <- length(unique(paste(
    treatment_data$ID, 
    treatment_data$Treatment, 
    treatment_data$Cage, 
    sep="_"
  )))
  cat("  Treatment", treatment, "has", unique_mice, "mice\n")
}

# Verify each model's reported sample size against the correct value
cat("\nVerifying reported sample sizes:\n")

if (results_auc$summary$data_description$subjects == correct_sample_size) {
  cat("SUCCESS: AUC model correctly reports", correct_sample_size, "subjects\n")
} else {
  cat("FAILURE: AUC model reports", results_auc$summary$data_description$subjects, 
      "subjects instead of the expected", correct_sample_size, "\n")
}

if (results_lme$summary$data_description$subjects == correct_sample_size) {
  cat("SUCCESS: lme4 model correctly reports", correct_sample_size, "subjects\n")
} else {
  cat("FAILURE: lme4 model reports", results_lme$summary$data_description$subjects, 
      "subjects instead of the expected", correct_sample_size, "\n")
} 