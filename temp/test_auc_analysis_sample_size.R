# Test script to verify that tumor_auc_analysis correctly identifies sample sizes
# with the combo_treatment_synthetic_data
library(mouseExperiment)

# Load the synthetic dataset
data(combo_treatment_synthetic_data)

# Process the data to calculate volume
processed_data <- combo_treatment_synthetic_data  # Volume is already calculated

# Summarize dataset to check our expectations
cat("Dataset Summary:\n")
cat("Total number of unique mice:", length(unique(processed_data$Mouse_ID)), "\n")
cat("Total number of unique ID-Treatment-Cage combinations:", 
    length(unique(paste(processed_data$ID, processed_data$Treatment, processed_data$Cage, sep="_"))), "\n")

# Count mice per treatment
cat("\nMice per treatment group:\n")
for (treatment in unique(processed_data$Treatment)) {
  treatment_data <- subset(processed_data, Treatment == treatment)
  treatment_mice <- unique(treatment_data$Mouse_ID)
  cat("  Treatment", treatment, "has", length(treatment_mice), "mice\n")
  
  # Count mice per cage within this treatment
  for (cage in unique(treatment_data$Cage)) {
    cage_data <- subset(treatment_data, Cage == cage)
    cage_mice <- unique(cage_data$Mouse_ID)
    cat("    Cage", cage, "has", length(cage_mice), "mice\n")
  }
}

# Run tumor_auc_analysis
cat("\nRunning tumor_auc_analysis...\n")
auc_results <- tumor_auc_analysis(
  df = processed_data,
  time_column = "Day",
  volume_column = "Volume", 
  treatment_column = "Treatment",
  id_column = "ID",
  cage_column = "Cage"
)

# Check the sample sizes in the results
cat("\nSample sizes in AUC summary:\n")
print(auc_results$auc_summary[, c("Treatment", "N")])

# Check if the total sample size matches our expectation
expected_total <- length(unique(paste(processed_data$ID, processed_data$Treatment, processed_data$Cage, sep="_")))
actual_total <- sum(auc_results$auc_summary$N)

cat("\nExpected total sample size:", expected_total, "\n")
cat("Actual total sample size in results:", actual_total, "\n")

# Check if we have the expected 8 mice per treatment
expected_per_treatment <- 8
cat("\nChecking if we have the expected", expected_per_treatment, "mice per treatment:\n")
for (treatment in unique(processed_data$Treatment)) {
  actual_n <- auc_results$auc_summary$N[auc_results$auc_summary$Treatment == treatment]
  if (actual_n == expected_per_treatment) {
    cat("  SUCCESS: Treatment", treatment, "has the expected", expected_per_treatment, "mice\n")
  } else {
    cat("  FAILURE: Treatment", treatment, "has", actual_n, "mice instead of the expected", expected_per_treatment, "\n")
  }
}

# Examine the individual AUC data to verify each mouse is correctly identified
cat("\nVerifying individual mice identification:\n")
auc_data <- auc_results$auc_data

# Count mice per treatment in the AUC data
for (treatment in unique(auc_data$Treatment)) {
  treatment_auc_data <- subset(auc_data, Treatment == treatment)
  cat("  Treatment", treatment, "has", nrow(treatment_auc_data), "mice in AUC results\n")
  
  # Check if we have records for both cages
  if ("Cage" %in% colnames(auc_data)) {
    cages_in_treatment <- unique(treatment_auc_data$Cage)
    cat("    Cages found:", paste(cages_in_treatment, collapse=", "), "\n")
    
    for (cage in cages_in_treatment) {
      cage_auc_data <- subset(treatment_auc_data, Cage == cage)
      cat("      Cage", cage, "has", nrow(cage_auc_data), "mice in AUC results\n")
    }
  }
}

# Final conclusion
if (actual_total == expected_total) {
  cat("\nSUCCESS: tumor_auc_analysis correctly identified all", expected_total, "mice\n")
} else {
  cat("\nFAILURE: tumor_auc_analysis identified", actual_total, "mice instead of the expected", expected_total, "\n")
} 