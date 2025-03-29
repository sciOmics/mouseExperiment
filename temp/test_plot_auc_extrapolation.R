# Test script to debug the extrapolated points issue in plot_auc
library(mouseExperiment)

# Load the synthetic dataset
data(combo_treatment_synthetic_data)

# Run tumor_auc_analysis to get AUC data
cat("Running tumor_auc_analysis to generate AUC data...\n")
auc_results <- tumor_auc_analysis(
  df = combo_treatment_synthetic_data,
  time_column = "Day",
  volume_column = "Volume", 
  treatment_column = "Treatment",
  id_column = "ID",
  cage_column = "Cage"
)

# Examine the Extrapolated column in the AUC data
cat("\nChecking Extrapolated column in AUC data:\n")
cat("- Column exists:", "Extrapolated" %in% colnames(auc_results$auc_data), "\n")
if ("Extrapolated" %in% colnames(auc_results$auc_data)) {
  cat("- TRUE values:", sum(auc_results$auc_data$Extrapolated), "\n")
  cat("- FALSE values:", sum(!auc_results$auc_data$Extrapolated), "\n")
  
  # Check by treatment group
  cat("\nExtrapolated points by treatment group:\n")
  for (treatment in unique(auc_results$auc_data$Treatment)) {
    treatment_data <- subset(auc_results$auc_data, Treatment == treatment)
    cat("  Treatment", treatment, ":\n")
    cat("    - Total mice:", nrow(treatment_data), "\n")
    cat("    - Extrapolated:", sum(treatment_data$Extrapolated), "\n")
    cat("    - Non-extrapolated:", sum(!treatment_data$Extrapolated), "\n")
  }
}

# Check the class of the Extrapolated column
cat("\nClass of Extrapolated column:", class(auc_results$auc_data$Extrapolated), "\n")

# Create a version with inverted extrapolation flags for comparison
cat("\nCreating a version with corrected extrapolation flags for testing...\n")
test_data <- auc_results$auc_data
if ("Extrapolated" %in% colnames(test_data)) {
  # Invert the flags for the HDACi + PD1 group as a test
  test_data$Extrapolated_Original <- test_data$Extrapolated
  test_data$Extrapolated <- ifelse(test_data$Treatment == "HDACi + PD1", 
                                   !test_data$Extrapolated, 
                                   test_data$Extrapolated)
}

# Create plot with original data
cat("\nCreating plot with original data...\n")
p1 <- plot_auc(
  auc_results$auc_data,
  title = "Original Data: Possible Incorrect Extrapolated Points",
  extrapolated_column = "Extrapolated"
)

# Create plot with test data
cat("\nCreating plot with test data (inverted flags for HDACi + PD1)...\n")
p2 <- plot_auc(
  test_data,
  title = "Test Data: Inverted Extrapolation Flags for HDACi + PD1",
  extrapolated_column = "Extrapolated"
)

# Output both plots for comparison
cat("\nSaving plots to PDF for comparison...\n")
pdf("temp/plot_auc_extrapolation_test.pdf", width = 12, height = 8)
print(p1)
print(p2)
dev.off()

cat("\nTest completed. Check 'temp/plot_auc_extrapolation_test.pdf' for visual comparison.\n")
cat("Compare the original data plot with the test plot to see if inverting the extrapolation\n")
cat("flags for the HDACi + PD1 group makes the visualization more accurate.\n") 