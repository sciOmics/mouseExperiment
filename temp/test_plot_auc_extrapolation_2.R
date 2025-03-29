# Test script to verify the fix for the extrapolated points issue in plot_auc
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

# Create a version with manually set extrapolation flags for testing
# For demonstration purposes, we'll mark half of each treatment group as extrapolated
cat("\nCreating a version with manually set extrapolation flags for testing...\n")
test_data <- auc_results$auc_data

# Set half of each treatment group to be extrapolated
for (treatment in unique(test_data$Treatment)) {
  treatment_rows <- which(test_data$Treatment == treatment)
  n_rows <- length(treatment_rows)
  
  # Mark approximately half as extrapolated
  half_point <- ceiling(n_rows / 2)
  extrapolated_rows <- treatment_rows[1:half_point]
  
  # Set extrapolation flags
  test_data$Extrapolated[extrapolated_rows] <- TRUE
}

# Verify our modification
cat("\nVerifying manual modifications to extrapolation flags:\n")
for (treatment in unique(test_data$Treatment)) {
  treatment_data <- subset(test_data, Treatment == treatment)
  cat("  Treatment", treatment, ":\n")
  cat("    - Total mice:", nrow(treatment_data), "\n")
  cat("    - Extrapolated:", sum(treatment_data$Extrapolated), "\n")
  cat("    - Non-extrapolated:", sum(!treatment_data$Extrapolated), "\n")
}

# Create plot with test data
cat("\nCreating plot with test data (manually set extrapolation flags)...\n")
p <- plot_auc(
  test_data,
  title = "Test Data: Manually Set Extrapolation Flags",
  extrapolated_column = "Extrapolated"
)

# Save the plot to a PDF file for visual inspection
cat("\nSaving plot to PDF for visual inspection...\n")
pdf("temp/plot_auc_extrapolation_fixed.pdf", width = 10, height = 8)
print(p)
dev.off()

cat("\nTest completed. Check 'temp/plot_auc_extrapolation_fixed.pdf' for visual verification.\n")
cat("The plot should now correctly show approximately half of each treatment group as\n")
cat("extrapolated (open circles) and half as non-extrapolated (filled circles).\n") 