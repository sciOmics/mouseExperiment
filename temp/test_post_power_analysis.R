# Test script to verify the new Treatment column in post_power_analysis
library(mouseExperiment)

# Load the synthetic dataset
data(combo_treatment_synthetic_data)

# The synthetic data is already in the correct format with Volume and Day columns
# Use it directly without transformations
processed_data <- combo_treatment_synthetic_data

# Run the post_power_analysis function with AUC method
cat("Running post_power_analysis with AUC method...\n")
auc_power_results <- post_power_analysis(
  data = processed_data,
  method = "auc",
  effect_sizes = c(0.5, 0.8, 1.0),
  time_column = "Day",
  volume_column = "Volume",
  treatment_column = "Treatment",
  id_column = "ID"
)

# Check if power_analysis has the Treatment column
cat("\nChecking power_analysis result structure:\n")
cat("- Columns in power_analysis:", paste(colnames(auc_power_results$power_analysis), collapse=", "), "\n")
cat("- Does power_analysis have Treatment column:", "Treatment" %in% colnames(auc_power_results$power_analysis), "\n")

# Check unique treatments in the power_analysis
if ("Treatment" %in% colnames(auc_power_results$power_analysis)) {
  cat("- Unique treatments:", paste(unique(auc_power_results$power_analysis$Treatment), collapse=", "), "\n")
}

# Print the power_analysis data frame
cat("\nPower analysis results:\n")
print(auc_power_results$power_analysis)

# Check if the power curve plot includes treatment groups
cat("\nChecking if power curve plot uses treatment groups...\n")
plot_data <- ggplot2::ggplot_build(auc_power_results$plots$power_curve)$data[[1]]
cat("- Plot has multiple groups/colors:", length(unique(plot_data$group)) > 1, "\n")

# Save power curve plot to PDF
cat("\nSaving power curve plot to PDF...\n")
pdf("temp/power_analysis_treatment_plot.pdf", width = 10, height = 8)
print(auc_power_results$plots$power_curve)
dev.off()

cat("\nTest completed successfully. Check 'temp/power_analysis_treatment_plot.pdf' for visual verification.\n") 