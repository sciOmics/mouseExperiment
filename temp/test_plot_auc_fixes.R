# Test script to verify the fixes in the plot_auc function
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

# Create a plot using the default parameters
cat("\nCreating plot with default parameters...\n")
cat("- show_mean should be TRUE by default\n")
cat("- Mean should be displayed as a horizontal bar (not diamond)\n")
cat("- Extrapolated points should be correctly marked\n")
cat("- Mean bars should be in the correct group columns\n")
p1 <- plot_auc(
  auc_results$auc_data,
  title = "Default Plot: Should show horizontal mean bars"
)

# Create a plot with error bars
cat("\nCreating plot with error bars...\n")
p2 <- plot_auc(
  auc_results$auc_data,
  error_bar_type = "SEM",
  title = "Plot with SEM error bars"
)

# Display the plots
cat("\nPrinting plots to PDF...\n")
pdf("temp/plot_auc_test.pdf", width = 12, height = 8)
print(p1)
print(p2)
dev.off()

cat("\nTest completed. Plot saved to temp/plot_auc_test.pdf\n")
cat("Please check the visualization for:\n")
cat("1. Horizontal mean bars (not diamonds)\n")
cat("2. Mean bars in the correct columns\n")
cat("3. Proper handling of extrapolated points\n")
cat("4. Consistent display of all mice per treatment group\n") 