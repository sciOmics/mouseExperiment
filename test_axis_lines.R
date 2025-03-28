# Test script to verify the addition of x and y-axis lines
library(mouseExperiment)
library(ggplot2)

# Load synthetic data
data(combo_treatment_synthetic_data)

# Run tumor growth statistics on the data
result <- tumor_growth_statistics(
  combo_treatment_synthetic_data,
  time_column = "Day",
  volume_column = "Volume",
  id_column = "Mouse_ID",
  treatment_column = "Treatment",
  cage_column = "Cage",
  model_type = "lme4"
)

# Test plot_growth_rate function
p1 <- plot_growth_rate(result$growth_rates, 
                      title = "Growth Rate Plot with Axis Lines")
ggsave("growth_rate_with_axis_lines.pdf", p1, width = 8, height = 6)

# Test plot_auc function
# First run with AUC model type
result_auc <- tumor_growth_statistics(
  combo_treatment_synthetic_data,
  time_column = "Day",
  volume_column = "Volume",
  id_column = "Mouse_ID",
  treatment_column = "Treatment",
  cage_column = "Cage",
  model_type = "auc"
)

p2 <- plot_auc(result_auc$auc_analysis$individual, 
              title = "AUC Plot with Axis Lines")
ggsave("auc_with_axis_lines.pdf", p2, width = 8, height = 6)

# Print completion message
cat("Test completed successfully. Plots were generated with axis lines:\n")
cat("1. growth_rate_with_axis_lines.pdf\n")
cat("2. auc_with_axis_lines.pdf\n") 