library(mouseExperiment)

# Create a test dataset for survival statistics
set.seed(123)
treatments <- c("Control", "PD1", "HDACi", "HDACi + PD1")

# Generate survival data with HDACi + PD1 having intentionally problematic values
# to test the robustness of our forest plot function
test_data <- data.frame(
  Group = treatments,
  HR = c(1.0, 0.7, 0.5, NA),  # NA value for HDACi + PD1
  CI_Lower = c(1.0, 0.4, 0.3, NA),
  CI_Upper = c(1.0, 1.2, 0.8, NA),
  P_Value = c(NA, 0.04, 0.01, 0.001),
  Note = c("Reference group", "", "", ""),
  Events = c(5, 4, 3, 2),
  Total = c(8, 8, 8, 8)
)

# Test the plot_forest function directly
cat("Testing plot_forest function from plot_forest.R with problematic data:\n")
pdf("temp/fixed_forest_plot.pdf", width = 10, height = 6)
plot <- plot_forest(test_data)
print(plot)
dev.off()

cat("Test completed successfully.\n") 