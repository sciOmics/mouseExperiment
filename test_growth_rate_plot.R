# Test script for plot_growth_rate function
library(mouseExperiment)
library(ggplot2)

# Load the synthetic data
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

# Look at the structure of the growth_rates
str(result$growth_rates)

# Get unique treatment groups for reference
unique_treatments <- unique(result$growth_rates$Treatment)
print(unique_treatments)

# Create basic growth rate plot
p1 <- plot_growth_rate(result$growth_rates)
ggsave("basic_growth_rate_plot.pdf", p1, width = 8, height = 6)

# Create growth rate plot with SD error bars
p2 <- plot_growth_rate(result$growth_rates, error_bar_type = "SD")
ggsave("sd_growth_rate_plot.pdf", p2, width = 8, height = 6)

# Create growth rate plot with custom colors
colors <- c("Control" = "gray", "HDACi" = "blue", 
           "aPD1" = "green", "HDACi + PD1" = "red")
p3 <- plot_growth_rate(result$growth_rates, colors = colors, title = "Custom Growth Rate Plot")
ggsave("custom_growth_rate_plot.pdf", p3, width = 8, height = 6)

# Create growth rate plot with custom order
p4 <- plot_growth_rate(result$growth_rates, 
                      group_order = c("Control", "aPD1", "HDACi", "HDACi + PD1"),
                      show_mean = TRUE,
                      error_bar_type = "CI")
ggsave("ordered_growth_rate_plot.pdf", p4, width = 8, height = 6)

# Print completion message
cat("Test completed successfully. Four plots were generated:\n")
cat("1. basic_growth_rate_plot.pdf\n")
cat("2. sd_growth_rate_plot.pdf\n")
cat("3. custom_growth_rate_plot.pdf\n")
cat("4. ordered_growth_rate_plot.pdf\n") 