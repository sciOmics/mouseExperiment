# Example script for testing dose-response relationship
# This is a complete example that shows how to test if there is a
# dose-response relationship in mouse tumor data

# Load the package
library(mouseExperiment)

# Read dose-response data
data_path <- "/Users/david/Library/Mobile Documents/com~apple~CloudDocs/Programming/mouseExperiment/data/dose_levels_synthetic_data.csv"
dose_data <- read.csv(data_path)

# Calculate tumor volumes (required for analysis)
dose_data <- calculate_volume(dose_data)

# Calculate days from start date
dose_data <- calculate_dates(dose_data, 
                            start_date = "24-Mar", 
                            date_format = "%d-%b", 
                            year = 2023)

# Test for dose-response relationship at the last time point
# (this analyzes the final tumor measurements)
results <- dose_response_statistics(dose_data,
                                  dose_column = "Dose",
                                  treatment_column = "Treatment",
                                  volume_column = "Volume")

# The linear regression p-value directly tests if there's a 
# linear relationship between dose and tumor volume
# A p-value < 0.05 indicates a significant dose-response
cat("\n=================================================\n")
cat("Test for linear dose-response relationship:\n")
cat("p-value =", format.pval(results$dose_effect_test$linear_pvalue, digits = 3), "\n")
cat("Slope =", round(results$dose_effect_test$slope, 4), "\n")

if (results$dose_effect_test$linear_pvalue < 0.05) {
  cat("RESULT: Significant dose-response relationship detected.\n")
  if (results$dose_effect_test$slope < 0) {
    cat("As dose increases, tumor volume decreases (treatment is effective).\n")
  } else {
    cat("As dose increases, tumor volume increases (dose may be harmful).\n")
  }
} else {
  cat("RESULT: No significant linear dose-response relationship detected.\n")
}
cat("=================================================\n")

# Optionally analyze at a specific time point
# For example, day 14 (or the closest day to 14 in the data)
days <- unique(dose_data$Day)
day14_approx <- days[which.min(abs(days - 14))]

day14_results <- dose_response_statistics(dose_data,
                                        dose_column = "Dose",
                                        treatment_column = "Treatment",
                                        volume_column = "Volume",
                                        time_point = day14_approx)

cat("\n=================================================\n")
cat("Test for dose-response relationship at day", day14_approx, ":\n")
cat("p-value =", format.pval(day14_results$dose_effect_test$linear_pvalue, digits = 3), "\n")
cat("Slope =", round(day14_results$dose_effect_test$slope, 4), "\n")

if (day14_results$dose_effect_test$linear_pvalue < 0.05) {
  cat("RESULT: Significant dose-response relationship detected at day", day14_approx, ".\n")
  if (day14_results$dose_effect_test$slope < 0) {
    cat("As dose increases, tumor volume decreases (treatment is effective).\n")
  } else {
    cat("As dose increases, tumor volume increases (dose may be harmful).\n")
  }
} else {
  cat("RESULT: No significant linear dose-response relationship at day", day14_approx, ".\n")
}
cat("=================================================\n")

# Plot the dose-response relationship
plot_day <- day14_approx # Choose which day to plot
plot_data <- dose_data[dose_data$Day == plot_day, ]

# Simple scatter plot with regression line
par(mfrow = c(1, 1), mar = c(5, 5, 4, 2))
plot(plot_data$Dose, plot_data$Volume, 
     xlab = "Dose", ylab = "Tumor Volume",
     main = paste("Dose-Response on Day", plot_day),
     pch = 16, cex = 1.5)

# Add regression line
abline(lm(Volume ~ Dose, data = plot_data), col = "red", lwd = 2)

# Add text with p-value
p_val <- format.pval(day14_results$dose_effect_test$linear_pvalue, digits = 3)
legend("topright", legend = paste("p =", p_val), bty = "n")