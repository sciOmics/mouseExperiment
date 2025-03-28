library(mouseExperiment)

# Create a synthetic dataset with unequal variances between groups
set.seed(123)

# Create a dataset with unequal variances between treatment groups
n_per_group <- 8
treatments <- c("Control", "Drug A", "Drug B", "Combination")

# Create different time points
times <- c(0, 3, 7, 10, 14, 17, 21, 24, 28)

# Function to generate volume data with different growth patterns and variances
generate_volume <- function(treatment, time_points, n_subjects) {
  subjects <- paste0("M", sprintf("%02d", 1:n_subjects))
  
  # Different growth rates and variability for each treatment
  if (treatment == "Control") {
    # Fast growth with medium variance
    base_size <- 0.15
    growth_rate <- 0.15
    noise_sd <- 0.1
  } else if (treatment == "Drug A") {
    # Medium growth with high variance
    base_size <- 0.15
    growth_rate <- 0.1
    noise_sd <- 0.2  # Higher variability
  } else if (treatment == "Drug B") {
    # Medium growth with low variance
    base_size <- 0.15
    growth_rate <- 0.08
    noise_sd <- 0.05  # Lower variability
  } else if (treatment == "Combination") {
    # Slow growth with medium variance
    base_size <- 0.15
    growth_rate <- 0.05
    noise_sd <- 0.1
  }
  
  # Generate data for all subjects in this treatment group
  result <- data.frame()
  
  for (subject in subjects) {
    # Individual variation in growth rate
    subject_growth <- growth_rate * (1 + rnorm(1, 0, 0.2))
    
    # Generate volumes for each time point
    for (t in time_points) {
      # Calculate expected volume with exponential growth
      expected_volume <- base_size * exp(subject_growth * t)
      
      # Add noise
      actual_volume <- expected_volume * (1 + rnorm(1, 0, noise_sd))
      
      # Ensure volume is positive
      actual_volume <- max(0.01, actual_volume)
      
      # Add to result
      result <- rbind(result, data.frame(
        Mouse_ID = subject,
        Day = t,
        Treatment = treatment,
        Volume = actual_volume,
        ID = subject,
        Cage = paste0("C", ceiling(which(subjects == subject) / 2))
      ))
    }
  }
  
  return(result)
}

# Generate data for each treatment
control_data <- generate_volume("Control", times, n_per_group)
drugA_data <- generate_volume("Drug A", times, n_per_group)
drugB_data <- generate_volume("Drug B", times, n_per_group)
combo_data <- generate_volume("Combination", times, n_per_group)

# Combine all data
all_data <- rbind(control_data, drugA_data, drugB_data, combo_data)

# Test tumor_growth_statistics with AUC model
cat("Testing tumor_growth_statistics with AUC model...\n")
results <- tumor_growth_statistics(
  df = all_data,
  time_column = "Day",
  volume_column = "Volume",
  treatment_column = "Treatment",
  id_column = "Mouse_ID",
  cage_column = "Cage",
  model_type = "auc"
)

# Print posthoc results to verify Welch's t-tests are being used
cat("\n===== Posthoc Test Results =====\n")
if (!is.null(results$posthoc) && is.list(results$posthoc)) {
  cat("Method:", results$posthoc$method, "\n\n")
  print(results$posthoc$pairwise)
} else {
  cat("Posthoc tests not available\n")
}

# Print the analysis summary section showing the methods used
cat("\n===== Analysis Methods =====\n")
if (!is.null(results$summary) && !is.null(results$summary$methods)) {
  cat("Statistical test:", results$summary$methods$statistical_test, "\n")
  cat("Posthoc method:", results$summary$methods$posthoc_method, "\n")
}

# Check the variances in each group to show why Welch's test is appropriate
cat("\n===== Group Variances =====\n")
auc_data <- results$auc_analysis$individual
variances <- tapply(auc_data$AUC, auc_data$Treatment, var)
for (treatment in names(variances)) {
  cat(treatment, "variance:", variances[treatment], "\n")
}

# Create variance ratio to demonstrate heterogeneity
max_var <- max(variances)
min_var <- min(variances)
variance_ratio <- max_var / min_var
cat("\nMax/Min variance ratio:", variance_ratio, "\n")
cat("(Ratio > 4 typically indicates substantial heterogeneity of variance)\n")

# Optional: create a PDF with a diagnostic plot
pdf("temp/auc_values_by_group.pdf", width = 10, height = 6)
boxplot(AUC ~ Treatment, data = results$auc_analysis$individual, 
        main = "AUC Values by Treatment Group",
        ylab = "Area Under the Curve")
stripchart(AUC ~ Treatment, data = results$auc_analysis$individual,
          method = "jitter", vertical = TRUE, add = TRUE,
          pch = 16, col = "red", alpha = 0.5)
dev.off()

cat("\nTest completed. Results saved to temp/auc_values_by_group.pdf\n") 