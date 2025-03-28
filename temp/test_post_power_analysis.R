# Test script for post_power_analysis function with AUC method
library(mouseExperiment)
library(ggplot2)

# Set seed for reproducibility
set.seed(123)

# Create a synthetic dataset with 2 treatment groups
create_test_data <- function(n_per_group = 10, noise_level = 0.1) {
  treatment_groups <- c("Control", "Treatment")
  time_points <- seq(0, 21, by = 3)
  
  # Define growth rates for each group
  control_rate <- 0.15   # Higher growth rate (worse)
  treatment_rate <- 0.08 # Lower growth rate (better)
  
  # Create an empty data frame
  all_data <- data.frame()
  
  # Generate data for each treatment group
  for (treatment_idx in 1:length(treatment_groups)) {
    treatment <- treatment_groups[treatment_idx]
    growth_rate <- if (treatment == "Control") control_rate else treatment_rate
    
    for (subject_idx in 1:n_per_group) {
      subject_id <- paste0(substr(treatment, 1, 1), subject_idx)
      
      # Initial volume with some randomness
      initial_volume <- 50 + rnorm(1, 0, 10)
      
      # Generate data for each time point
      for (time_idx in 1:length(time_points)) {
        time <- time_points[time_idx]
        
        # Exponential growth model with noise
        volume <- initial_volume * exp(growth_rate * time) * (1 + rnorm(1, 0, noise_level))
        
        # Add to the data frame
        all_data <- rbind(all_data, data.frame(
          ID = subject_id,
          Treatment = treatment,
          Day = time,
          Volume = volume,
          Cage = paste0("Cage", ceiling(subject_idx/2))
        ))
      }
    }
  }
  
  return(all_data)
}

# Generate the test data
cat("Creating synthetic tumor growth data...\n")
tumor_data <- create_test_data(n_per_group = 10)

# Print sample of data
cat("\nSample of synthetic data:\n")
print(head(tumor_data))

# Run post_power_analysis with AUC method
cat("\nRunning post_power_analysis with AUC method...\n")
tryCatch({
  power_results <- post_power_analysis(
    data = tumor_data,
    method = "auc",
    effect_sizes = c(0.5, 0.8, 1.0, 1.5, 2.0)
  )
  
  # Print power analysis results
  cat("\nPower analysis results:\n")
  print(power_results$power_analysis)
  
  # Print sample size recommendations
  cat("\nSample size recommendations:\n")
  print(power_results$sample_size_recommendations)
  
  # Save power curve plot if available
  if (!is.null(power_results$plots$power_curve)) {
    pdf("temp/power_curve_plot.pdf", width = 8, height = 6)
    print(power_results$plots$power_curve)
    dev.off()
    cat("\nPower curve plot saved to temp/power_curve_plot.pdf\n")
  }
  
  cat("\nTest PASSED: post_power_analysis with AUC method completed successfully.\n")
}, error = function(e) {
  cat("\nTest FAILED: Error in post_power_analysis -", e$message, "\n")
}) 