# Test script for extrapolation_points parameter in plot_auc and tumor_auc_analysis
library(mouseExperiment)
library(ggplot2)

# Set seed for reproducibility
set.seed(123)

# Create synthetic data with some mice having incomplete data
create_synthetic_tumor_data <- function(n_per_group = 10) {
  treatments <- c("Control", "Drug A", "Drug B", "Combination")
  time_points <- seq(0, 21, by = 3)
  
  all_data <- data.frame()
  
  for (treatment in treatments) {
    for (id in 1:(n_per_group)) {
      # Determine if this mouse will have incomplete data
      # More likely in control group, less likely in treatment groups
      incomplete <- FALSE
      dropout_time <- max(time_points)
      
      if (treatment == "Control" && runif(1) < 0.4) {
        incomplete <- TRUE
        # Drop out at day 9, 12, or 15
        dropout_time <- sample(c(9, 12, 15), 1)
      } else if (treatment != "Combination" && runif(1) < 0.2) {
        incomplete <- TRUE
        # Drop out at day 15 or 18
        dropout_time <- sample(c(15, 18), 1)
      }
      
      # Generate data points for this mouse up to dropout time
      mouse_data <- data.frame()
      
      for (day in time_points) {
        if (day <= dropout_time) {
          # Generate volume based on treatment
          # Control: rapid growth
          # Drug A: moderate growth
          # Drug B: slow growth
          # Combination: minimal growth
          base_volume <- 100 + day * 10  # Base volume increases with time
          
          growth_factor <- switch(treatment,
                                  "Control" = 5,
                                  "Drug A" = 3,
                                  "Drug B" = 2,
                                  "Combination" = 1)
          
          # Add randomness
          noise <- rnorm(1, mean = 0, sd = 20)
          
          volume <- base_volume + (day * growth_factor * 10) + noise
          volume <- max(volume, 50)  # Minimum tumor size
          
          mouse_data <- rbind(mouse_data, data.frame(
            ID = paste0(treatment, "_", id),
            Treatment = treatment,
            Day = day,
            Volume = volume
          ))
        }
      }
      
      all_data <- rbind(all_data, mouse_data)
    }
  }
  
  return(all_data)
}

# Generate synthetic data
tumor_data <- create_synthetic_tumor_data()

# Print summary of data points per mouse to verify we have incomplete data
data_points_summary <- aggregate(Day ~ ID + Treatment, data = tumor_data, FUN = function(x) length(x))
names(data_points_summary)[3] <- "NumDataPoints"
print("Summary of data points per mouse:")
print(table(data_points_summary$Treatment, data_points_summary$NumDataPoints))

# Test with default extrapolation_points (3)
cat("\n\nRunning tumor_auc_analysis with default extrapolation_points (3)...\n")
auc_results_default <- tumor_auc_analysis(tumor_data)

# Examine extrapolation information
cat("\nExtrapolation summary with default extrapolation_points (3):\n")
print(auc_results_default$summary[, c("Treatment", "N_Extrapolated", "N_Total", "Pct_Extrapolated")])

# Save the default plot
pdf("temp/auc_default_extrapolation.pdf", width = 10, height = 8)
print(auc_results_default$plot)
dev.off()

# Test with higher extrapolation_points (5)
cat("\n\nRunning tumor_auc_analysis with extrapolation_points = 5...\n")
auc_results_strict <- tumor_auc_analysis(tumor_data, extrapolation_points = 5)

# Examine extrapolation information
cat("\nExtrapolation summary with stricter extrapolation_points (5):\n")
print(auc_results_strict$summary[, c("Treatment", "N_Extrapolated", "N_Total", "Pct_Extrapolated")])

# Save the stricter plot
pdf("temp/auc_strict_extrapolation.pdf", width = 10, height = 8)
print(auc_results_strict$plot)
dev.off()

# Direct test of plot_auc function
cat("\n\nTesting plot_auc function directly...\n")

# Get the individual AUC data from previous analysis and ensure it has the correct structure
auc_data <- auc_results_default$auc_data

# Make sure the data frame has the required structure for plot_auc
if (!"Group" %in% colnames(auc_data)) {
  auc_data$Group <- auc_data$Treatment
}

cat("Data structure for plot_auc:\n")
str(auc_data)

# Create plots with different extrapolation_points settings
p1 <- plot_auc(
  auc_data = auc_data,
  title = "AUC with extrapolation_points = 3",
  extrapolation_points = 3,
  extrapolated_column = "Extrapolated",
  show_mean = TRUE
)

p2 <- plot_auc(
  auc_data = auc_data,
  title = "AUC with extrapolation_points = 5",
  extrapolation_points = 5,
  extrapolated_column = "Extrapolated",
  show_mean = TRUE
)

# Save both plots
pdf("temp/plot_auc_comparison.pdf", width = 12, height = 10)
print(p1)
print(p2)
dev.off()

cat("\nTest completed. Check the PDF files in the temp directory for results.\n") 