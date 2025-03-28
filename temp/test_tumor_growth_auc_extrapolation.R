# Test script to verify AUC extrapolation functionality between tumor_growth_statistics and plot_auc
library(mouseExperiment)
library(ggplot2)

# Set seed for reproducibility
set.seed(123)

# Create a synthetic dataset with incomplete data for some subjects
create_test_data <- function(n_per_group = 8) {
  treatments <- c("Control", "Drug A", "Drug B", "Combination")
  time_points <- seq(0, 21, by = 3)
  
  all_data <- data.frame()
  
  for (treatment in treatments) {
    for (id in 1:(n_per_group)) {
      # Determine if this subject will have incomplete data
      dropout <- FALSE
      dropout_day <- max(time_points)
      
      if (treatment == "Control" && id <= 3) {
        # Some controls have incomplete data (early dropouts)
        dropout <- TRUE
        dropout_day <- 12
      } else if (treatment == "Drug A" && id <= 2) {
        # Some Drug A subjects have incomplete data
        dropout <- TRUE
        dropout_day <- 15
      }
      
      # Generate data points for this subject
      subject_data <- data.frame()
      
      for (day in time_points) {
        if (day <= dropout_day) {
          # Generate volume based on treatment
          if (treatment == "Control") {
            # Accelerating growth - quadratic
            volume <- 100 + 5*day + 1.2*day^2
          } else if (treatment == "Drug A") {
            # Linear growth
            volume <- 100 + 15*day
          } else if (treatment == "Drug B") {
            # Slow growth
            volume <- 100 + 10*day
          } else { # Combination
            # Very slow growth
            volume <- 100 + 5*day
          }
          
          # Add noise
          volume <- volume + rnorm(1, 0, 20)
          
          # Add to subject data
          subject_data <- rbind(subject_data, data.frame(
            ID = paste0(treatment, "_", id),
            Treatment = treatment,
            Cage = paste0("Cage", ((id-1) %/% 2) + 1), # 2 mice per cage
            Day = day,
            Volume = max(50, volume)  # Minimum size of 50
          ))
        }
      }
      
      all_data <- rbind(all_data, subject_data)
    }
  }
  
  return(all_data)
}

# Generate test data
test_data <- create_test_data()

# Print summary of data points per ID
data_points <- aggregate(Day ~ ID + Treatment, data = test_data, FUN = length)
names(data_points)[3] <- "NumPoints"
cat("Data points per subject:\n")
print(table(data_points$Treatment, data_points$NumPoints))

# Test 1: Use tumor_growth_statistics with model_type="auc" and extrapolation_points=0
cat("\n\nTest 1: tumor_growth_statistics with model_type='auc' and extrapolation_points=0\n")
results_no_extrap <- tumor_growth_statistics(
  test_data,
  model_type = "auc",
  extrapolation_points = 0,
  verbose = TRUE
)

# Check if AUC data has Extrapolated column
cat("\nChecking for Extrapolated column in results_no_extrap$auc_analysis$individual:\n")
str(results_no_extrap$auc_analysis$individual)

# Test 2: Use tumor_growth_statistics with model_type="auc" and extrapolation_points=3
cat("\n\nTest 2: tumor_growth_statistics with model_type='auc' and extrapolation_points=3\n")
results_with_extrap <- tumor_growth_statistics(
  test_data,
  model_type = "auc",
  extrapolation_points = 3,
  verbose = TRUE
)

# Check extrapolation values
cat("\nSubjects with extrapolation when extrapolation_points=3:\n")
extrap_summary <- aggregate(Extrapolated ~ Treatment, 
                           data = results_with_extrap$auc_analysis$individual, 
                           FUN = function(x) c(n_extrap = sum(x), total = length(x)))
print(extrap_summary)

# Create plots to visualize the differences
cat("\nCreating visualization of AUC data with and without extrapolation...\n")

# Plot 1: Without extrapolation
p1 <- plot_auc(
  auc_data = results_no_extrap$auc_analysis$individual,
  title = "AUC without extrapolation",
  extrapolated_column = "Extrapolated",
  show_mean = TRUE,
  error_bar_type = "SEM"
)

# Plot 2: With extrapolation
p2 <- plot_auc(
  auc_data = results_with_extrap$auc_analysis$individual,
  title = "AUC with extrapolation (using last 3 points)",
  extrapolated_column = "Extrapolated",
  show_mean = TRUE,
  error_bar_type = "SEM"
)

# Save plots
pdf("temp/tumor_growth_auc_extrapolation_test.pdf", width = 10, height = 12)
print(p1)
print(p2)
dev.off()

cat("\nTest completed. Check the PDF file for visualization of results.\n") 