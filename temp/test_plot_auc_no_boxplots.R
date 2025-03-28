# Test script for plot_auc function without boxplots
library(mouseExperiment)
library(ggplot2)

# Set seed for reproducibility
set.seed(123)

# Create synthetic test data
create_test_data <- function() {
  treatments <- c("Control", "Drug A", "Drug B", "Combination")
  n_per_group <- c(12, 10, 10, 8)  # Different sample sizes per group
  
  auc_data <- data.frame()
  
  # Different means and variances per group
  means <- c(150, 110, 90, 65)
  sds <- c(25, 15, 12, 10)
  
  # Generate AUC values for each group with some randomness
  # and flag some points as extrapolated
  for (i in 1:length(treatments)) {
    treatment <- treatments[i]
    n <- n_per_group[i]
    
    # Decide which subjects have extrapolated data
    # For demonstration, we'll make a few in each group extrapolated
    extrapolated <- rep(FALSE, n)
    if (i <= 2) {  # Only in Control and Drug A groups
      extrapolated[sample(1:n, floor(n/4))] <- TRUE
    }
    
    # Generate AUC values with some randomness
    auc_values <- rnorm(n, means[i], sds[i])
    
    # Create data frame for this group
    group_data <- data.frame(
      Subject = paste0(treatment, "_", 1:n),
      Group = treatment,
      AUC = auc_values,
      Extrapolated = extrapolated
    )
    
    # Add to main data frame
    auc_data <- rbind(auc_data, group_data)
  }
  
  return(auc_data)
}

# Generate test data
auc_data <- create_test_data()

# Print summary of data
cat("Summary of test data:\n")
print(table(auc_data$Group))
print(table(auc_data$Group, auc_data$Extrapolated))

# Create color scheme
custom_colors <- c(
  "Control" = "#1f77b4",      # blue
  "Drug A" = "#ff7f0e",       # orange
  "Drug B" = "#2ca02c",       # green
  "Combination" = "#d62728"   # red
)

# Create and save multiple plots to demonstrate parameters
pdf("temp/plot_auc_no_boxplots_test.pdf", width = 10, height = 8)

# Plot 1: Basic plot - points only, no mean or error bars
p1 <- plot_auc(
  auc_data,
  title = "Points Only - No Mean or Error Bars",
  colors = custom_colors,
  extrapolated_column = "Extrapolated"
)
print(p1)

# Plot 2: With mean but no error bars
p2 <- plot_auc(
  auc_data,
  title = "With Mean Lines (show_mean = TRUE)",
  colors = custom_colors,
  show_mean = TRUE,
  error_bar_type = "none",
  extrapolated_column = "Extrapolated"
)
print(p2)

# Plot 3: With mean and standard error (SEM) error bars
p3 <- plot_auc(
  auc_data,
  title = "With Mean and SEM Error Bars",
  colors = custom_colors,
  show_mean = TRUE,
  error_bar_type = "SEM",
  extrapolated_column = "Extrapolated"
)
print(p3)

# Plot 4: With mean and standard deviation (SD) error bars
p4 <- plot_auc(
  auc_data,
  title = "With Mean and SD Error Bars",
  colors = custom_colors,
  show_mean = TRUE,
  error_bar_type = "SD",
  extrapolated_column = "Extrapolated"
)
print(p4)

# Plot 5: With mean and confidence interval (CI) error bars
p5 <- plot_auc(
  auc_data,
  title = "With Mean and 95% CI Error Bars",
  colors = custom_colors,
  show_mean = TRUE,
  error_bar_type = "CI",
  extrapolated_column = "Extrapolated",
  point_size = 3,
  jitter_width = 0.25,
  caption = "Open circles represent extrapolated values, error bars show 95% confidence intervals"
)
print(p5)

dev.off()

cat("\nTest completed. Check the PDF file for visualization of parameter effects.\n") 