# Test script for plot_auc parameters (colors, point_size, jitter_width)
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
pdf("temp/plot_auc_parameters_test.pdf", width = 10, height = 8)

# Plot 1: Basic plot with default parameters
p1 <- plot_auc(auc_data)
print(p1)

# Plot 2: With custom colors
p2 <- plot_auc(
  auc_data,
  title = "AUC with Custom Colors",
  colors = custom_colors,
  extrapolated_column = "Extrapolated"
)
print(p2)

# Plot 3: With larger point size
p3 <- plot_auc(
  auc_data,
  title = "AUC with Larger Points",
  colors = custom_colors,
  point_size = 4,
  extrapolated_column = "Extrapolated"
)
print(p3)

# Plot 4: With wider jitter
p4 <- plot_auc(
  auc_data,
  title = "AUC with Wider Jitter",
  colors = custom_colors,
  jitter_width = 0.4,
  extrapolated_column = "Extrapolated"
)
print(p4)

# Plot 5: Combining all parameters with mean and error bars
p5 <- plot_auc(
  auc_data,
  title = "AUC with All Custom Parameters",
  colors = custom_colors,
  point_size = 3.5,
  jitter_width = 0.25,
  show_mean = TRUE,
  error_bar_type = "SEM",
  extrapolated_column = "Extrapolated",
  caption = "Shows all restored parameters: colors, point_size, and jitter_width"
)
print(p5)

dev.off()

cat("\nTest completed. Check the PDF file for visualization of parameter effects.\n") 