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
          # Generate volume based on treatment with different growth patterns
          # Control: accelerating growth
          # Drug A: linear growth
          # Drug B: decelerating growth
          # Combination: minimal growth
          
          base_volume <- 100  # Initial tumor size
          
          # Growth patterns by treatment
          if (treatment == "Control") {
            # Accelerating growth - quadratic
            volume <- base_volume + 10 * day + 1.5 * day^2
          } else if (treatment == "Drug A") {
            # Linear growth - constant rate
            volume <- base_volume + 15 * day
          } else if (treatment == "Drug B") {
            # Decelerating growth - logarithmic
            volume <- base_volume + 40 * log(day + 1)
          } else { # Combination
            # Very slow growth - almost plateau
            volume <- base_volume + 5 * day + 1 * log(day + 1)
          }
          
          # Add randomness
          noise <- rnorm(1, mean = 0, sd = 20)
          
          volume <- volume + noise
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

# Manually examine a dropped-out mouse to understand its growth pattern
control_dropouts <- data_points_summary[data_points_summary$Treatment == "Control" & 
                                      data_points_summary$NumDataPoints < 8, "ID"]
if (length(control_dropouts) > 0) {
  dropout_mouse <- control_dropouts[1]
  dropout_data <- tumor_data[tumor_data$ID == dropout_mouse, ]
  
  cat("\nData for a dropout mouse (", dropout_mouse, "):\n")
  print(dropout_data)
  
  # Plot the actual data points
  p_dropout <- ggplot(dropout_data, aes(x = Day, y = Volume)) +
    geom_point(size = 3) +
    geom_line() +
    labs(title = paste("Growth pattern for", dropout_mouse)) +
    theme_classic()
  
  pdf("temp/dropout_mouse_data.pdf", width = 8, height = 6)
  print(p_dropout)
  dev.off()
}

# Test with different extrapolation_points settings
cat("\nRunning tumor_auc_analysis with extrapolation_points = 2 (using last 2 points)...\n")
auc_results_2pts <- tumor_auc_analysis(tumor_data, extrapolation_points = 2)

cat("\nRunning tumor_auc_analysis with extrapolation_points = 3 (using last 3 points)...\n")
auc_results_3pts <- tumor_auc_analysis(tumor_data, extrapolation_points = 3)

cat("\nRunning tumor_auc_analysis with extrapolation_points = 4 (using last 4 points)...\n")
auc_results_4pts <- tumor_auc_analysis(tumor_data, extrapolation_points = 4)

# Extract AUC values for incomplete mice to compare extrapolation differences
incomplete_mice <- data_points_summary$ID[data_points_summary$NumDataPoints < 8]
extrapolation_comparison <- data.frame()

for (mouse_id in incomplete_mice) {
  # Extract the treatment from the ID
  treatment <- tumor_data$Treatment[tumor_data$ID == mouse_id][1]
  
  # Find the mouse's AUC values with different extrapolation settings
  auc_2pts <- auc_results_2pts$auc_data$AUC[auc_results_2pts$auc_data$Subject == mouse_id]
  auc_3pts <- auc_results_3pts$auc_data$AUC[auc_results_3pts$auc_data$Subject == mouse_id]
  auc_4pts <- auc_results_4pts$auc_data$AUC[auc_results_4pts$auc_data$Subject == mouse_id]
  
  # Check if mouse was extrapolated
  extrapolated_2pts <- auc_results_2pts$auc_data$Extrapolated[auc_results_2pts$auc_data$Subject == mouse_id]
  extrapolated_3pts <- auc_results_3pts$auc_data$Extrapolated[auc_results_3pts$auc_data$Subject == mouse_id]
  extrapolated_4pts <- auc_results_4pts$auc_data$Extrapolated[auc_results_4pts$auc_data$Subject == mouse_id]
  
  # Add to comparison dataframe
  extrapolation_comparison <- rbind(extrapolation_comparison, data.frame(
    ID = mouse_id,
    Treatment = treatment,
    Points = data_points_summary$NumDataPoints[data_points_summary$ID == mouse_id],
    AUC_2pts = auc_2pts,
    AUC_3pts = auc_3pts,
    AUC_4pts = auc_4pts,
    Extrapolated_2pts = extrapolated_2pts,
    Extrapolated_3pts = extrapolated_3pts,
    Extrapolated_4pts = extrapolated_4pts
  ))
}

cat("\nExtrapolation comparison for incomplete mice:\n")
print(extrapolation_comparison)

# Create a plot to visualize the differences in AUC by extrapolation points
extrapolation_long <- reshape2::melt(extrapolation_comparison, 
                                    id.vars = c("ID", "Treatment", "Points"),
                                    measure.vars = c("AUC_2pts", "AUC_3pts", "AUC_4pts"),
                                    variable.name = "Extrapolation",
                                    value.name = "AUC")

# Clean up the Extrapolation names
extrapolation_long$Extrapolation <- factor(extrapolation_long$Extrapolation,
                                         levels = c("AUC_2pts", "AUC_3pts", "AUC_4pts"),
                                         labels = c("2 points", "3 points", "4 points"))

# Create a visualization of the different extrapolation results
p_compare <- ggplot(extrapolation_long, aes(x = Extrapolation, y = AUC, group = ID, color = Treatment)) +
  geom_line() +
  geom_point(size = 3) +
  facet_wrap(~Treatment, scales = "free_y") +
  labs(title = "AUC Values by Extrapolation Points Setting",
       subtitle = "Showing how using different numbers of points affects the extrapolation",
       x = "Number of Points Used for Extrapolation",
       y = "AUC Value") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the comparison plot
pdf("temp/extrapolation_comparison.pdf", width = 10, height = 8)
print(p_compare)
dev.off()

# Test with direct visualization
cat("\nCreating visualizations to show extrapolation differences...\n")

# Create plots with different extrapolation_points settings
p1 <- plot_auc(
  auc_data = auc_results_2pts$auc_data,
  title = "AUC with extrapolation_points = 2 (using last 2 points)",
  extrapolated_column = "Extrapolated",
  show_mean = TRUE
)

p2 <- plot_auc(
  auc_data = auc_results_3pts$auc_data,
  title = "AUC with extrapolation_points = 3 (using last 3 points)",
  extrapolated_column = "Extrapolated",
  show_mean = TRUE
)

p3 <- plot_auc(
  auc_data = auc_results_4pts$auc_data,
  title = "AUC with extrapolation_points = 4 (using last 4 points)",
  extrapolated_column = "Extrapolated",
  show_mean = TRUE
)

# Save all plots
pdf("temp/extrapolation_visualization.pdf", width = 12, height = 15)
print(p1)
print(p2)
print(p3)
dev.off()

cat("\nTest completed. Check the PDF files in the temp directory for results.\n") 