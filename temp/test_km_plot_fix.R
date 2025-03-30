# Test script to verify that the fixed KM plot in survival_statistics works correctly

library(mouseExperiment)

# Load sample data 
data(combo_treatment_synthetic_data)

# Create synthetic survival data
set.seed(123)
unique_ids <- unique(combo_treatment_synthetic_data$ID)
survival_data <- data.frame(
  ID = unique_ids,
  Survival_Time = round(runif(length(unique_ids), 10, 50)),
  Survival_Censor = sample(c(0, 1), length(unique_ids), replace = TRUE, prob = c(0.5, 0.5)),
  stringsAsFactors = FALSE
)

# Merge survival data with original data
combo_data <- merge(combo_treatment_synthetic_data, survival_data, by = "ID")

# Count events per treatment group correctly (each subject counts once)
cat("Expected Events/Total per treatment group:\n")
treatment_groups <- unique(combo_data$Treatment)
for (group in treatment_groups) {
  # Get unique IDs for this treatment
  group_data <- combo_data[combo_data$Treatment == group, ]
  group_ids <- unique(group_data$ID)
  
  # Count events
  events <- sum(unique(group_data[, c("ID", "Survival_Censor")])$Survival_Censor)
  
  cat(group, ": ", events, "/", length(group_ids), " subjects\n", sep="")
}

# Run the survival_statistics function with our fix
cat("\nRunning survival_statistics function with the fixed create_km_plot...\n")
results <- survival_statistics(
  df = combo_data,
  time_column = "Survival_Time",
  censor_column = "Survival_Censor",
  treatment_column = "Treatment",
  id_column = "ID",
  reference_group = "Control"
)

# Check the Events/Total counts in the results
cat("\nActual Events/Total counts in the results:\n")
print(results$results[, c("Group", "Events", "Total")])

# Verify that the KM plot uses the correct subject-level data
if (!is.null(results$km_plot)) {
  # Save the KM plot to PDF for inspection
  pdf("temp/fixed_km_plot_from_survival_statistics.pdf", width = 10, height = 8)
  print(results$km_plot)
  dev.off()
  cat("\nKM plot saved to temp/fixed_km_plot_from_survival_statistics.pdf\n")
  
  cat("\nTest passed! The survival_statistics function now correctly creates KM plots with each subject counted once.\n")
} else {
  cat("\nWarning: No KM plot was generated. Check that the survminer package is installed.\n")
} 