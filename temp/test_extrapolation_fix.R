# Test script to verify that the fixed extrapolation logic in tumor_growth_statistics works correctly
library(mouseExperiment)

# Load the synthetic dataset
data(combo_treatment_synthetic_data)

# Get the maximum day in the dataset
max_day <- max(combo_treatment_synthetic_data$Day)
cat("Maximum day in the dataset:", max_day, "\n")

# Create a modified dataset with some subjects missing the last day
modified_data <- combo_treatment_synthetic_data

# Remove day 26 for selected subjects to force extrapolation
# For the first 2 mice in Control and HDACi groups
control_subjects <- unique(subset(modified_data, Treatment == "Control")$ID)[1:2]
hdaci_subjects <- unique(subset(modified_data, Treatment == "HDACi")$ID)[1:2]

cat("\nRemoving day", max_day, "data for subjects:", 
    paste(c(control_subjects, hdaci_subjects), collapse=", "), "\n")

# Filter out the rows for day 26 for these subjects
modified_data <- modified_data[!(modified_data$Day == max_day & 
                               (modified_data$ID %in% control_subjects | 
                                  modified_data$ID %in% hdaci_subjects)), ]

# Verify the data modification
cat("\nVerifying data modification:\n")
for (id in c(control_subjects, hdaci_subjects)) {
  subject_data <- subset(modified_data, ID == id)
  cat("Subject", id, "max day:", max(subject_data$Day), "\n")
}

# Run tumor_growth_statistics with extrapolation_points > 0
cat("\nRunning tumor_growth_statistics with extrapolation...\n")
results <- tumor_growth_statistics(
  df = modified_data,
  time_column = "Day",
  volume_column = "Volume", 
  treatment_column = "Treatment",
  id_column = "ID",
  cage_column = "Cage",
  model_type = "auc",
  extrapolation_points = 1,  # Request extrapolation
  verbose = TRUE  # Show detailed messages
)

# Check the Last_Day values in the auc_analysis$individual data frame
cat("\nVerifying Last_Day values in results:\n")
auc_data <- results$auc_analysis$individual
days <- sort(unique(auc_data$Last_Day))
cat("Unique Last_Day values:", paste(days, collapse = ", "), "\n")
cat("Maximum Last_Day:", max(auc_data$Last_Day), "\n")

# Check if the maximum Last_Day is less than or equal to the true maximum day
cat("Maximum Last_Day <= true maximum day:", max(auc_data$Last_Day) <= max_day, "\n")

# Check if the Last_Day values are integers
cat("All Last_Day values are integers:", all(auc_data$Last_Day == round(auc_data$Last_Day)), "\n")

# Check which subjects were extrapolated
cat("\nExtrapolation status by subject:\n")
extrapolation_summary <- table(auc_data$Treatment, auc_data$Extrapolated)
print(extrapolation_summary)

# Check if there are any extrapolated points
n_extrapolated <- sum(auc_data$Extrapolated)
cat("\nNumber of extrapolated subjects:", n_extrapolated, "\n")

# For extrapolated subjects, verify their original last observation day
if (n_extrapolated > 0) {
  extrapolated_subjects <- auc_data[auc_data$Extrapolated, ]
  cat("\nDetails for extrapolated subjects:\n")
  cat("ID  Treatment  Last_Day  Extrapolated\n")
  for (i in 1:nrow(extrapolated_subjects)) {
    cat(sprintf("%s  %s  %d  %s\n", 
                extrapolated_subjects$ID[i],
                extrapolated_subjects$Treatment[i],
                extrapolated_subjects$Last_Day[i],
                extrapolated_subjects$Extrapolated[i]))
  }
  
  # Check that the extrapolated subjects match those we removed data for
  extrapolated_ids <- extrapolated_subjects$ID
  cat("\nChecking if extrapolated subjects match those we removed data for:\n")
  cat("Expected subjects with missing day 26:", paste(c(control_subjects, hdaci_subjects), collapse=", "), "\n")
  cat("Actual extrapolated subjects:", paste(extrapolated_ids, collapse=", "), "\n")
  
  all_expected_extrapolated <- all(c(control_subjects, hdaci_subjects) %in% extrapolated_ids)
  only_expected_extrapolated <- all(extrapolated_ids %in% c(control_subjects, hdaci_subjects))
  
  cat("All expected subjects were extrapolated:", all_expected_extrapolated, "\n")
  cat("Only expected subjects were extrapolated:", only_expected_extrapolated, "\n")
}

cat("\nTest completed.\n") 