# Test script to verify that cage information is included in the composite IDs note
# and that the auc_data dataframe includes the Cage column
library(mouseExperiment)

# Load the synthetic dataset
data(combo_treatment_synthetic_data)

# Run tumor_growth_statistics with AUC model type
cat("Running tumor_growth_statistics with AUC model...\n")
results_auc <- tumor_growth_statistics(
  df = combo_treatment_synthetic_data,
  time_column = "Day",
  volume_column = "Volume", 
  treatment_column = "Treatment",
  id_column = "ID",
  cage_column = "Cage",
  model_type = "auc"
)

# Check if the notes mention cage information in composite IDs
cat("\nVerifying that the notes mention cage information in composite IDs:\n")
notes <- results_auc$summary$notes
for (i in 1:length(notes)) {
  cat(i, ": ", notes[i], "\n")
}

cat("\nChecking for mention of 'cage' in composite ID note: ", 
    any(grepl("cage", tolower(notes))), "\n")

# Check if the auc_data dataframe includes a Cage column
cat("\nVerifying that the auc_data dataframe includes a Cage column:\n")
auc_data_columns <- names(results_auc$auc_analysis$individual)
cat("Columns in auc_data: ", paste(auc_data_columns, collapse = ", "), "\n")
cat("'Cage' column exists: ", "Cage" %in% auc_data_columns, "\n")

# Check the first few rows of the auc_data dataframe
cat("\nFirst few rows of the auc_data dataframe:\n")
head_data <- head(results_auc$auc_analysis$individual)
print(head_data)

# Count unique mice by treatment and cage
cat("\nCounting unique mice by treatment and cage:\n")
for (treatment in unique(combo_treatment_synthetic_data$Treatment)) {
  treatment_data <- subset(combo_treatment_synthetic_data, Treatment == treatment)
  for (cage in unique(treatment_data$Cage)) {
    cage_data <- subset(treatment_data, Cage == cage)
    unique_mice <- length(unique(cage_data$ID))
    cat("  Treatment:", treatment, "- Cage:", cage, "has", unique_mice, "unique mice\n")
  }
}

cat("\nTest completed.\n") 