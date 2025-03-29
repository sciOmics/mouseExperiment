# Test script to verify the fixed Kaplan-Meier plot in survival_statistics
library(mouseExperiment)

# Load the synthetic dataset
data(combo_treatment_synthetic_data)

# Create a clean copy of the data for testing
test_data <- combo_treatment_synthetic_data

# Add a Survival_Censor column if it doesn't exist
if (!("Survival_Censor" %in% colnames(test_data))) {
  cat("Adding synthetic Survival_Censor column...\n")
  
  # Get unique mice
  mouse_ids <- unique(paste(test_data$ID, test_data$Cage, test_data$Treatment, sep="_"))
  
  # Create a lookup table for censor status
  # Randomly censor ~30% of mice, mark the rest as events (1)
  set.seed(123) # For reproducibility
  censor_lookup <- data.frame(
    id_cage_treatment = mouse_ids,
    censor_status = sample(c(0, 1), length(mouse_ids), replace = TRUE, 
                          prob = c(0.3, 0.7))
  )
  
  # Create a composite ID in the dataset for joining
  test_data$id_cage_treatment <- paste(test_data$ID, test_data$Cage, 
                                     test_data$Treatment, sep="_")
  
  # Merge censor status
  test_data$Survival_Censor <- censor_lookup$censor_status[
    match(test_data$id_cage_treatment, censor_lookup$id_cage_treatment)
  ]
  
  # Clean up
  test_data$id_cage_treatment <- NULL
}

# Ensure we have the necessary columns
cat("Checking required columns in dataset:\n")
cat("- Survival_Censor column exists:", "Survival_Censor" %in% colnames(test_data), "\n")
cat("- Treatment column exists:", "Treatment" %in% colnames(test_data), "\n")
cat("- Day column exists:", "Day" %in% colnames(test_data), "\n")
cat("- ID column exists:", "ID" %in% colnames(test_data), "\n")
cat("- Cage column exists:", "Cage" %in% colnames(test_data), "\n")

# Show treatment groups
cat("\nTreatment groups in dataset:\n")
cat(paste("- ", unique(test_data$Treatment), "\n"))

# Summarize Survival_Censor values
cat("\nSurvival_Censor summary:\n")
censor_table <- table(test_data$Survival_Censor, test_data$Treatment)
print(censor_table)

# Run survival_statistics function
cat("\nRunning survival_statistics function...\n")
survival_results <- survival_statistics(
  df = test_data,
  time_column = "Day",
  censor_column = "Survival_Censor",
  treatment_column = "Treatment",
  cage_column = "Cage",
  id_column = "ID",
  reference_group = "Control"
)

# Check if survival_results contains the km_plot
cat("\nChecking if Kaplan-Meier plot was successfully created:\n")
cat("- km_plot exists in results:", "km_plot" %in% names(survival_results), "\n")

if ("km_plot" %in% names(survival_results)) {
  cat("- km_plot is not NULL:", !is.null(survival_results$km_plot), "\n")
  
  # Save the plot to a PDF file for visual inspection
  if (!is.null(survival_results$km_plot)) {
    cat("\nSaving Kaplan-Meier plot to PDF file...\n")
    pdf("temp/km_plot_test.pdf", width = 10, height = 8)
    print(survival_results$km_plot)
    dev.off()
    cat("Plot saved to 'temp/km_plot_test.pdf'\n")
  }
}

# Check other elements of the results
cat("\nOther elements in the survival_statistics results:\n")
cat(paste("- ", names(survival_results), "\n"))

cat("\nTest completed successfully.\n") 