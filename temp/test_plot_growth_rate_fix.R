# Test script to verify the fix in the plot_growth_rate function
# that accounts for mice with the same ID in different cages
library(mouseExperiment)

# Load the synthetic dataset
data(combo_treatment_synthetic_data)

# Run tumor_growth_statistics to get growth rates data
cat("Running tumor_growth_statistics to generate growth rates data...\n")
results <- tumor_growth_statistics(
  df = combo_treatment_synthetic_data,
  time_column = "Day",
  volume_column = "Volume", 
  treatment_column = "Treatment",
  id_column = "ID",
  cage_column = "Cage"
)

# Extract the growth rates data
growth_data <- results$growth_rates

# Check the source data for unique mice
cat("\nChecking source data for unique mice:\n")
for (treatment in unique(combo_treatment_synthetic_data$Treatment)) {
  treatment_data <- subset(combo_treatment_synthetic_data, Treatment == treatment)
  n_ids <- length(unique(treatment_data$ID))
  n_composite_ids <- length(unique(paste(treatment_data$ID, treatment_data$Treatment, treatment_data$Cage, sep="_")))
  
  cat("Treatment", treatment, "has:\n")
  cat("  - Unique IDs:", n_ids, "\n")
  cat("  - Unique ID-Treatment-Cage combinations:", n_composite_ids, "\n")
}

# Check the processed growth_data
cat("\nChecking growth_data structure:\n")
cat("Total rows in growth_data:", nrow(growth_data), "\n")
cat("Columns in growth_data:", paste(colnames(growth_data), collapse=", "), "\n")

# Create a plot using the fixed plot_growth_rate function
cat("\nCreating growth rate plot...\n")
p <- plot_growth_rate(
  growth_data,
  title = "Tumor Growth Rates with Cage Information Included"
)

# Save the plot to a PDF file for visual inspection
cat("\nSaving plot to PDF for visual inspection...\n")
pdf("temp/plot_growth_rate_test.pdf", width = 10, height = 8)
print(p)
dev.off()

cat("\nTest completed. Check 'temp/plot_growth_rate_test.pdf' for visual verification.\n")
cat("The plot_growth_rate function should now display all 8 mice per treatment group.\n") 