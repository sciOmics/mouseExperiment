# Example script for analyzing combo treatment data
# This script demonstrates how to use the enhanced date handling features

# Load the package
library(mouseExperiment)

# Read combo treatment data
data_path <- "/Users/david/Library/Mobile Documents/com~apple~CloudDocs/Programming/mouseExperiment/data/combo_treatment_synthetic_data.csv"
combo_data <- read.csv(data_path)

# Print information about the data
cat("Combo treatment data format:\n")
cat("First date in data:", as.character(combo_data$Date[1]), "\n")
cat("Column names:", paste(colnames(combo_data), collapse = ", "), "\n")
cat("Treatment groups:", paste(unique(combo_data$Treatment), collapse = ", "), "\n\n")

# Calculate tumor volumes (required for analysis)
combo_data <- calculate_volume(combo_data)

# Calculate days from start date - the function will auto-detect the MM/DD/YYYY format
cat("Calculating days using the auto-detected date format:\n")
combo_data <- calculate_dates(combo_data, 
                             start_date = "03/24/2025")

# Show the resulting data with Day column
cat("\nFirst 5 rows of processed data:\n")
print(head(combo_data, 5))

# Create a simple plot showing tumor growth by treatment group
cat("\nCreating tumor growth plot by treatment group...\n")
result <- plot_tumor_growth(combo_data, 
                           group_column = "Treatment",
                           volume_column = "Volume")

# You can also analyze survival differences between treatment groups
cat("\nAnalyzing survival differences between treatment groups...\n")
result <- survival_statistics(combo_data, 
                             time_column = "Day",
                             censor_column = "Survival_Censor",
                             group_column = "Treatment")

cat("\nAnalysis complete!\n")