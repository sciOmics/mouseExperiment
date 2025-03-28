# Improved fix for the tumor_growth_statistics function
library(mouseExperiment)

# Create a function that allows us to debug the issue by tracing the execution
debug_composite_id_creation <- function() {
  # Create a small dataset to show the issue
  test_data <- data.frame(
    Mouse_ID = paste0("M", rep(1:4, each = 2, times = 4)),
    ID = rep(1:4, each = 2, times = 4),
    Treatment = rep(c("Control", "Drug A", "Drug B", "Combination"), each = 8),
    Cage = rep(rep(c("Cage1", "Cage2"), each = 4), times = 4),
    Day = rep(0:1, times = 32),
    Volume = runif(64, 0.1, 0.5)
  )
  
  # Show the dataset structure
  cat("Test dataset structure:\n")
  print(head(test_data, 10))
  
  # Current implementation
  composite_id <- paste(test_data$ID, test_data$Treatment, test_data$Cage, sep = "_")
  unique_composite_ids <- unique(composite_id)
  cat("\nNumber of unique composite IDs:", length(unique_composite_ids), "\n")
  cat("First few unique composite IDs:\n")
  print(head(unique_composite_ids))
  
  # The issue: When we use the same ID across different cages and treatments,
  # we need to make sure we create a truly unique composite ID that works
  # for our specific dataset structure
  
  # Print a summary of unique ID-Cage-Treatment combinations
  id_cage_treatment <- unique(test_data[, c("ID", "Cage", "Treatment")])
  cat("\nUnique ID-Cage-Treatment combinations:", nrow(id_cage_treatment), "\n")
  
  # Count how many times each ID appears across all cages and treatments
  id_counts <- table(test_data$ID)
  cat("\nID occurrence counts:\n")
  print(id_counts)
  
  # Print the first few rows of unique ID-Cage-Treatment combinations
  cat("\nSample of ID-Cage-Treatment combinations:\n")
  print(head(id_cage_treatment, 10))
  
  # Problem diagnosis: The function isn't distinguishing between cage values when
  # looking at the same ID and Treatment, causing mice to be collapsed together
  
  cat("\nDetailed investigation of the issue:\n")
  # Inspect how composite_id relates to unique Mouse_IDs
  mouse_id_to_composite <- data.frame(
    Mouse_ID = test_data$Mouse_ID, 
    ID = test_data$ID,
    Cage = test_data$Cage, 
    Treatment = test_data$Treatment,
    CompositeID = paste(test_data$ID, test_data$Treatment, test_data$Cage, sep = "_")
  )
  mouse_id_to_composite <- mouse_id_to_composite[!duplicated(mouse_id_to_composite[c("Mouse_ID", "CompositeID")]),]
  cat("Mapping of Mouse_ID to CompositeID:\n")
  print(mouse_id_to_composite)
}

# Run the debug function
cat("=== Debugging the composite_id creation ===\n\n")
debug_composite_id_creation()

# Propose the improved fix
cat("\n\n=== IMPROVED FIX ===\n\n")
cat("The problem is that we need to use a more robust method to identify\n")
cat("unique mice in our dataset. Here's the improved fix for tumor_growth_statistics.R:\n\n")
cat("1. Correctly handle composite_id creation at line ~491 by using a better approach:\n\n")
cat("   Original code:\n")
cat("   # Calculate AUC for each individual (Composite ID ensures unique subject-treatment combinations)\n")
cat("   composite_id <- paste(auc_df[[id_column]], auc_df[[treatment_column]], sep = \"_\")\n\n")
cat("   First attempt at fix:\n")
cat("   # Calculate AUC for each individual using ID + Treatment + Cage for unique identification\n")
cat("   composite_id <- paste(auc_df[[id_column]], auc_df[[treatment_column]], auc_df[[cage_column]], sep = \"_\")\n\n")
cat("   Improved fix (more robust approach):\n")
cat("   # For each unique ID-Treatment-Cage combination, create a unique identifier\n")
cat("   # This ensures proper distinction of mice even when they share the same ID but are in different cages\n")
cat("   # First find all unique combinations\n")
cat("   unique_combinations <- unique(auc_df[, c(id_column, treatment_column, cage_column)])\n")
cat("   # Create a mapping of these combinations to sequential numbers\n")
cat("   unique_combinations$unique_id <- 1:nrow(unique_combinations)\n")
cat("   # Merge back with the original data to assign the correct unique ID to each row\n")
cat("   auc_df_with_id <- merge(auc_df, unique_combinations, by=c(id_column, treatment_column, cage_column))\n")
cat("   # Use this unique_id for processing\n")
cat("   composite_id <- auc_df_with_id$unique_id\n\n")
cat("2. Update the loop that iterates through unique_id to handle the proper identifier:\n\n")
cat("   for (unique_id in unique(composite_id)) {\n")
cat("     # Extract data for this unique ID\n")
cat("     subject_data <- auc_df_with_id[auc_df_with_id$unique_id == unique_id,]\n")
cat("     # Get the original ID, treatment, and cage for this unique subject\n")
cat("     actual_id <- subject_data[[id_column]][1]  # Just use the first row's ID\n")
cat("     treatment <- subject_data[[treatment_column]][1]  # First row's treatment\n\n")
cat("3. The final fix is to properly create and handle these unique identifiers,\n")
cat("   ensuring each mouse (unique by ID-Treatment-Cage) is correctly processed.\n\n")
cat("This approach is more robust as it doesn't rely on string concatenation\n")
cat("for identification, but creates a truly unique identifier for each\n")
cat("ID-Treatment-Cage combination in the dataset.\n")

cat("\n=== Implementation Steps ===\n\n")
cat("1. Back up the current version of tumor_growth_statistics.R\n")
cat("2. Implement the improved fix in the AUC calculation section\n")
cat("3. Update the corresponding code in the extraction of subject data\n")
cat("4. Test with datasets having duplicate IDs across cages\n")
cat("5. Update the documentation to explain the approach used\n") 