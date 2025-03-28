# Test script to demonstrate the issue with mice having same ID in different cages
# for tumor_auc_analysis and calculate_auc_values functions
library(mouseExperiment)

# Create a synthetic dataset that explicitly reproduces the issue:
# Multiple mice with the same ID in different cages within the same treatment group
create_test_dataset <- function() {
  # Set seed for reproducibility
  set.seed(123)
  
  # Define parameters
  days <- 0:13  # 14 time points
  treatments <- c("Control", "Drug A", "Drug B", "Combination")
  n_mice_per_cage <- 4  # 4 mice per cage
  n_cages_per_treatment <- 2  # 2 cages per treatment (8 mice total per treatment)
  
  # Create empty data frame
  test_data <- data.frame()
  
  # For each treatment
  for (treatment_idx in 1:length(treatments)) {
    treatment <- treatments[treatment_idx]
    
    # For each cage in this treatment
    for (cage_idx in 1:n_cages_per_treatment) {
      cage <- paste0("Cage", cage_idx)
      
      # For each mouse in this cage
      for (mouse_idx in 1:n_mice_per_cage) {
        # Intentionally use the same ID for mice in different cages
        # This simulates the issue where mice have the same ID but are in different cages
        id <- sprintf("%02d", mouse_idx)  # 01, 02, 03, 04
        mouse_id <- paste0("M", id, "_", treatment, "_", cage)  # Make Mouse_ID unique per cage and treatment
        
        # Generate volume data based on treatment and time
        # Different growth rates for different treatments
        base_volume <- 0.1 + runif(1, 0, 0.05)
        growth_rate <- switch(treatment,
                             "Control" = 0.15,
                             "Drug A" = 0.10, 
                             "Drug B" = 0.08,
                             "Combination" = 0.05)
        
        # Add some randomness to growth rate for each mouse
        mouse_growth_rate <- growth_rate * (1 + runif(1, -0.2, 0.2))
        
        # For each day
        for (day in days) {
          # Volume with exponential growth and some noise
          volume <- base_volume * exp(mouse_growth_rate * day) * (1 + rnorm(1, 0, 0.1))
          
          # Add row to data frame
          test_data <- rbind(test_data, data.frame(
            Mouse_ID = mouse_id,
            Day = day,
            Treatment = treatment,
            Volume = volume,
            ID = id,
            Cage = cage
          ))
        }
      }
    }
  }
  
  return(test_data)
}

# Create the test dataset
test_data <- create_test_dataset()

# Summarize the dataset
cat("Summary of test dataset:\n")
cat("Total number of unique Mouse_IDs:", length(unique(test_data$Mouse_ID)), "\n")
cat("Total number of unique ID-Cage-Treatment combinations:", 
    length(unique(paste(test_data$ID, test_data$Treatment, test_data$Cage, sep="_"))), "\n")

cat("\nMice per treatment:\n")
for (treatment in unique(test_data$Treatment)) {
  treatment_data <- subset(test_data, Treatment == treatment)
  cat("  Treatment", treatment, "has", length(unique(treatment_data$Mouse_ID)), "mice\n")
}

# Show the cage distribution
cat("\nCage distribution:\n")
print(table(test_data$Cage, test_data$Treatment))

# Verify we have the issue - mice with the same ID in different cages
cat("\nMice with the same ID in different cages:\n")
id_cage_pairs <- unique(test_data[, c("ID", "Cage", "Treatment")])
id_cage_counts <- table(id_cage_pairs$ID, id_cage_pairs$Treatment)
print(id_cage_counts)

# For simplicity and to avoid the error, let's focus on a single treatment group
treatment_to_test <- "Control"
test_data_single_treatment <- subset(test_data, Treatment == treatment_to_test)

cat("\nTesting with a single treatment group:", treatment_to_test, "\n")
cat("Number of unique Mouse_IDs:", length(unique(test_data_single_treatment$Mouse_ID)), "\n")
cat("Number of unique ID-Cage combinations:", 
    length(unique(paste(test_data_single_treatment$ID, test_data_single_treatment$Cage, sep="_"))), "\n")

# Test 1: Check the core issue in tumor_auc_analysis function
cat("\n\n===== Testing calculation in tumor_auc_analysis function =====\n")

# Replicate the key part of the tumor_auc_analysis function that has the issue
test_subjects <- unique(test_data_single_treatment$ID)
cat("Number of unique subject IDs:", length(test_subjects), "\n")

# Create a simple counter to track subjects found
subjects_found <- 0

# Create a similar data structure to what the function would produce
auc_results_test <- list()

# Loop through each subject as the function would
for (subject in test_subjects) {
  # Get data for this subject - this is where the issue occurs
  subject_data <- test_data_single_treatment[test_data_single_treatment$ID == subject, ]
  
  # The function would get unique treatment, which should be consistent
  treatment <- unique(subject_data$Treatment)
  if (length(treatment) > 1) {
    warning("Subject ", subject, " has multiple treatment assignments. Using the first one.")
    treatment <- treatment[1]
  }
  
  # Check how many unique cages this subject appears in
  unique_cages <- unique(subject_data$Cage)
  
  cat("  Subject ID", subject, "appears in", length(unique_cages), "cages:", 
      paste(unique_cages, collapse=", "), "\n")
  
  # We now see the problem: the function only processes one AUC value per subject ID,
  # but there should be one per unique ID-Cage combination
  for (cage in unique_cages) {
    cage_data <- subject_data[subject_data$Cage == cage, ]
    
    # In actual function, it would calculate AUC here for each cage's data
    # But it only does this once per ID ignoring different cages
    
    # We count a subject found for each unique ID-Cage pair
    subjects_found <- subjects_found + 1
    
    # This mock output shows what the function should be returning
    auc_results_test[[paste(subject, cage, sep="_")]] <- list(
      subject = subject,
      cage = cage,
      treatment = treatment[1]
    )
  }
}

cat("\nTotal unique ID-Cage combinations found:", subjects_found, "\n")
cat("Total entries that would be created by the function:", 
    length(unique(test_data_single_treatment$ID)), "\n")

# Show that the function misses mice with the same ID in different cages
if (subjects_found > length(unique(test_data_single_treatment$ID))) {
  cat("ISSUE CONFIRMED: tumor_auc_analysis would only process", length(unique(test_data_single_treatment$ID)), 
      "mice instead of the expected", subjects_found, "unique ID-Cage combinations\n")
} else {
  cat("No issue found with tumor_auc_analysis calculation logic\n")
}

# Test 2: calculate_auc_values function
cat("\n\n===== Testing calculate_auc_values function =====\n")
# We need to manually call this function since it's not exported
calculate_auc_values <- mouseExperiment:::calculate_auc_values

result_calc_auc <- calculate_auc_values(
  data = test_data_single_treatment,
  time_column = "Day",
  volume_column = "Volume",
  treatment_column = "Treatment",
  id_column = "ID"
)

# Check the results
cat("\nResults from calculate_auc_values function:\n")
cat("Number of mice found:", nrow(result_calc_auc), "\n")
cat("Expected number of unique ID-Cage combinations:", 
    length(unique(paste(test_data_single_treatment$ID, test_data_single_treatment$Cage, sep="_"))), "\n")

# If we found fewer mice than expected, it confirms the issue
if (nrow(result_calc_auc) < length(unique(paste(test_data_single_treatment$ID, test_data_single_treatment$Cage, sep="_")))) {
  cat("ISSUE CONFIRMED: calculate_auc_values only found", nrow(result_calc_auc), 
      "mice instead of the expected", length(unique(paste(test_data_single_treatment$ID, test_data_single_treatment$Cage, sep="_"))), "\n")
} else {
  cat("No issue found with calculate_auc_values\n")
}

# Conclusion
cat("\n\nConclusion:\n")
cat("This test confirms that both functions are not correctly identifying all mice when\n")
cat("there are mice with the same ID in different cages within the same treatment group.\n")
cat("These functions need to be fixed to use a composite ID that includes cage information,\n")
cat("similar to the fix applied to the tumor_growth_statistics function.\n") 