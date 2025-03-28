# Test script for the composite_id fix in tumor_growth_statistics function
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
id_cage_counts <- table(id_cage_pairs$ID)
for (id in names(id_cage_counts)) {
  id_data <- subset(id_cage_pairs, ID == id)
  if (nrow(id_data) > 1) {
    cat("  ID", id, "appears in", nrow(id_data), "different cage-treatment combinations:\n")
    print(id_data)
  }
}

# Extract the current implementation of the tumor_growth_statistics function
cat("\nRunning AUC analysis with modified tumor_growth_statistics function...\n")
result <- tumor_growth_statistics(
  df = test_data,
  model_type = "auc"
)

# Check the AUC analysis results
cat("\nAUC analysis results:\n")
cat("Number of mice found in AUC analysis:", nrow(result$auc_analysis$individual), "\n")
cat("Expected number of unique ID-Cage-Treatment combinations:", 
    length(unique(paste(test_data$ID, test_data$Treatment, test_data$Cage, sep="_"))), "\n")

# Count mice per treatment in AUC analysis
cat("\nMice per treatment in AUC analysis:\n")
treatment_counts <- table(result$auc_analysis$individual$Treatment)
print(treatment_counts)

# Check the number of mice detected per treatment group
cat("\nDetected mice per treatment vs. expected:\n")
for (treatment in unique(test_data$Treatment)) {
  treatment_data <- subset(test_data, Treatment == treatment)
  expected_count <- length(unique(paste(treatment_data$ID, treatment_data$Cage, sep="_")))
  
  treatment_auc <- subset(result$auc_analysis$individual, Treatment == treatment)
  actual_count <- nrow(treatment_auc)
  
  cat("  Treatment", treatment, ": Found", actual_count, "mice, Expected", expected_count, "mice\n")
}

# Display the first few rows of AUC data
cat("\nFirst few rows of AUC data:\n")
print(head(result$auc_analysis$individual))

# Conclusion
cat("\nFix validation result:\n")
expected_count <- length(unique(paste(test_data$ID, test_data$Treatment, test_data$Cage, sep="_")))
if (nrow(result$auc_analysis$individual) == expected_count) {
  cat("SUCCESS: The fix correctly identifies all", expected_count, "unique mice.\n")
  cat("  Found:", nrow(result$auc_analysis$individual), "mice\n")
  cat("  Expected:", expected_count, "mice\n")
} else {
  cat("ERROR: The fix did not correctly identify all mice.\n")
  cat("  Found:", nrow(result$auc_analysis$individual), "mice\n")
  cat("  Expected:", expected_count, "mice\n")
}

cat("\nTest completed.\n") 