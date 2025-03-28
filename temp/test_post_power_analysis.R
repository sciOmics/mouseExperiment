# Test script for post_power_analysis function
# This script tests:
# 1. The AUC method with proper data validation
# 2. The simulation method implementation

library(mouseExperiment)

# Set seed for reproducibility
set.seed(123)

# Create a synthetic dataset that ensures multiple treatment groups
generate_test_data <- function(n_mice_per_group = 8) {
  # Define treatment groups - we need at least 2 for the AUC method
  treatment_groups <- c("Control", "Drug A", "Drug B", "Combination")
  n_groups <- length(treatment_groups)
  
  # Define time points
  time_points <- 0:14
  n_time_points <- length(time_points)
  
  # Create empty data frame
  test_data <- data.frame(
    Mouse = character(),
    Day = numeric(),
    Treatment = character(),
    Volume = numeric(),
    Cage = numeric(),  # Add cage column for simulation method
    stringsAsFactors = FALSE
  )
  
  # Generate data for each treatment group
  for (t_idx in 1:n_groups) {
    treatment <- treatment_groups[t_idx]
    
    # Parameters for growth curves (different for each treatment)
    baseline <- 100
    growth_rate <- switch(treatment,
                         "Control" = 0.20,
                         "Drug A" = 0.15,
                         "Drug B" = 0.12,
                         "Combination" = 0.08)
    
    # Variance parameters (different for each treatment)
    # We'll make these heterogeneous to test robustness
    between_subj_sd <- switch(treatment,
                            "Control" = 15,
                            "Drug A" = 10,
                            "Drug B" = 12,
                            "Combination" = 8)
    
    within_subj_sd <- switch(treatment,
                           "Control" = 12,
                           "Drug A" = 8,
                           "Drug B" = 7,
                           "Combination" = 5)
    
    # Generate data for each mouse
    for (m in 1:n_mice_per_group) {
      mouse_id <- paste0(treatment, "_Mouse", m)
      
      # Assign to cage (2 mice per cage)
      cage_id <- floor((m-1)/2) + 1 + (t_idx-1) * ceiling(n_mice_per_group/2)
      
      # Generate random subject effect
      subject_effect <- rnorm(1, 0, between_subj_sd)
      
      # Generate volume for each time point
      for (t in time_points) {
        # Exponential growth model with random effects
        mean_volume <- baseline * exp(growth_rate * t) + subject_effect
        volume <- max(1, rnorm(1, mean_volume, within_subj_sd))
        
        # Add to data frame
        test_data <- rbind(test_data, data.frame(
          Mouse = mouse_id,
          Day = t,
          Treatment = treatment,
          Volume = volume,
          Cage = cage_id,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  return(test_data)
}

# Generate test data
test_data <- generate_test_data(n_mice_per_group = 6)

# View the structure of the data
str(test_data)

# Calculate AUC values for each mouse directly
calculate_auc_data <- function(data, time_col = "Day", volume_col = "Volume", id_col = "Mouse", treatment_col = "Treatment") {
  # Get unique mice
  unique_mice <- unique(data[[id_col]])
  
  # Initialize results data frame
  auc_data <- data.frame(
    ID = character(),
    Treatment = character(),
    AUC = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Calculate AUC for each mouse
  for (mouse in unique_mice) {
    # Get data for this mouse
    mouse_data <- data[data[[id_col]] == mouse, ]
    
    # Sort by time
    mouse_data <- mouse_data[order(mouse_data[[time_col]]), ]
    
    # Get treatment
    treatment <- unique(mouse_data[[treatment_col]])
    if (length(treatment) > 1) {
      warning("Mouse ", mouse, " has multiple treatments. Using first one.")
      treatment <- treatment[1]
    }
    
    # Calculate AUC using trapezoid method
    time_points <- mouse_data[[time_col]]
    volumes <- mouse_data[[volume_col]]
    
    if (length(time_points) < 2) {
      warning("Not enough time points for mouse ", mouse)
      next
    }
    
    auc_value <- 0
    for (i in 2:length(time_points)) {
      # Area of trapezoid = (y1 + y2) * (x2 - x1) / 2
      auc_value <- auc_value + (volumes[i-1] + volumes[i]) * (time_points[i] - time_points[i-1]) / 2
    }
    
    # Add to results
    auc_data <- rbind(auc_data, data.frame(
      ID = mouse,
      Treatment = treatment,
      AUC = auc_value,
      stringsAsFactors = FALSE
    ))
  }
  
  # Change column names to match expected
  auc_data_formatted <- data.frame(
    ID = auc_data$ID,
    AUC = auc_data$AUC,
    stringsAsFactors = FALSE
  )
  auc_data_formatted[[treatment_col]] <- auc_data$Treatment
  
  # Return auc data
  return(auc_data_formatted)
}

# Calculate AUC values
auc_data <- calculate_auc_data(test_data)

# Create data with AUC already calculated
model_results <- list(
  auc_analysis = list(
    individual = auc_data
  )
)

# Test 1: AUC method
cat("\n--- Testing post_power_analysis with AUC method ---\n")
auc_results <- tryCatch({
  post_power_analysis(
    data = model_results, # Use pre-calculated AUC data
    method = "auc",
    time_column = "Day",
    volume_column = "Volume",
    treatment_column = "Treatment",
    id_column = "Mouse",
    effect_sizes = c(0.5, 0.8, 1.0, 1.5)
  )
}, error = function(e) {
  cat("ERROR in AUC method:", e$message, "\n")
  return(NULL)
})

# Print AUC results if successful
if (!is.null(auc_results)) {
  cat("\nAUC Power Analysis Results:\n")
  print(auc_results$power_analysis)
  
  cat("\nSample Size Recommendations:\n")
  print(auc_results$sample_size_recommendations)
}

# Test 2: Simulation method
cat("\n--- Testing post_power_analysis with simulation method ---\n")
sim_results <- tryCatch({
  post_power_analysis(
    data = test_data,
    method = "simulation",
    time_column = "Day",
    volume_column = "Volume",
    treatment_column = "Treatment",
    id_column = "Mouse",
    effect_sizes = c(0.5, 0.8, 1.0),
    n_simulations = 100  # Use a smaller number for faster testing
  )
}, error = function(e) {
  cat("ERROR in simulation method:", e$message, "\n")
  return(NULL)
})

# Print simulation results if successful
if (!is.null(sim_results)) {
  cat("\nSimulation Power Analysis Results:\n")
  print(sim_results$power_analysis)
  
  cat("\nSample Size Recommendations:\n")
  print(sim_results$sample_size_recommendations)
}

# Create a plot if either method worked
if (!is.null(auc_results) && !is.null(auc_results$plots$power_curve)) {
  pdf("temp/auc_power_curve.pdf", width = 8, height = 6)
  print(auc_results$plots$power_curve)
  dev.off()
  cat("\nSaved AUC power curve plot to temp/auc_power_curve.pdf\n")
}

if (!is.null(sim_results) && !is.null(sim_results$plots$power_curve)) {
  pdf("temp/sim_power_curve.pdf", width = 8, height = 6)
  print(sim_results$plots$power_curve)
  dev.off()
  cat("\nSaved simulation power curve plot to temp/sim_power_curve.pdf\n")
}

# Test edge case: Small sample size (just 1 mouse per group)
cat("\n--- Testing with small sample size (should fail validation) ---\n")
small_data <- generate_test_data(n_mice_per_group = 1)
small_auc_data <- calculate_auc_data(small_data)
small_model_results <- list(
  auc_analysis = list(
    individual = small_auc_data
  )
)

small_result <- tryCatch({
  post_power_analysis(
    data = small_model_results,
    method = "auc",
    time_column = "Day", 
    volume_column = "Volume",
    treatment_column = "Treatment",
    id_column = "Mouse"
  )
}, error = function(e) {
  cat("Expected error with small sample:", e$message, "\n")
  return(NULL)
})

cat("\n--- Tests completed ---\n") 