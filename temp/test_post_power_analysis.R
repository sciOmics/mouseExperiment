library(mouseExperiment)

# Test post_power_analysis function directly
cat("Testing post_power_analysis function...\n\n")

# Create minimal synthetic test data
set.seed(123)
n_mice_per_group <- 5
treatments <- c("Control", "Drug A", "Drug B")
days <- 1:10

# Create mouse IDs for each treatment group
mice_ids <- unlist(lapply(1:length(treatments), function(t) {
  paste0(treatments[t], "_M", 1:n_mice_per_group)
}))

# Create cages (2 mice per cage)
cages <- rep(1:ceiling(length(mice_ids)/2), each = 2)[1:length(mice_ids)]

# Create mice data frame
mice_df <- data.frame(
  ID = mice_ids,
  Treatment = rep(treatments, each = n_mice_per_group),
  Cage = cages,
  stringsAsFactors = FALSE
)

# Create complete dataset with all time points
test_data <- expand.grid(
  ID = mice_ids,
  Day = days,
  stringsAsFactors = FALSE
)

# Merge to add treatment and cage info
test_data <- merge(test_data, mice_df, by = "ID")

# Add volume data with treatment effects
test_data$Volume <- with(test_data, {
  base <- 100 + 10 * Day
  
  # Add treatment effects
  effect <- ifelse(Treatment == "Control", 0,
                 ifelse(Treatment == "Drug A", -20, -40))
  
  # Add random noise
  base + effect + rnorm(nrow(test_data), 0, 10)
})

# Print a sample of the data
cat("Sample of test data:\n")
print(head(test_data))

# Test with AUC method
cat("\nTesting with AUC method...\n")
result_auc <- tryCatch({
  post_power_analysis(
    data = test_data,
    effect_sizes = c(0.5, 0.8, 1.0, 1.2), 
    method = "auc"
  )
}, error = function(e) {
  cat("Error with AUC method:", e$message, "\n")
  return(NULL)
})

cat("\nAUC Result is NULL:", is.null(result_auc), "\n")

if (!is.null(result_auc)) {
  cat("\nAUC Power Analysis Results:\n")
  print(result_auc$power_analysis)
  
  cat("\nAUC Sample Size Recommendations:\n")
  print(result_auc$sample_size_recommendations)
}

# Test with parametric method
cat("\nTesting with parametric method...\n")
result_param <- tryCatch({
  post_power_analysis(
    data = test_data,
    effect_sizes = c(0.5, 0.8, 1.0, 1.2), 
    method = "parametric"
  )
}, error = function(e) {
  cat("Error with parametric method:", e$message, "\n")
  return(NULL)
})

cat("\nParametric Result is NULL:", is.null(result_param), "\n")

if (!is.null(result_param)) {
  cat("\nParametric Power Analysis Results:\n")
  print(result_param$power_analysis)
  
  cat("\nParametric Sample Size Recommendations:\n")
  print(result_param$sample_size_recommendations)
}

cat("\nTest completed.\n") 