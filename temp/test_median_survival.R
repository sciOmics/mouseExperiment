# Test script for median survival calculation fix
library(mouseExperiment)
library(survival)

# Set seed for reproducibility
set.seed(42)

# Create a synthetic dataset with varying event rates
create_test_data <- function() {
  # Define treatment groups with different survival profiles
  n_per_group <- 12
  
  # Create data frame to store results
  test_data <- data.frame()
  
  # Group 1: High event rate (>50% events) - should have calculable median
  # Group 2: Medium event rate (=50% events) - should have calculable median
  # Group 3: Low event rate (<50% events) - "Not reached" median
  
  # Survival times for each group (different distributions)
  time_distributions <- list(
    "High_Events" = function(n) {rexp(n, rate=0.1) + 10},  # Short survival
    "Medium_Events" = function(n) {rexp(n, rate=0.05) + 20},  # Medium survival
    "Low_Events" = function(n) {rexp(n, rate=0.02) + 30}   # Long survival
  )
  
  # Censoring rates for each group
  # Higher values mean more censoring (fewer events)
  censoring_rates <- list(
    "High_Events" = 0.2,    # 80% events (high)
    "Medium_Events" = 0.5,  # 50% events (medium)
    "Low_Events" = 0.7      # 30% events (low)
  )
  
  # Generate data for each group
  groups <- names(time_distributions)
  for (group_idx in 1:length(groups)) {
    group <- groups[group_idx]
    
    # Generate survival times from distribution
    times <- time_distributions[[group]](n_per_group)
    times <- round(times, 1)  # Round to 1 decimal place
    
    # Generate censoring based on probability
    censor_prob <- censoring_rates[[group]]
    censored <- rbinom(n_per_group, 1, censor_prob) == 1  # TRUE if censored
    event <- !censored  # Event indicator (1 = event, 0 = censored)
    
    # Create IDs
    ids <- paste0(substr(group, 1, 1), 1:n_per_group)
    
    # Create cage assignments (2 mice per cage)
    cages <- paste0("Cage", ceiling((1:n_per_group) / 2))
    
    # Add to data frame
    for (i in 1:n_per_group) {
      test_data <- rbind(test_data, data.frame(
        ID = ids[i],
        Treatment = group,
        Day = times[i],
        Survival_Censor = as.integer(event[i]),  # 1 = event, 0 = censored
        Cage = cages[i]
      ))
    }
  }
  
  return(test_data)
}

# Generate test data
test_data <- create_test_data()

# Print summary of data by group showing censoring
summary_by_group <- function(data) {
  cat("\nData summary by treatment group:\n")
  for (group in unique(data$Treatment)) {
    group_data <- data[data$Treatment == group, ]
    n_subjects <- nrow(group_data)
    n_events <- sum(group_data$Survival_Censor)
    event_rate <- n_events / n_subjects
    
    cat(sprintf("%s: %d subjects, %d events (%.1f%% event rate)\n", 
                group, n_subjects, n_events, event_rate * 100))
  }
}

# Print data summary
summary_by_group(test_data)

# Run survival analysis
cat("\nRunning survival analysis with the fixed median calculation:\n")
surv_results <- survival_statistics(
  df = test_data,
  time_column = "Day",
  censor_column = "Survival_Censor",
  treatment_column = "Treatment",
  cage_column = "Cage",
  id_column = "ID"
)

# Output results
cat("\nVerifying median survival results:\n")
print(surv_results$results[, c("Group", "Median_Survival", "Events", "Total", "Event_Rate")])

# Create Kaplan-Meier curves to visualize
cat("\nCreating Kaplan-Meier curves to verify:\n")
fit <- survival::survfit(Surv(Day, Survival_Censor) ~ Treatment, data = test_data)
print(fit)

# Save KM plot
if (requireNamespace("ggplot2", quietly = TRUE) && 
    requireNamespace("survminer", quietly = TRUE)) {
  km_plot <- survminer::ggsurvplot(
    fit, 
    data = test_data,
    pval = TRUE,
    conf.int = TRUE,
    risk.table = TRUE,
    title = "Testing Median Survival Calculation",
    xlab = "Days",
    ggtheme = ggplot2::theme_bw()
  )
  
  # Save plot to PDF file
  pdf("temp/median_survival_test.pdf", width = 10, height = 8)
  print(km_plot)
  dev.off()
  cat("\nKaplan-Meier plot saved to temp/median_survival_test.pdf\n")
} 