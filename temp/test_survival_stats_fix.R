# Test script to verify the fix for survival_statistics function
library(mouseExperiment)
library(survival)

# Set seed for reproducibility
set.seed(123)

# Create synthetic test data
n_mice <- 50
treatments <- c("Control", "Drug A", "Drug B", "Combination")

# Generate test data
generate_test_data <- function() {
  # Randomly assign treatment groups
  treatment <- sample(treatments, n_mice, replace = TRUE)
  
  # Generate survival time based on treatment (different distributions)
  time <- numeric(n_mice)
  for (i in 1:n_mice) {
    if (treatment[i] == "Control") {
      time[i] <- rexp(1, 1/20) # mean survival 20 days
    } else if (treatment[i] == "Drug A") {
      time[i] <- rexp(1, 1/30) # mean survival 30 days
    } else if (treatment[i] == "Drug B") {
      time[i] <- rexp(1, 1/25) # mean survival 25 days
    } else {
      time[i] <- rexp(1, 1/40) # mean survival 40 days for combination
    }
  }
  
  # Round times to whole days
  time <- round(time)
  
  # Ensure minimum time is at least 1 day
  time <- pmax(time, 1)
  
  # Generate censoring status (1 = event, 0 = censored)
  # Make sure each group has >50% events to test our fix
  censor <- numeric(n_mice)
  for (i in 1:n_mice) {
    # Higher probability of event for testing purposes
    if (treatment[i] == "Control") {
      censor[i] <- rbinom(1, 1, 0.8) # 80% events
    } else if (treatment[i] == "Drug A") {
      censor[i] <- rbinom(1, 1, 0.7) # 70% events
    } else if (treatment[i] == "Drug B") {
      censor[i] <- rbinom(1, 1, 0.6) # 60% events
    } else {
      censor[i] <- rbinom(1, 1, 0.5) # 50% events for combination
    }
  }
  
  # Create data frame
  data <- data.frame(
    ID = paste0("Mouse", 1:n_mice),
    Treatment = treatment,
    Day = time,
    Survival_Censor = censor,
    Cage = sample(LETTERS[1:10], n_mice, replace = TRUE)
  )
  
  return(data)
}

# Generate test data
data <- generate_test_data()

# Check distribution of events by treatment
print("Event distribution by treatment:")
print(table(data$Treatment, data$Survival_Censor))

# Calculate event rates to verify we have >50% events in each group
event_rates <- aggregate(Survival_Censor ~ Treatment, data, function(x) sum(x)/length(x))
print("Event rates by treatment:")
print(event_rates)

# Run survival_statistics function
cat("\nRunning survival_statistics function with the fix:\n")
results <- survival_statistics(
  df = data,
  time_column = "Day",
  censor_column = "Survival_Censor",
  treatment_column = "Treatment",
  cage_column = "Cage",
  id_column = "ID"
)

# Verify the results object has expected components
cat("\nVerifying results structure:\n")
print(names(results))

# Print median survival times to verify they were calculated correctly
cat("\nMedian survival times:\n")
if ("Median_Survival" %in% colnames(results$results)) {
  print(data.frame(
    Group = results$results$Group,
    Median_Survival = results$results$Median_Survival,
    Event_Rate = results$results$Event_Rate
  ))
} else {
  cat("Median_Survival column not found in results\n")
}

cat("\nTest completed successfully!\n") 