library(mouseExperiment)

# Create a synthetic dataset with survival data
set.seed(123)

# Create groups
treatments <- c("Control", "PD1", "HDACi", "HDACi + PD1")
n_per_group <- 8
n_total <- length(treatments) * n_per_group

# Create mouse IDs
mouse_ids <- paste0("M", sprintf("%02d", 1:n_total))

# Generate survival times with HDACi + PD1 having the best survival
control_times <- rexp(n_per_group, rate = 0.03) + 10
pd1_times <- rexp(n_per_group, rate = 0.025) + 12
hdaci_times <- rexp(n_per_group, rate = 0.02) + 15
combo_times <- rexp(n_per_group, rate = 0.01) + 20

# Combine times
times <- c(control_times, pd1_times, hdaci_times, combo_times)

# Generate censoring (1 = event occurred, 0 = censored)
# Higher probability of events in control, lower in combination
control_events <- rbinom(n_per_group, 1, 0.8)
pd1_events <- rbinom(n_per_group, 1, 0.7)
hdaci_events <- rbinom(n_per_group, 1, 0.6)
combo_events <- rbinom(n_per_group, 1, 0.5)

# Combine events
events <- c(control_events, pd1_events, hdaci_events, combo_events)

# Censor times beyond 40 days
times[times > 40] <- 40
events[times == 40] <- 0  # Mark these as censored

# Create the dataset
survival_data <- data.frame(
  Mouse_ID = rep(mouse_ids, each = 1),
  Time = times,
  Event = events,
  Treatment = rep(treatments, each = n_per_group),
  ID = rep(1:n_total, each = 1)
)

# Add cage information (4 mice per cage, randomized)
cages <- rep(1:(n_total/4), each = 4)
# Shuffle the cages to simulate random cage assignment
set.seed(456)
cage_order <- sample(cages)
survival_data$Cage <- cage_order

# Save the data for later use
write.csv(survival_data, "temp/survival_synthetic_test_data.csv", row.names = FALSE)

# Print sample of the data
print(head(survival_data))

# Run survival statistics
results <- survival_statistics(
  df = survival_data,
  time_column = "Time",
  censor_column = "Event",
  treatment_column = "Treatment",
  id_column = "Mouse_ID",
  cage_column = "Cage",
  reference_group = "Control"
)

# Print results
print(results$results)

# Plot results
if (!is.null(results$forest_plot)) {
  pdf("temp/forest_plot_test.pdf", width = 10, height = 6)
  print(results$forest_plot)
  dev.off()
}

if (!is.null(results$km_plot)) {
  pdf("temp/km_plot_test.pdf", width = 10, height = 6)
  print(results$km_plot)
  dev.off()
}

# Print method used
cat("\nMethod Used:", results$method_used, "\n")

# Summary of event counts by treatment
event_summary <- tapply(survival_data$Event, survival_data$Treatment, function(x) {
  c(events = sum(x), total = length(x), percentage = sum(x)/length(x)*100)
})
print(event_summary) 