library(mouseExperiment)

# Create simple test data
set.seed(123)
n_per_group <- 8
treatments <- c("Control", "PD1", "HDACi", "HDACi + PD1")

# Generate survival data with good separation to use standard Cox model
test_data <- data.frame(
  ID = 1:(n_per_group * length(treatments)),
  Treatment = rep(treatments, each = n_per_group),
  Day = sample(30:60, n_per_group * length(treatments), replace = TRUE),
  Survival_Censor = c(
    rep(1, n_per_group * 2),  # All Control and PD1 have events
    sample(c(0, 1), n_per_group, replace = TRUE, prob = c(0.7, 0.3)),  # Mostly censored for HDACi
    sample(c(0, 1), n_per_group, replace = TRUE, prob = c(0.8, 0.2))   # Mostly censored for HDACi + PD1
  ),
  Cage = rep(1:8, length(treatments))
)

# Run survival analysis
results <- survival_statistics(
  df = test_data,
  reference_group = "Control"
)

# Check output structure
cat("Output structure:\n")
cat("----------------\n")
cat(paste(names(results), collapse = ", "), "\n\n")

# Check method_used
cat("Method used:", results$method_used, "\n\n")

# Verify forest_plot is not in output
if ("forest_plot" %in% names(results)) {
  cat("Warning: forest_plot is still in the output\n")
} else {
  cat("Success: forest_plot is no longer in the output\n")
}

# Test creating forest plot manually
cat("\nCreating forest plot manually:\n")
forest_plot_obj <- forest_plot(results$results)
cat("Forest plot created successfully\n")

# Save the forest plot to a file
pdf("temp/manual_forest_plot.pdf", width = 10, height = 6)
print(forest_plot_obj)
invisible(dev.off())
cat("Forest plot saved to temp/manual_forest_plot.pdf\n") 