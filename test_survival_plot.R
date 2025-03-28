# Test script for plot_survival function
library(mouseExperiment)
library(ggplot2)
library(survival)
library(survminer)

# Load the synthetic data
data(combo_treatment_synthetic_data)

# Add survival data for testing
set.seed(123)
survival_data <- combo_treatment_synthetic_data[!duplicated(combo_treatment_synthetic_data$Mouse_ID), ]
survival_data$Survival_Time <- runif(nrow(survival_data), 10, 60)
survival_data$Survival_Censor <- sample(c(0, 1), nrow(survival_data), replace = TRUE, prob = c(0.3, 0.7))

# Create survival plot with risk table
p1 <- plot_survival(
  survival_data,
  time_column = "Survival_Time",
  censor_column = "Survival_Censor",
  treatment_column = "Treatment",
  id_column = "Mouse_ID",
  cage_column = "Cage",
  show_risk_table = TRUE,
  title = "Survival Plot with Risk Table (Fixed)",
  subtitle = "Group names in risk table should be black"
)

# Save the plot to PDF
pdf("survival_plot_fixed.pdf", width = 10, height = 8)
print(p1)
dev.off()

# Print completion message
cat("Test completed successfully. Plot saved as 'survival_plot_fixed.pdf'\n")
cat("Please verify that the risk table group names are black, not colored.\n") 