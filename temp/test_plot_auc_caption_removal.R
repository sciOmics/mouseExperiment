# Test script to verify that the default caption has been removed from the plot_auc function

# Load required libraries
library(mouseExperiment)
library(ggplot2)

# Create some sample data
set.seed(123)
auc_data <- data.frame(
  ID = paste0("Mouse_", 1:20),
  Group = rep(c("Control", "Treatment"), each = 10),
  AUC = c(rnorm(10, 100, 15), rnorm(10, 70, 10)),
  Extrapolated = c(rep(FALSE, 15), rep(TRUE, 5))
)

# Test 1: Plot with extrapolation but no explicit caption
# The default caption "Open circles represent extrapolated values" should not appear
plot1 <- plot_auc(
  auc_data = auc_data,
  extrapolated_column = "Extrapolated",
  show_mean = TRUE,
  error_bar_type = "SEM"
)

# Save the plot to a PDF for inspection
pdf("temp/plot_auc_no_caption.pdf", width = 8, height = 6)
print(plot1)
dev.off()

# Test 2: Plot with extrapolation and an explicit caption
# The provided caption should appear
plot2 <- plot_auc(
  auc_data = auc_data,
  extrapolated_column = "Extrapolated",
  show_mean = TRUE,
  error_bar_type = "SEM",
  caption = "This is a custom caption"
)

# Save the plot to a PDF for inspection
pdf("temp/plot_auc_custom_caption.pdf", width = 8, height = 6)
print(plot2)
dev.off()

cat("Test plots saved to:\n")
cat("  - temp/plot_auc_no_caption.pdf (should have no caption)\n")
cat("  - temp/plot_auc_custom_caption.pdf (should show custom caption)\n") 