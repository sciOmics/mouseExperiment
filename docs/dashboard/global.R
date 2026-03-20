# Global setup for mouseExperimentDashboard
# This file is sourced once when the Shiny app starts

# Check if mouseExperiment is available
if (!requireNamespace("mouseExperiment", quietly = TRUE)) {
  stop(
    "mouseExperiment package is required but not installed.\n",
    "Load it from local path first:\n",
    "  devtools::load_all('../mouseExperiment')\n",
    "Or install it with:\n",
    "  devtools::install_local('../mouseExperiment')"
  )
}

# Load required libraries
library(shiny)
library(shinycssloaders)
library(shinyFeedback)
library(shinyjs)
library(DT)
library(dplyr)
library(ggplot2)
library(scales)
library(mouseExperiment)

# Source helper functions
source("R/helpers.R")

# Application constants
APP_VERSION <- "0.1.0 (Development)"
MAX_FILE_SIZE <- 50 * 1024^2  # 50 MB
SUPPORTED_FORMATS <- c("csv", "rda", "RData")

# Color scheme for plots
DASHBOARD_COLORS <- list(
  primary = "#0073B7",
  success = "#00A65A",
  info = "#00BCD4",
  warning = "#FF9800",
  danger = "#DD4B39"
)

# Sample datasets available to users
SAMPLE_DATASETS <- list(
  "combo_treatment" = list(
    name = "Combination Treatment Study",
    package = "mouseExperiment",
    description = "Drug combination efficacy with cage effects"
  ),
  "dose_levels" = list(
    name = "Dose-Response Study",
    package = "mouseExperiment",
    description = "Dose escalation study with multiple treatment levels"
  )
)

# Numeric precision for display
DECIMAL_PLACES <- 3

# Enable shinyjs
shinyjs::useShinyjs()
