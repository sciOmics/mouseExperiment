# mouseExperimentDashboard

An interactive web-based dashboard for the `mouseExperiment` R package. Provides a point-and-click interface for tumor growth analysis without requiring R programming knowledge.

## Installation

```r
# Step 1: Load mouseExperiment package from local path
devtools::load_all("../mouseExperiment")

# Step 2: Install dashboard dependencies
install.packages(c("shiny", "shinycssloaders", "shinyFeedback", "shinyjs", "DT",
                   "dplyr", "ggplot2", "scales", "tidyr"))

# Step 3: Load the dashboard package
devtools::load_all()
```

## Usage

```r
# Run the dashboard
mouseExperimentDashboard::runDashboard()
```

The dashboard will open in your default web browser at `http://localhost:3838`

## Features

**Current (MVP):**
- Data upload (CSV, RDA formats)
- Automatic column detection and mapping
- Data preview and validation
- Tumor growth statistical analysis (LME4 or AUC methods)
- Interactive visualizations
- Download results (plots, tables)

**Coming Soon:**
- Survival analysis
- Synergy analysis
- Power analysis
- Report generation

## Requirements

- R >= 4.0
- mouseExperiment >= 0.3.0
- Shiny >= 1.7.0

## Development

```r
# Load the package in development mode
devtools::load_all()

# Run the app
runDashboard()

# Run tests
devtools::test()
```

## License

MIT License - See LICENSE file for details
