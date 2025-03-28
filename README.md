# mouseExperiment

An R package for analyzing mouse tumor growth experiments, including survival analysis and synergy calculations.

## Features

- Tumor growth analysis using mixed-effects models or area under the curve (AUC)
- Survival analysis with Kaplan-Meier curves and hazard ratios
- Synergy analysis using Bliss Independence and Combination Index methods
- Comprehensive visualization tools for all analyses
- Support for cage effects and repeated measures
- Synthetic datasets for testing and demonstration

## Installation

```r
# Install from CRAN
install.packages("mouseExperiment")

# Or install the development version from GitHub
# devtools::install_github("yourusername/mouseExperiment")
```

## Usage

### Tumor Growth Analysis

```r
library(mouseExperiment)

# Load example data
data(combo_treatment_synthetic_data)

# Analyze tumor growth
results <- tumor_growth_statistics(
  data = combo_treatment_synthetic_data,
  time_column = "Day",
  volume_column = "Volume",
  id_column = "Mouse_ID",
  treatment_column = "Treatment",
  cage_column = "Cage",
  model_type = "lme4"
)

# Visualize growth rates
plot_growth_rate(results$growth_rates)
```

### Survival Analysis

```r
# Load survival data
data(survival_synthetic_data)

# Perform survival analysis
surv_results <- survival_statistics(
  data = survival_synthetic_data,
  time_column = "Time",
  event_column = "Event",
  treatment_column = "Treatment",
  id_column = "Mouse_ID"
)

# Create survival plot
plot_survival(surv_results$survival_data)
```

### Synergy Analysis

```r
# Calculate Bliss synergy
bliss_results <- bliss_independence(
  data = combo_treatment_synthetic_data,
  time_column = "Day",
  volume_column = "Volume",
  id_column = "Mouse_ID",
  treatment_column = "Treatment",
  cage_column = "Cage"
)

# Visualize Bliss synergy
plot_bliss(bliss_results)

# Calculate Combination Index
synergy_results <- synergy_stats(
  data = combo_treatment_synthetic_data,
  time_column = "Day",
  volume_column = "Volume",
  id_column = "Mouse_ID",
  treatment_column = "Treatment",
  cage_column = "Cage"
)

# Visualize Combination Index
plot_combination_index(synergy_results$synergy_summary)
```

## Documentation

For detailed documentation and examples, see the package vignette:

```r
browseVignettes("mouseExperiment")
```

## Functions

### Tumor Growth Analysis
- `tumor_growth_statistics()`: Analyze tumor growth using mixed-effects models or AUC
- `plot_growth_rate()`: Visualize tumor growth rates
- `plot_auc()`: Visualize area under the curve results

### Survival Analysis
- `survival_statistics()`: Perform survival analysis
- `plot_survival()`: Create Kaplan-Meier survival curves
- `forest_plot()`: Visualize hazard ratios
- `print_results()`: Print formatted results

### Synergy Analysis
- `bliss_independence()`: Calculate Bliss synergy
- `synergy_stats()`: Calculate Combination Index
- `plot_bliss()`: Visualize Bliss synergy
- `plot_combination_index()`: Visualize Combination Index

## Datasets

The package includes several synthetic datasets for testing and demonstration:

- `combo_treatment_synthetic_data`: Tumor growth data with combination treatments
- `survival_synthetic_data`: Survival data with multiple treatment groups
- `single_treatment_synthetic_data`: Tumor growth data with single treatments

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use this package in your research, please cite it as:

```r
citation("mouseExperiment")
```