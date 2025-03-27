# mouseExperiment

An R package for analyzing and visualizing data from mouse tumor growth experiments.

## Features

This package provides tools to:

- Calculate tumor volumes from measurement data
- Analyze tumor growth kinetics with mixed-effects models
- Perform survival analysis with hazard ratios
- Visualize tumor growth with customizable plots
- Analyze dose-response relationships
- Investigate drug synergy effects
- Conduct statistical analyses
- Perform post-hoc power analysis
- Generate Area Under the Curve (AUC) analysis and visualization

## Installation

This package is for internal use by Insight BioAnalytics personnel.

```r
# Install from local file
install.packages("path/to/mouseExperiment_x.y.z.tar.gz", repos = NULL, type = "source")

# Or with devtools
devtools::install_local("path/to/mouseExperiment")
```

## Key Functions

### Data Processing
- `calculate_volume()` - Calculate tumor volume from caliper measurements
- `calculate_dates()` - Convert dates to experimental days
- `tumor_auc_analysis()` - Calculate and analyze Area Under the Curve for tumor growth data

### Visualization
- `plot_tumor_growth()` - Create tumor growth curves with optional treatment day indicators
- `plot_survival()` - Generate Kaplan-Meier survival curves
- `plot_caterpillar()` - Create coefficient plots for mixed-effects models
- `plot_treatments()` - Visualize treatment administration schedules
- `plot_auc()` - Create comprehensive visualizations of Area Under the Curve (AUC) data

### Statistical Analysis
- `tumor_growth_statistics()` - Analyze tumor growth with mixed-effects models
- `survival_statistics()` - Calculate hazard ratios with confidence intervals
- `dose_response_statistics()` - Analyze dose-response relationships
- `analyze_drug_synergy()` - Test for synergistic effects between treatments
- `analyze_drug_synergy_over_time()` - Examine synergy effects across timepoints
- `post_power_analysis()` - Calculate achieved power from experimental results

## Example Usage

```r
library(mouseExperiment)

# Load demo data
data(combo_treatment_synthetic_data)

# Calculate tumor volumes
data <- calculate_volume(combo_treatment_synthetic_data)

# Calculate experimental days from calendar dates
data <- calculate_dates(data, start_date = "03/24/2025", date_format = "%m/%d/%Y")

# Plot tumor growth with treatment indicators
growth_plot <- plot_tumor_growth(data, treatment_days = c(0, 7, 14, 21))

# Visualize treatment schedule
data(combo_treatment_schedule)
treatment_plot <- plot_treatments(combo_treatment_schedule)

# Run tumor growth statistics
results <- tumor_growth_statistics(data)
print(results$anova)

# Plot survival curves
plot_survival(data, show_legend = TRUE, legend_title = "Treatment Group")

# Calculate hazard ratios
survival_results <- survival_statistics(data, reference_group = "Control")
print(survival_results$results)
survival_results$forest_plot

# Analyze dose-response relationship
data(dose_levels_synthetic_data)
dose_data <- calculate_volume(dose_levels_synthetic_data)
dose_data <- calculate_dates(dose_data, start_date = "24-Mar", date_format = "%d-%b", year = 2025)
dose_results <- dose_response_statistics(dose_data, dose_column = "Dose")
dose_results$plots$scatter

# Investigate drug synergy
synergy_results <- analyze_drug_synergy(data,
                                     group_column = "Treatment",
                                     volume_column = "Volume",
                                     drug_a = "Drug A",
                                     drug_b = "Drug B", 
                                     combo = "Combo",
                                     control = "Control")
print(synergy_results$combo_index)

# Perform post-hoc power analysis
power_results <- post_power_analysis(data, 
                                  effect_size = 0.8,
                                  power_target = 0.8)
print(power_results$achieved_power)

# Calculate AUC and analyze results
auc_results <- tumor_auc_analysis(data)
print(auc_results$auc_summary)
print(auc_results$auc_comparisons)

# Create custom AUC visualizations
example_auc_data <- data.frame(
  ID = paste0("Mouse", 1:20),
  Group = rep(c("Control", "Treatment A", "Treatment B", "Treatment C"), each=5),
  AUC = c(runif(5, 10, 15), runif(5, 8, 12), runif(5, 5, 10), runif(5, 3, 8))
)
auc_plot <- plot_auc(example_auc_data)
```

## Statistical Methodology

### Tumor Growth Statistics

The `tumor_growth_statistics()` function uses linear mixed-effects models to analyze tumor growth data:

1. **Model specification**: The default model is `log(Volume) ~ Day * Group + (1|ID)`, which:
   - Uses log-transformed volume to account for exponential growth
   - Tests for effects of time (Day), treatment (Group), and their interaction
   - Includes random intercepts for individual subjects (mice)

2. **Statistical tests**:
   - Type III ANOVA for fixed effects
   - Pairwise comparisons with Bonferroni correction
   - Diagnostic plots for model validation

3. **Interpretation**:
   - Significant Day effect: Tumor volume changes over time
   - Significant Group effect: Treatment affects overall tumor volume
   - Significant interaction: Treatment affects growth rate

### Survival Statistics

The `survival_statistics()` function uses Cox proportional hazards models:

1. **Model specification**: The default model is `Surv(time, status) ~ Group`, which:
   - Models time-to-event data with censoring
   - Compares survival between treatment groups

2. **Statistical outputs**:
   - Hazard ratios with confidence intervals
   - p-values for each comparison
   - Forest plot visualization

3. **Interpretation**:
   - Hazard ratio < 1: Reduced risk compared to reference group
   - Hazard ratio > 1: Increased risk compared to reference group
   - Confidence intervals not crossing 1 indicate statistical significance

### Dose-Response Analysis

The `dose_response_statistics()` function uses several modeling approaches:

1. **Linear regression**: Tests for a linear relationship between dose and response

2. **ANOVA with post-hoc tests**: Compares responses between different dose levels

3. **Polynomial trend analysis**: Tests for non-linear responses (quadratic, cubic)

4. **Non-linear modeling**: Fits 4-parameter logistic models when appropriate

### AUC Analysis

The `tumor_auc_analysis()` function calculates and analyzes Area Under the Curve:

1. **AUC calculation**: Uses trapezoidal method to calculate area under tumor growth curve

2. **Statistical analysis**: Compares AUC values between treatment groups using ANOVA

3. **Interpretation**:
   - Lower AUC indicates more effective tumor control
   - Provides a single value summary of treatment effect over entire experiment

### Drug Synergy Analysis

The `analyze_drug_synergy()` function evaluates combination treatments:

1. **Additivity model**: Calculates expected combined effect assuming independence

2. **Synergy testing**: Compares observed combination effect to expected effect

3. **Interpretation**:
   - Observed > Expected: Synergistic effect
   - Observed = Expected: Additive effect
   - Observed < Expected: Antagonistic effect

## Dependencies

This package requires several R packages:
- `ggplot2` and `ggpubr` for plotting
- `survival` and `survminer` for survival analysis
- `lme4`, `lmerTest`, and `emmeans` for statistical modeling
- `drc` for dose-response modeling
- `dplyr` for data manipulation
- `anytime` for date parsing
- `rlang` for tidy evaluation

## Documentation

For more detailed documentation and examples, see the vignettes:

```r
browseVignettes("mouseExperiment")
```

## Recent Improvements

- **New Features**: Added AUC analysis with `tumor_auc_analysis()` and `plot_auc()` functions
- **Code organization**: Improved error handling for `tumor_growth_statistics()`
- **Documentation**: Enhanced package documentation with roxygen2 and properly documented all datasets
- **Reliability**: Fixed model type argument handling in statistical functions
- **Visualization**: Implemented improved plotting for AUC data with multiple visualization options
- **Datasets**: Added comprehensive examples to all dataset documentation
- **Package Structure**: Added LazyData support for easier access to included datasets

## License

This package is proprietary software owned by Insight BioAnalytics. All rights reserved.
Not for distribution or use without explicit permission.