# mouseExperiment

An R package for analyzing and visualizing mouse tumor growth experiments. This package offers comprehensive tools for processing, analyzing, and visualizing tumor growth data from preclinical mouse experiments, including both frequentist and Bayesian statistical approaches.

## Features

- Calculate tumor volume from length and width measurements
- Track experiment timeline with date calculations
- Plot tumor growth curves by treatment group
- Generate survival plots with Kaplan-Meier curves
- Perform statistical analysis of tumor growth data using linear mixed-effects models
- Support for Bayesian statistical modeling of tumor growth
- Analyze survival data using Cox proportional hazards models

## Installation

```r
# Install from GitHub
devtools::install_github("https://github.com/SciOmicsLab/mouseExperiment")
```

## Dependencies

This package requires several R packages:
- **ggplot2** and **ggpubr** for plotting
- **survival** and **survminer** for survival analysis
- **lme4**, **lmerTest**, and **emmeans** for statistical modeling
- **brms** and **bridgesampling** for Bayesian analysis
- **anytime** for date handling
- **performance** for model diagnostics
- **future** for parallel computing

## Detailed Function Documentation

### Data Processing Functions

#### `calculate_volume()`
**Purpose**: Calculates tumor volume using different formulas common in preclinical oncology research.

**Arguments**:
- `df`: A data frame containing tumor measurements
- `length_column`: The name of the column storing tumor length measurements (default: "Length")
- `width_column`: The name of the column storing tumor width measurements (default: "Width")
- `height_column`: The name of the column storing tumor height measurements (optional, for 3D formulas)
- `formula`: The formula to use for volume calculation (default: "ellipsoid")
- `in_place`: Logical, whether to modify the input data frame (TRUE) or return a new data frame (FALSE, default)

**Output**: Returns a data frame with a "Volume" column appended, calculated using the specified formula

**Supported Formulas**:
1. **Ellipsoid** (default): V = (length × width² × π) / 6
   - Most accurate for ovoid tumors when true height isn't measured
   - Widely used in xenograft studies

2. **Modified Ellipsoid**: V = (length × width²) / 2
   - Simplified formula omitting π/3
   - Used in many older studies
   - The formula previously used as default in this package

3. **Ellipsoid with 3 axes**: V = (length × width × height × π) / 6
   - Most accurate when all three dimensions can be measured
   - Requires height measurement or uses width as approximation

4. **Cylinder**: V = (π × width² × length) / 4
   - Used for more cylindrical tumors

5. **Sphere**: V = (π × width³) / 6
   - For spherical tumors where only diameter is measured

6. **Box**: V = length × width × height
   - Simple approximation when tumors have more rectangular shape
   - Requires height measurement or uses width as approximation

**Note**: The function ensures that the longer dimension is used as "length" and the shorter as "width", regardless of column labeling. When 3D formulas are used and the height is not provided, the width is used as an approximation for height (a common practice when height cannot be measured).

**Usage**:
```r
# Default ellipsoid formula
df <- calculate_volume(synthetic_data)

# Using the modified ellipsoid formula
df <- calculate_volume(synthetic_data, formula = "modified_ellipsoid")

# Using a 3D formula with height measurements
df <- calculate_volume(synthetic_data, height_column = "Height", formula = "ellipsoid_3axis")
```

#### `calculate_dates()`
**Purpose**: Calculates days from a starting point for each measurement date.

**Arguments**:
- `df`: A data frame containing longitudinal tumor measurements
- `start_date`: Date the tumors were injected
- `date_column`: The name of the column with dates of measurements (default: "Date")
- `in_place`: Logical, whether to modify the input data frame (TRUE) or return a new data frame (FALSE, default)

**Output**: Returns a data frame with a "Day" column appended, representing days since injection.

### Visualization Functions

#### `plot_tumor_growth()`
**Purpose**: Creates a plot of tumor growth over time by treatment group.

**Arguments**:
- `df`: A data frame containing tumor measurements
- `volume_column`: The name of the column storing tumor volume measurements (default: "Volume")
- `day_column`: The name of the column with number of days (default: "Day")
- `group_column`: The name of the column with the group indicator (default: "Group")
- `ID_column`: The name of the column with the individual mouse identifier (default: "ID")
- `group_summary_line`: Boolean, whether to include group average lines (default: TRUE)

**Output**: A ggplot object showing individual tumor growth trajectories colored by treatment group, with optional group means.

#### `plot_survival()`
**Purpose**: Generates Kaplan-Meier survival curves for different treatment groups.

**Arguments**:
- `df`: A data frame containing survival data
- `time_column`: The name of the column representing time to event (default: "Day")
- `censor_column`: The name of the column indicating censoring (1 = event occurred, 0 = censored) (default: "Survival_Censor")
- `group_column`: The name of the column representing different groups for comparison (default: "Group")

**Output**: A Kaplan-Meier survival plot with risk tables, showing survival probability over time for each treatment group.

**Statistical Details**: 
- Uses the Kaplan-Meier method to estimate survival function
- Includes log-rank test p-value to assess differences between groups
- Displays number of subjects at risk in a table below the main plot

#### `plot_caterpillar()`
**Purpose**: Creates caterpillar plots (coefficient plots) for frequentist mixed-effects models.

**Arguments**:
- `model`: A linear mixed-effects model object from lme4::lmer() or a similar function
- `title`: A title for the plot (default: "Fixed Effects with 95% Confidence Intervals")
- `colors`: A vector of colors for different coefficient groups
- `show_intercept`: Logical, whether to include the intercept in the plot (default: TRUE)
- `ci_level`: Confidence level for intervals (default: 0.95 for 95% CI)

**Output**: A ggplot object showing fixed effects coefficients with confidence intervals.

**Details**:
- Extracts coefficients and standard errors from model output
- Groups coefficients by type (intercept, main effects, interactions)
- Shows 95% confidence intervals (or other specified level)
- Includes a vertical reference line at zero
- Color-codes different coefficient types

**Usage**:
```r
# After running tumor growth statistics
results <- tumor_growth_statistics(df)

# Create and display caterpillar plot
cat_plot <- plot_caterpillar(results$model)
print(cat_plot)
```

#### `plot_caterpillar_bayes()`
**Purpose**: Creates caterpillar plots for Bayesian models from the brms package.

**Arguments**:
- `model`: A brmsfit object from brms::brm()
- `title`: A title for the plot (default: "Posterior Distributions with 95% Credible Intervals")
- `colors`: A vector of colors for different parameter groups
- `show_intercept`: Logical, whether to include the intercept in the plot (default: TRUE)
- `ci_level`: Credible interval level (default: 0.95 for 95% CI)
- `show_random`: Logical, whether to show group-level (random) effects (default: FALSE)

**Output**: A ggplot object showing posterior distributions of model parameters with credible intervals.

**Details**:
- Extracts posterior medians and credible intervals from Bayesian model
- Groups parameters by type (intercept, main effects, interactions, random effects)
- Shows 95% credible intervals (or other specified level)
- Includes a vertical reference line at zero
- Option to include random effects (can be useful but may make plots very large)

**Usage**:
```r
# After running Bayesian analysis
bayes_results <- tumor_growth_statistics_bayes(df)

# Create and display Bayesian caterpillar plot
bayes_cat_plot <- plot_caterpillar_bayes(bayes_results$model)
print(bayes_cat_plot)

# Show random effects too (may be many parameters)
bayes_cat_plot_with_random <- plot_caterpillar_bayes(
  bayes_results$model, 
  show_random = TRUE
)
```

### Statistical Analysis Functions

#### `tumor_growth_statistics()`
**Purpose**: Performs statistical analysis of tumor growth data using mixed-effects models.

**Arguments**:
- `df`: A data frame containing tumor growth data
- `time_column`: The name of the column for time points (default: "Day")
- `volume_column`: The name of the column for tumor volume (default: "Volume")
- `group_column`: The name of the column for treatment groups (default: "Group")
- `id_column`: The name of the column for individual subjects (default: "ID")

**Output**: A list containing:
- `model`: The fitted linear mixed-effects model (using log-transformed volume)
- `anova`: Type III ANOVA table for fixed effects
- `posthoc`: Pairwise comparisons between groups (with Bonferroni correction)
- `plots`: Diagnostic and tumor growth visualization plots

**Statistical Details**:
- Applies log-transformation to volume data to handle exponential growth patterns
- Uses a linear mixed-effects model with random intercepts for each mouse
- Tests for time, treatment, and interaction effects using Type III ANOVA
- Performs post-hoc comparisons with Bonferroni adjustment for multiple testing
- Provides diagnostic plots to check model assumptions (normality, homoscedasticity)

**Interpretation**:
- The ANOVA table shows significance of main effects and interactions
- Post-hoc tests reveal specific group differences at different time points
- The log-transformation helps meet the assumption of normally distributed residuals

#### `survival_statistics()`
**Purpose**: Analyzes survival data using Cox proportional hazards models.

**Arguments**:
- `df`: A data frame containing the data to be analyzed
- `time_column`: The name of the column for time to event (default: "Day")
- `censor_column`: The name of the column for censoring indicator (default: "Survival_Censor")
- `group_column`: The name of the column for treatment groups (default: "Group")
- `id_column`: The name of the column for individual IDs (default: "ID")
- `reference_group`: The level to use as reference for hazard ratios (default: NULL, uses first alphabetically)

**Output**: A list containing:
- `model`: The fitted Cox proportional hazards model
- `results`: Data frame with hazard ratios, confidence intervals, and p-values
- `forest_plot`: A forest plot visualizing hazard ratios with confidence intervals

**Statistical Details**:
- Implements a Cox proportional hazards model with frailty term for individual mice
- Calculates hazard ratios relative to a reference group
- Provides 95% confidence intervals for hazard ratios
- Creates a forest plot for visual comparison of treatment effects

**Interpretation**:
- Hazard ratios < 1 indicate reduced risk (improved survival) compared to reference group
- Hazard ratios > 1 indicate increased risk (worse survival) compared to reference group
- P-values indicate statistical significance of the difference from reference group
- The forest plot provides visual assessment of treatment effects on survival

#### `tumor_growth_statistics_bayes()`
**Purpose**: Performs Bayesian analysis of tumor growth using multi-core processing.

**Arguments**:
- `df`: A data frame containing tumor growth data
- `time_column`: The name of the column for time points (default: "Day")
- `volume_column`: The name of the column for tumor volume (default: "Volume")
- `group_column`: The name of the column for treatment groups (default: "Group")
- `id_column`: The name of the column for individual subjects (default: "ID")
- `cores`: Number of CPU cores to use for parallel processing (default: 4)

**Output**: A list containing:
- `model`: The fitted Bayesian model object
- `posterior_samples`: Data frame with posterior samples of model parameters
- `bayes_factor`: Computed Bayes factor for model comparison (if available)
- `pp_check`: Posterior predictive check plot (if available)

**Statistical Details**:
- Uses Bayesian hierarchical modeling with the brms package
- Specifies informative priors for fixed and random effects
- Employs MCMC for parameter estimation with multiple chains
- Performs posterior predictive checks
- Calculates Bayes factors for model comparison when possible

**Interpretation**:
- Posterior distributions provide probabilistic estimates of parameters
- Credible intervals offer range of plausible parameter values
- Posterior predictive checks assess model fit
- Bayes factors quantify evidence for one model over another

## Usage Examples

```r
library(mouseExperiment)

# Load example data
data(synthetic_data)

# Calculate tumor volume
df <- calculate_volume(synthetic_data)

# Add experiment timeline 
df <- calculate_dates(df, start_date = "2022-02-24")

# Visualize tumor growth
tumor_plot <- plot_tumor_growth(df)
print(tumor_plot)

# Analyze tumor growth with mixed-effects model (log-transformed)
growth_stats <- tumor_growth_statistics(df)
print(growth_stats$anova)  # View ANOVA results
print(growth_stats$posthoc)  # View pairwise comparisons

# Visualize and analyze survival data
surv_plot <- plot_survival(df)
print(surv_plot)

# Get hazard ratios and forest plot
surv_stats <- survival_statistics(df, reference_group = "A")
print(surv_stats$results)

# Optional: Bayesian analysis for more complex modeling
# Note: This requires the brms package and may take longer to run
bayes_results <- tumor_growth_statistics_bayes(df, cores = 4)
```

## Experimental Design Considerations

This package is designed to work with common experimental designs in oncology research:
- **Treatment groups**: Different experimental conditions (e.g., control, drug treatments)
- **Longitudinal measurements**: Repeated measurements of the same mice over time
- **Endpoint data**: Survival time and censoring information

The statistical methods account for:
- The non-independence of repeated measurements from the same subject
- The exponential nature of tumor growth (via log transformation)
- Right-censored survival data (via Kaplan-Meier estimation and Cox models)
- Between-subject variability (via random effects)

## Contributions

Contributions are welcome. Please feel free to submit a Pull Request.

## License

This package is licensed under the MIT license.
