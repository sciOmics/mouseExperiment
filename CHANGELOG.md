# Changelog

All notable changes to the mouseExperiment package will be documented in this file.

## [Unreleased]

### Added
- New function `plot_combination_index()` to visualize Combination Index over time with synergy indicators
- Implemented `post_power_analysis()` function for calculating statistical power and sample size recommendations based on experimental data

### Fixed
- Fixed `post_power_analysis()` function to properly handle AUC analysis with better error checking and data validation
- Implemented simulation-based power analysis method in `post_power_analysis()` function
- Fixed AUC calculation in `post_power_analysis()` to properly validate data and handle treatments with small sample sizes
- Fixed contrasts issue in `post_power_analysis()` for AUC method by adding proper validation of treatment groups
- Fixed `plot_forest()` to properly handle both `HR` and `Hazard_Ratio` column naming from `survival_statistics()` results
- Fixed `plot_forest()` to properly display confidence intervals and handle special treatment group names like "HDACi + PD1"
- Improved error handling in `plot_forest()` to handle missing hazard ratios and confidence intervals
- Improved compatibility between column naming conventions in various plotting functions
- Fixed Events/Total count calculation in `survival_statistics()` to correctly count unique subjects and their events per treatment group
- Fixed hazard ratio calculation in `survival_statistics()` to correctly set the reference group and handle coefficient names
- Fixed `survival_statistics()` to properly return the `method_used` object in the results
- Fixed `print_results()` function to properly handle survival statistics summary table

### Changed
- Renamed `forest_plot()` function to `plot_forest()` for more consistent function naming throughout the package
- Removed `forest_plot` object from `survival_statistics()` output as this functionality is now handled by the separate `plot_forest()` function
- Moved `plot_forest()` function (previously `forest_plot()`) from `plot_caterpillar.R` to its own file `plot_forest.R` for better code organization
- Enhanced summary output in `tumor_growth_statistics()` to provide detailed description of statistical tests and methods used
- Improved growth rate calculation description in result summaries to clarify the log-transformation and interpretation
- Modified posthoc tests in `tumor_growth_statistics()` for AUC analysis (model_type = "auc") to use Welch's t-tests instead of standard t-tests, which better handles unequal variances between treatment groups

## [0.2.3] - 2025-03-27

### Added
- Enhanced growth rate calculation description in the tumor_growth_statistics function summary
- Added plot_bliss function for visualizing tumor growth inhibition and Bliss synergy over time
- Added plot_combination_index function for visualizing Combination Index values with colored background regions indicating synergy, additivity, and antagonism

### Changed
- Fixed event and total counts calculation in survival_statistics
- Improved implementation of print_results to properly handle survival statistics summary table
- Improved integration between tumor_growth_statistics and plot_auc functions
- Enhanced summary output in tumor_growth_statistics to provide detailed description of statistical tests and methods used
- Removed debugging output from the survival_statistics function for cleaner code
- Modified plot_survival to use default risk table styling from survminer for better readability

### Fixed
- Fixed bug in print_results function that was incorrectly indexing the summary table

## [0.1.0] - 2025-03-27

### Added
- Initial release of the mouseExperiment package
- Added `plot_auc` function to create visualizations for Area Under the Curve (AUC) data
- Added `tumor_growth_statistics` function for analyzing tumor growth data
- Properly documented all datasets with roxygen2 comments
- Added `plot_growth_rate` function to visualize tumor growth rates from tumor_growth_statistics output
- Functions for tumor growth statistics
- Functions for survival analysis
- Functions for combination index calculation
- Plotting functions for AUC, growth rates, and survival curves
- Example synthetic datasets

### Fixed
- Fixed issue in `plot_auc` function by replacing `annotate_figure` with a simpler approach using `ggpubr::ggarrange`
- Improved error handling in `tumor_growth_statistics` function
- Fixed model type argument matching to properly handle "lme4" and "auc" types
- Fixed dataset exports by adding proper documentation and using LazyData: true
- Added examples for all datasets
- Fixed issue with `tumor_growth_statistics` function returning NULL for model object, ANOVA results, and post-hoc tests
- Fixed posthoc comparison error in `tumor_growth_statistics` function by changing dynamic variable reference to direct variable name
- Fixed AUC calculation in `tumor_growth_statistics` to use original untransformed volume data instead of log-transformed data
- Fixed AUC calculation to use composite IDs that combine subject ID and treatment group, ensuring correct AUC values for each unique subject-treatment combination
- Fixed `plot_auc` function error where `is.named()` function was called but not defined
- Fixed event counting bug in `survival_statistics` function where event totals were incorrectly calculated when generating the formatted summary table
- Fixed hazard ratio calculation in `survival_statistics` function to properly count unique mice rather than summing across all timepoints, ensuring correct Events/Total counts
- Fixed display issues in `plot_survival` function's risk table by removing colored text from group names while preserving colors in the main plot and legend

### Changed
- Simplified `plot_auc` function to focus on a single scatterplot showing individual data points by treatment group
- Added support in `plot_auc` for differentiating between extrapolated and non-extrapolated data points (open vs. filled circles)
- Added configurable error bars to `plot_auc` (none, SEM, SD, or 95% CI)
- Added option to show/hide mean lines in `plot_auc`
- Removed dependency on ggpubr for `plot_auc` function
- Improved integration between `tumor_growth_statistics` and `plot_auc`:
  - Updated `calculate_auc` function to directly accept time and volume vectors
  - Modified AUC result structure for seamless integration with the new plot_auc function
- Enhanced `tumor_growth_statistics` function:
  - Improved summary output to provide detailed description of statistical tests and methods used instead of duplicating ANOVA results
  - Added detailed explanation of growth rate calculation in summary output, including mathematical interpretation
- Enhanced visualization in `plot_growth_rate` and `plot_auc` functions by adding x and y-axis lines for improved readability
- Improved plot aesthetics in `plot_growth_rate` and `plot_auc` functions by replacing theme_minimal with theme_classic to remove the grid background