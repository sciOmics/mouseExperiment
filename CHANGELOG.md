# Changelog

All notable changes to the mouseExperiment package will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- New function `plot_combination_index()` to visualize Combination Index over time with synergy indicators
- Implemented `post_power_analysis()` function for calculating statistical power and sample size recommendations based on experimental data

### Fixed
- Fixed `post_power_analysis()` function to properly handle AUC analysis with better error checking and data validation
- Fixed contrasts error in `post_power_analysis()` AUC method by ensuring treatment groups are properly converted to factors
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
- Fixed median survival calculation in `survival_statistics()` to correctly handle cases where more than 50% of subjects have events, preventing incorrect "Not Reached" messages when median survival should be calculable
- Fixed reporting of median survival times for groups with >50% events, now properly calculating and displaying median survival times rather than showing "Not calculated" message

### Changed
- Renamed `forest_plot()` function to `plot_forest()` for more consistent function naming throughout the package
- Removed `forest_plot` object from `survival_statistics()` output as this functionality is now handled by the separate `plot_forest()` function
- Moved `plot_forest()` function (previously `forest_plot()`) from `plot_caterpillar.R` to its own file `plot_forest.R` for better code organization
- Enhanced summary output in `tumor_growth_statistics()` to provide detailed description of statistical tests and methods used
- Improved growth rate calculation description in result summaries to clarify the log-transformation and interpretation
- Modified posthoc tests in `tumor_growth_statistics()` for AUC analysis (model_type = "auc") to use Welch's t-tests instead of standard t-tests, which better handles unequal variances between treatment groups
- Restored `colors`, `point_size`, and `jitter_width` parameters to the `plot_auc()` function for greater customization of visualization output

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
- Initial release with core functionality
- Added `calculate_volume` function to transform length and width measurements into tumor volume
- Added `calculate_dates` function for converting experimental days into calendar dates
- Added `tumor_growth_statistics` function for analyzing tumor growth data
- Added `survival_statistics` function for analyzing survival data
- Added visualization functions including `plot_tumor_growth()`, `plot_survival()`, `plot_auc()`, and `forest_plot()`
- Added `post_power_analysis` function for estimating statistical power
- Added synthetic datasets for examples and testing

### Changed
- Changed hazard ratio calculation in `survival_statistics()` to properly handle the reference group
- Changed column naming from EventsTotal to Events.Total for consistency with other column names
- Changed posthoc tests for AUC analysis in `tumor_growth_statistics()` to use Welch's t-tests instead of standard t-tests to handle unequal variances
- Enhanced `plot_auc()` function to show extrapolated data points as open circles and automatically detect extrapolation status from the data
- Updated `tumor_growth_statistics()` function to include an `Extrapolated` column in AUC analysis results for seamless integration with `plot_auc()`
- Improved `tumor_auc_analysis()` function to implement linear extrapolation using the last N points specified by `extrapolation_points`

### Fixed
- Fixed median survival calculation in `survival_statistics()` to correctly handle cases with more than 50% censored data