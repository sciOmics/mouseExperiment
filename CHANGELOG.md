# Changelog

All notable changes to the mouseExperiment package will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- New function `plot_combination_index()` to visualize Combination Index over time with synergy indicators
- Implemented `post_power_analysis()` function for calculating statistical power and sample size recommendations based on experimental data
- Added `Treatment` column to the power_analysis object returned by `post_power_analysis` function to clearly identify which treatment group each power estimate refers to
- Enhanced power curve visualization in `post_power_analysis` to distinguish between treatment groups using colors

### Fixed
- Fixed `survival_statistics()` function to correctly generate Kaplan-Meier plots, resolving the "object of type 'symbol' is not subsettable" error
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
- Fixed error in `survival_statistics()` when calling `print_results()` by correctly passing necessary parameters
- Fixed `tumor_growth_statistics()` to correctly report the number of subjects in summary's data_description by using composite IDs that include cage information
- Fixed `plot_growth_rate` function to correctly identify and display all mice (8 per group instead of 4) in treatment groups by incorporating cage information
- Fixed `plot_auc` function to correctly display extrapolated points and ensure mean bars are shown in the correct columns
- Fixed `plot_auc` function to properly distinguish between extrapolated and non-extrapolated data points by improving shape mapping
- Changed default value of `show_mean` parameter in `plot_auc` function to TRUE
- Modified mean display in `plot_auc` to use horizontal bars only (no diamonds)

### Changed
- Renamed `forest_plot()` function to `plot_forest()` for more consistent function naming throughout the package
- Removed `forest_plot` object from `survival_statistics()` output as this functionality is now handled by the separate `plot_forest()` function
- Moved `plot_forest()` function (previously `forest_plot()`) from `plot_caterpillar.R` to its own file `plot_forest.R` for better code organization
- Enhanced summary output in `tumor_growth_statistics()` to provide detailed description of statistical tests and methods used
- Improved growth rate calculation description in result summaries to clarify the log-transformation and interpretation
- Modified posthoc tests in `tumor_growth_statistics()` for AUC analysis (model_type = "auc") to use Welch's t-tests instead of standard t-tests, which better handles unequal variances between treatment groups
- Restored `colors`, `point_size`, and `jitter_width` parameters to the `plot_auc()` function for greater customization of visualization output
- Removed boxplots from `plot_auc()` function to restore original functionality that focused on points with optional mean and error bars

## [0.3.0] - 2023-07-15

### Fixed
- Fixed `tumor_growth_statistics()` to properly account for all mice in each treatment group when mice have the same ID but are in different cages, by including the cage identifier in the composite ID creation
- Fixed `tumor_auc_analysis()` and `calculate_auc_values()` functions to properly account for all mice in each treatment group when mice have the same ID but are in different cages

### Added
- Added a cage_column parameter to `tumor_auc_analysis()` and `calculate_auc_values()` functions

### Changed
- Enhanced summary output in `tumor_growth_statistics()` to provide detailed description of statistical tests and methods used
- Modified posthoc tests in `tumor_growth_statistics()` for AUC analysis (model_type = "auc") to use Welch's t-tests instead of standard t-tests, which better handles unequal variances between treatment groups

## [0.2.0] - 2023-06-01

### Changed
- Enhanced growth rate calculation description in the tumor_growth_statistics function summary
- Added ability to customize plot aesthetics in plot_tumor_growth function
- Added support for additional statistical methods in tumor_growth_statistics
- Improved integration between tumor_growth_statistics and plot_auc functions
- Enhanced summary output in tumor_growth_statistics to provide detailed description of statistical tests and methods used

### Fixed
- Fixed an issue with extrapolation in calculate_dates function
- Corrected error handling in tumor_auc_analysis for single-treatment datasets

## [0.1.0] - 2023-05-01

### Added
- Initial release of the mouseExperiment package
- Added `calculate_volume()` function for calculating tumor volume from measurements
- Added `calculate_dates()` function for converting study days to calendar dates
- Added `tumor_growth_statistics` function for analyzing tumor growth data
- Added `tumor_auc_analysis()` function for area under the curve analysis
- Added visualization functions including `plot_tumor_growth()`, `plot_survival()`, `plot_auc()`, and `forest_plot()`
- Added `bliss_interaction()` and `plot_bliss()` functions for drug interaction analysis
- Added `post_power_analysis()` function for statistical power calculations
- Included sample datasets for testing and demonstration

### Changed
- Changed posthoc tests for AUC analysis in `tumor_growth_statistics()` to use Welch's t-tests instead of standard t-tests to handle unequal variances
- Updated `tumor_growth_statistics()` function to include an `Extrapolated` column in AUC analysis results for seamless integration with `plot_auc()`
- Enhanced `plot_auc()` function to show extrapolated data points as open circles and automatically detect extrapolation status from the data