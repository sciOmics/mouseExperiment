# Changelog

All notable changes to the mouseExperiment package will be documented in this file.

## [0.1.0] - 2024-03-27

### Added
- Initial release of the mouseExperiment package
- Added `plot_auc` function to create visualizations for Area Under the Curve (AUC) data
- Added `tumor_growth_statistics` function for analyzing tumor growth data
- Properly documented all datasets with roxygen2 comments
- Added `plot_growth_rate` function to visualize tumor growth rates from tumor_growth_statistics output

### Fixed
- Fixed issue in `plot_auc` function by replacing `annotate_figure` with a simpler approach using `ggpubr::ggarrange`
- Improved error handling in `tumor_growth_statistics` function
- Fixed model type argument matching to properly handle "lme4" and "auc" types
- Fixed dataset exports by adding proper documentation and using LazyData: true
- Added examples for all datasets 
- Fixed issue with `tumor_growth_statistics` function returning NULL for model object, ANOVA results, and post-hoc tests
- Fixed posthoc comparison error in `tumor_growth_statistics` function by changing dynamic variable reference to direct variable name 
- Fixed AUC calculation in `tumor_growth_statistics` to use original untransformed volume data instead of log-transformed data
- Enhanced integration between `tumor_growth_statistics` and `plot_auc` functions:
  - Added Group column to AUC data for compatibility with plot_auc
  - Added group_order parameter to plot_auc to control the order of treatment groups in plots 
  - Improved implementation of "auc" model_type
- Fixed AUC calculation to use composite IDs that combine subject ID and treatment group, ensuring correct AUC values for each unique subject-treatment combination 
- Fixed `plot_auc` function error where `is.named()` function was called but not defined
- Fixed event counting bug in `survival_statistics` function where event totals were incorrectly calculated when generating the formatted summary table
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