# Changelog

All notable changes to the mouseExperiment package will be documented in this file.

## [0.1.0] - 2024-03-27

### Added
- Initial release of the mouseExperiment package
- Added `plot_auc` function to create visualizations for Area Under the Curve (AUC) data
- Added `tumor_growth_statistics` function for analyzing tumor growth data
- Properly documented all datasets with roxygen2 comments

### Fixed
- Fixed issue in `plot_auc` function by replacing `annotate_figure` with a simpler approach using `ggpubr::ggarrange`
- Improved error handling in `tumor_growth_statistics` function
- Fixed model type argument matching to properly handle "lme4" and "auc" types
- Fixed dataset exports by adding proper documentation and using LazyData: true
- Added examples for all datasets 