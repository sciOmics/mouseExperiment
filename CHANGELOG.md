# Changelog

All notable changes to the mouseExperiment package will be documented in this file.

## [Unreleased]

### Added
- New function `plot_combination_index()` to visualize Combination Index over time with synergy indicators

### Fixed
- Fixed `forest_plot()` to properly handle both `HR` and `Hazard_Ratio` column naming from `survival_statistics()` results
- Improved compatibility between column naming conventions in various plotting functions
- Fixed Events/Total count calculation in `survival_statistics()` to correctly count unique subjects and their events per treatment group

### Changed
- Enhanced summary output in `tumor_growth_statistics()` to provide detailed description of statistical tests and methods used
- Improved growth rate calculation description in result summaries to clarify the log-transformation and interpretation

## [0.1.0] - 2025-03-27

### Added
- Initial release of the mouseExperiment package
- Functions for tumor growth statistics
- Functions for survival analysis
- Functions for combination index calculation
- Plotting functions for AUC, growth rates, and survival curves
- Example synthetic datasets