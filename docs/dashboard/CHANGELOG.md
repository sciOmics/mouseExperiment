# Changelog

All notable changes to the mouseExperiment package will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)

## [v26.04.005] 2026 April 20th

### Added
- **Toxicity tab** — new top-level analysis tab for body weight toxicity
  - Body Weight Trajectory: spaghetti or mean±SEM/SD with shaded ribbons
  - Mixed Model: longitudinal LME with tumor-weight adjustment, covariates (volume, sex, initial mass), auto-fallback random effects
  - Body Weight AUC: trapezoidal AUC of % weight change, boxplot + jitter, nadir analysis, pairwise comparisons
  - Weight Loss Survival: KM curves for time-to-threshold weight loss, log-rank, optional Cox PH, risk table
- **Efficacy-Toxicity tab** — new top-level analysis tab combining tumor and weight data
  - Therapeutic Window Metric (TWM): TGI/MaxWeightLoss% bar chart with noise-floor handling
  - Safety-Efficacy Plot: per-mouse bivariate scatter with group ellipses/centroids, 3 efficacy metrics (TGI, tumor AUC, log cell kill)
  - Total Benefit Score: B = AUC(efficacy) − λ×AUC(toxicity) with interactive λ slider
  - Weight-Corrected TGI: side-by-side comparison after excluding mice exceeding safety threshold
  - Correlation Heatmap: efficacy vs toxicity metrics
  - Combined Report: merged summary table with CSV export
- Weight and Sex column auto-detection in data upload
- Tab enable/disable: Toxicity requires weight column, Efficacy-Toxicity requires both weight and volume

## [v26.04.004] 2026 April 18th

### Added
- "Tumor Growth Inhibition vs Dose" plot in Dose Response module
  - Per-observation %TGI computed as `(1 - Volume / mean_control_volume) * 100`
  - Dashed reference line at 0% inhibition
- Trend line fit selector (shared by both Dose Response plots): Linear or Non-linear (LOESS)
- Point size slider for Dose Response scatter plots
- Confidence interval shading toggle for Dose Response trend lines

### Changed
- Renamed Dose Response "Scatter Plot" to "Tumor Volume vs Dose"
- Removed "Boxplot by Dose" and "Mean by Dose (Bar)" plot options from Dose Response module
- Appearance controls (point colour, line colour, point size, CI toggle) now shared across both Dose Response plots

### Fixed
- Dose Response module no longer requires running Tumor Growth Analysis first — `ensure_preprocessed()` now resolves dose, ID, and cage column mappings from settings

## [v26.04.003] 2026 April 17th

### Added
- Enabled Dose Response and Drug Synergy analysis tabs (previously disabled)
- Tumor Growth Inhibition (TGI) plot in Drug Synergy module
  - Days on x-axis, percent inhibition on y-axis
  - Per-mouse TGI computed from raw volume data against control group mean
  - Combination group displayed as dashed line; monotherapies as solid
  - Configurable error bars (SEM, SD, 95% CI, or None)
- Toggleable region shading and labels on Combination Index plot (Synergistic / Additive / Antagonistic bands)
- CI plot now rebuilt dynamically to support interactive option controls
- Mobile incompatibility overlay — on screens ≤ 768 px a full-screen message recommends using a desktop, with a "Continue anyway" dismiss button
- Viewport meta tag for proper mobile scaling

### Changed
- Updated Drug Synergy analysis to use proper Loewe Additivity model (Berenbaum, 1989)
  - Expected fractional effect: min(FE_A + FE_B, 1.0)
  - Combination Index: min(FE_A + FE_B, 1) / FE_combo
- Updated synergy info box, labels, and analytics tracking to reflect Loewe terminology

### Fixed
- Corrected synergy metrics info box text to match actual package output
- Fixed Combination Index formula: numerator now uses capped Loewe expected FE — previously the uncapped sum could inflate CI and falsely indicate antagonism when individual drug effects sum to > 1.0

## [v26.04.002] 2026 April 8th 

### Added
- Comprehensive Umami analytics event tracking across all modules via a new `umami-events.js` helper library
  - `trackAnalysisStart` / `trackAnalysisComplete` for tumor growth, survival, synergy, power, and dose-response modules
  - `trackDownload` for all CSV, PNG, and ZIP download handlers
  - `trackDataUpload` and `trackDemoDataLoaded` on file/demo load
  - `trackColumnMapping` when all required column mappings are set
  - `trackReferenceGroupSelected` and `trackError` for validation paths
  - `trackTabNavigation` for main navbar; `trackFeatureAccess` when a coming-soon tab is clicked
  - `model_type_selected` and `plot_type_selected` events for tumor growth model/plot switching
  - `session_saved`, `session_restored`, and `analysis_log` download events
  - `dark_mode_toggled` event fired directly from `darkmode.js`

### Fixed
- Loading new data now resets all downstream analysis results — tumor growth, survival, synergy, power, and dose-response results are cleared whenever a new file or demo dataset is loaded, preventing stale results from a previous dataset from remaining visible
