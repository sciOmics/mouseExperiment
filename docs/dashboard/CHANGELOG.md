# Changelog

All notable changes to the mouseExperiment package will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)

## [v26.05.005] 2026 May 5th

### Added
- **Body Weight AUC plot** — new "Show box plot" checkbox (default: on) to toggle the boxplot layer, leaving only individual data points when unchecked; new "Point size" slider (0.5–6, default: 2); dot colours now follow the Groups palette when the boxplot is hidden

### Fixed
- **Safety-Efficacy Scatter plot** — corrected column name mismatches between `efficacy_toxicity_bivariate()` return value and the `renderPlot` code (`Max_Pct_Weight_Loss`, `Efficacy_Mean`, `Toxicity_Mean`); previously returned `[object Object]`
- **Total Benefit Score plot** — corrected column reference from `Benefit` to `Benefit_Score`; previously returned `[object Object]`
- **Weight-Corrected vs Standard TGI plot** — corrected column references from `Uncorrected_TGI`/`Corrected_TGI` to `TGI_Uncorrected`/`TGI_Corrected` (merge-suffix order from `weight_corrected_tgi()`); previously returned "arguments imply differing number of rows"
- **Correlation Heatmap plot** — corrected AUC column references from `ea$AUC`/`ta$AUC` to `ea$Efficacy_AUC`/`ta$Toxicity_AUC`; previously returned "undefined columns selected"

## [v26.05.004] 2026 May 5th

### Fixed
- **Body Weight Trajectory "Individual mice" plot** — spaghetti lines now correctly group by `paste(Treatment, ID)` composite key; previously, animals sharing the same ID across treatment groups were joined into a single crossing line
- **Body Weight Trajectory, AUC, and Weight Loss KM plots** — added group recoloring and reordering controls (drag-and-drop order list + per-group colour pickers with palette presets) matching the Tumor Growth and Survival tab functionality

## [v26.05.003] 2026 May 4th

### Added
- **Survival tab — IACUC endpoint context and volume-threshold event derivation**
  - Page header now shows subtitle "(Time-to-Endpoint)" and a warning banner explaining that in most preclinical mouse studies the recorded endpoint is time-until-tumour-volume-threshold (not true death), with a note on competing risks and the recommended pragmatic approach
  - New **Endpoint** section in the Analysis Parameters sidebar:
    - When a Tumour Volume column is configured in Data Upload, a radio group lets users choose between "Use Event/Censor column" (existing behaviour) and "Derive from tumour volume threshold" (new)
    - Volume mode: `numericInput` for the IACUC volume limit (default 2000 mm³); automatically computes per-animal Time = first day volume ≥ threshold (Event = 1) or last observation day (Event = 0) if never exceeded
  - New **Natural Death** checkbox in the sidebar (default: checked); text explains the competing-risks tradeoff and the pragmatic recommendation
  - Summary output now shows "Event Source: Derived from tumour volume threshold (>= X mm³)" when in volume mode, replacing the "Event Column:" line
  - KM plot title and axis labels update dynamically when volume mode is active ("Time-to-Tumour-Volume-Threshold by Treatment", "Time to Endpoint (Days)", "Proportion Below Volume Threshold")
- **Data Upload — updated Event/Censor Status column descriptions**
  - `help_tip` now reads: "1 = endpoint reached (tumour exceeded IACUC volume limit, or animal euthanised/died), 0 = censored"
  - Below-input text notes that volume-threshold derivation is available in the Survival sidebar
  - Help page bullet updated to mention tumour volume threshold as the typical endpoint

### Changed
- Summary tab header renamed from "Survival summary" to "Time-to-Endpoint summary"
- Cox Model "Events" bullet updated to mention tumour volume threshold; Median_Survival bullet updated to say "reached the endpoint" rather than "experienced the event"
- Risk Table N_Censored bullet updated to say "without reaching the endpoint (e.g. still within tumour volume limit at study end)"

## [v26.05.002] 2026 May 4th

### Fixed
- **Body Weight AUC and Weight Loss KM plots showing only the control group** — `body_weight_auc()` and `weight_loss_threshold()` in the `mouseExperiment` package now use a composite `Treatment + ID` key internally when grouping observations per mouse. Previously, if the same ID value (e.g. an ear tag) appeared in multiple treatment groups, all mice were incorrectly attributed to the first group encountered (typically the vehicle control). The fix ensures each treatment-group × ID combination is processed as a separate animal.

## [v26.05.001] 2026 May 4th

### Changed
- **Toxicity tab** — merged the former "Efficacy-Toxicity" tab into the Toxicity tab
  - New **Statistical Model** dropdown in Analysis Parameters sidebar (matches Tumor Growth pattern):
    - "Body Weight Analysis" — mixed-effects model, AUC, and time-to-threshold KM curve
    - "Efficacy-Toxicity Analysis" — therapeutic window metric, safety-efficacy scatter, total benefit score, weight-corrected TGI, correlation heatmap, combined report
  - Replaced four separate sub-tabs (Trajectory, Mixed Model, AUC, KM) with a unified **Plots** tab containing a plot-type selector, plus a **Results** tab with all data tables
  - Plot controls (type selector, display style, scaling, export) consolidated into a single left-hand sidebar panel within the Plots tab, consistent with Tumor Growth module layout
  - Sidebar parameters are now model-specific: Body Weight options (covariates, REML/ML, threshold, baseline day) shown only when "Body Weight Analysis" is selected; Efficacy-Toxicity options (efficacy metric, lambda, safety threshold) shown only when "Efficacy-Toxicity Analysis" is selected
- Removed standalone **Efficacy-Toxicity** top-level tab from navigation bar

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
