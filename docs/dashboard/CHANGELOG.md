# Changelog

All notable changes to the mouseExperiment package will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)

## [v26.05.017] 2026 May 20th

### Added
- **Synergy — Bayesian Combination Synergy (brms)** — new `"Bayesian (brms)"` option in the Drug Synergy analysis mode (uses `bayesian_synergy()` v0.4.2). Sidebar controls: Transform (log / sqrt / none), Prior Strength preset, Chains, Iterations. Result tabs: TGI Summary (per-group posterior median ± 95 % CrI), Bliss Summary (excess posterior distribution + P(synergy)), Loewe Summary (CI posterior distribution + P(CI < 1) + interpretation string), Synergy Table (all six rows including Bliss-expected and Loewe-expected), MCMC Diagnostics, Synergy Plot (bar chart with Bliss/Loewe reference lines), Posterior Distribution (density overlays for Bliss excess and Loewe CI).
- **Toxicity — Bayesian Therapeutic Window (brms)** — new `"Bayesian TWM (brms)"` option in the Therapeutic Window Metric mode (uses `bayesian_therapeutic_window()` v0.4.1). Requires both a Bayesian tumor growth result (`return_model = TRUE`) and a Bayesian body weight result (`return_model = TRUE`) to already be available in the session. Result tabs: TWM Table (per-group TWM median + 95 % CrI), TGI Summary, Weight-Loss Summary, TWM Plot (forest plot with break-even reference line at TWM = 1), Efficacy vs Safety Scatter (TGI vs weight-loss % with TWM = 1 isoline).
- **Tumor Growth — Bayesian Dose-Response (brms)** — new `"Bayesian Hill/Emax (brms)"` option in the Dose-Response analysis mode (uses `bayesian_dose_response()` v0.4.0). Sidebar controls: Prior Strength, Chains, Iterations. Result tabs: DR Parameters (EC50/Emax/Hill posterior median + 95 % CrI), Dose-Response Summary, MCMC Diagnostics, Dose-Response Curve (posterior median ± CrI ribbon + EC50 line), Prior vs Posterior.

## [v26.05.016] 2026 May 19th

### Changed
- **Tumor Growth — Bayesian LMM: cage random-intercept and prior controls** — updated to reflect `bayesian_tumor_growth()` v0.3.7 changes:
  - Prior Strength selector now includes `"Skeptical"` as a new option and sets it as the default (replaces `"Weakly informative"`). A `"Manual"` option is reserved for future sidebar expansion.
  - New **Cage Effects** checkbox (default: on) in the Bayesian sidebar, wired to the `include_cage_effect` argument. Displayed only when a cage column is mapped in Data Upload.
  - **Prior vs Posterior** plot tab dynamically populated from `res$prior_posterior_plot` (available when `plots = TRUE`); shows density overlay of prior (grey) and posterior (blue) per treatment coefficient.
  - **Residuals vs Day** plot tab dynamically populated from `res$residuals_plot`; faceted by treatment group with loess smoother.

### Added
- **Survival — Bayesian AFT (brms)** — new `"Bayesian AFT (brms)"` option in the Survival Statistical Model selector (uses `bayesian_survival()` v0.3.8). Sidebar controls: Family (Weibull / Log-normal / Exponential / Gamma), Prior Strength, Chains, Iterations, Include Cage Frailty checkbox. Result tabs: Treatment Effects (Time Ratio + HR), Posterior Summary, MCMC Diagnostics, Survival Curves (posterior credible bands + KM), Prior vs Posterior.
- **Toxicity — Bayesian Body Weight (brms)** — new `"Bayesian Body Weight (brms)"` option under Toxicity → Body Weight mode (uses `bayesian_body_weight()` v0.3.9). Result tabs: Treatment Effects, Weight Loss Summary (posterior % change first-to-last day with 95 % CrI), Weight Trajectory (posterior median ± CrI ribbon + observed points), MCMC Diagnostics, Prior vs Posterior.

## [v26.05.015] 2026 May 14th

### Added
- **Tumor Growth — Bayesian LMM (brms) analysis mode** — new "Bayesian LMM (brms)" option in the Statistical Model selector. When selected, the analysis runs asynchronously in a background worker via `future`/`promises` so the app remains fully responsive during the 3–12 minute Stan compilation and MCMC sampling. A spinner notification is shown while the model runs.
  - Sidebar controls: Prior Strength preset (Weakly informative / Informative / Diffuse), MCMC chains (default 4), Iterations per chain (default 2 000), and an amber timing notice.
  - Random Effects selector is now visible for both `lme4` and `bayes` modes (previously lme4-only).
  - Two Bayesian-specific result tabs are dynamically inserted when Bayesian mode is active and removed when switching away: **Posterior Check** (posterior predictive density overlay via `brms::pp_check`) and **MCMC Diagnostics** (Rhat, Bulk_ESS, Tail_ESS, Converged flag table + trace plots).
  - The **ANOVA** tab displays an informative "not applicable" message in Bayesian mode (ANOVA has no meaning for posterior inference).
  - The **Diagnostics** tab (Q–Q and Residuals vs Fitted plots) redirects users to the Posterior Check tab in Bayesian mode.
  - Treatment Effects and Pairwise Comparisons tabs populate from the same column structure as the lme4 path (95 % CrI for Bayes, 95 % CI for lme4).
  - **Bayesian LMM Mode** section added to the Tumor Growth Info tab describing priors, outputs, and interpretation guidance.
- **Async execution infrastructure** — `future::plan(multisession, workers = 1L)` initialised at app startup in `app.R`; cleaned up on exit via `on.exit`. Prevents blocking the Shiny event loop during long-running Stan fits.
- `future (>= 1.33.0)`, `promises (>= 1.2.0)`, `brms (>= 2.19)`, `bayesplot (>= 1.10)` added to `Suggests` in DESCRIPTION.

## [v26.05.014] 2026 May 14th

### Changed
- **Power Analysis — LMM Simulation defaults and labels** — growth-rate and variance-component inputs now use log-scale defaults (`Control log-growth rate = 0.15`, `Treatment reduction = 0.10`, `Random intercept SD = 0.20`, `Random slope SD = 0.05`, `Residual SD = 0.10`) matching the corrected `apriori_power_simulation()` which generates log-normal data and fits `log(Volume)`.
- **Toxicity — Therapeutic Window Metric formula** — TWM denominator now reads "Mean Body Weight Loss" throughout the UI (info panel, tab description, plot subtitle) and the summary text correctly references `Mean_Pct_Weight_Loss` and `TGI` columns from the updated package.

### Fixed
- **Toxicity — TWM summary text** — the combined-report text block referenced non-existent columns `TGI_Pct` and `Max_Weight_Loss_Pct`; corrected to `TGI` and `Mean_Pct_Weight_Loss`.

## [v26.05.013] 2026 May 14th

### Removed
- **Power Analysis — Post-hoc (observed data) mode** — eliminated in line with the removal of `post_power_analysis()` from the `mouseExperiment` package. Post-hoc power is statistically uninformative (a monotone function of the p-value) and should not be offered as an analysis option. Removed: the "Post-hoc" radio option, post-hoc sidebar controls (method selector, simulation count, custom effect sizes), post-hoc results panel (Summary, Effect Sizes, Power Analysis, Sample Size Recommendations, Power Curves tabs), and all associated server outputs and download handlers (`pw_summary`, `pw_effect_sizes_table`, `pw_power_table`, `pw_sample_sizes_table`, `pw_plot_curves`, `download_pw_effect_sizes`, `download_pw_power_table`, `download_pw_sample_sizes`, `download_pw_plot`).

### Changed
- **Power Analysis — default mode** changed from "Post-hoc (observed data)" to "A priori — Analytic".
- **Power Analysis — Target Power** panel now shown only for "A priori — Analytic" mode (previously also shown for post-hoc mode).
- **Power Analysis — results area** is now a single unconditional tabset showing Scenario Table, Power Curves, and SD Sensitivity (previously two conditional panels switching between post-hoc and a priori tabs).

## [v26.05.012] 2026 May 5th

### Added
- **Power Analysis — A priori Analytic mode** — new sidebar mode alongside the existing post-hoc. Users enter Cohen's d values (one or several, comma-separated) or a raw mean difference + pooled SD. Supports two-sample t-test (two groups) or one-way ANOVA omnibus (k ≥ 3 groups via `pwr::pwr.anova.test`). Sub-modes: "Find N" (required N per group for each target power × alpha combination) and "Find Power" (achieved power at user-specified N values). A "Seed SD from data" button fills the pooled SD from the currently loaded dataset. For mean-difference inputs, an SD sensitivity table shows how required N changes as the assumed SD is varied ±20% and ±40%, communicating uncertainty in the pilot SD estimate.
- **Power Analysis — A priori LMM Simulation mode** — simulates longitudinal tumour-growth data from user-supplied growth parameters (baseline volume, control rate, treatment rate-reduction, random effects, residual SD, timepoint schedule). Fits `Volume ~ Treatment * Day + (Day|ID)` via `lme4::lmer` and uses a likelihood-ratio test for the Treatment × Day interaction; power = proportion of simulations where p < α. Evaluates power across a user-specified range of N per group. Progress bar updates during each (N, alpha) scenario.
- **Power Analysis — Analysis Mode selector** — `radioButtons` at the top of the sidebar switches between "Post-hoc (observed data)", "A priori — Analytic", and "A priori — LMM Simulation". Result panels are conditionally displayed; post-hoc tabs only show for post-hoc mode, a priori tabs (Scenario Table, Power Curves, SD Sensitivity) only show for a priori modes.

## [v26.05.011] 2026 May 5th

### Added
- **Power Analysis tab** — activated as a full tab; previously rendered as a greyed-out "Coming soon" placeholder even though the module was fully implemented. The tab is now enabled whenever data is loaded (same condition as Tumor Growth). The analysis itself requires preprocessed data from a prior Tumor Growth run; a warning is shown if that has not been done yet.

## [v26.05.010] 2026 May 5th

### Removed
- **Toxicity — Log Cell Kill efficacy metric** — removed from the ET mode sidebar radio buttons, the Info tab Efficacy Metrics list, the scatter-plot label switch, and the Shared Parameters description; matches the removal in mouseExperiment 0.3.2

## [v26.05.009] 2026 May 5th

### Added
- **Toxicity — BW Trajectory combined view** — new display mode "Individual + Group Mean" overlays individual mouse trajectories (faint) with bold group-mean lines; default display changed to "Individual mice"
- **Toxicity — BW AUC plot controls** — replaced "Show box plot" checkbox with TG-style centre-bar selector (Mean / Median / None), error-bar selector (SEM / SD / 95% CI / None), and jitter-width slider; y-axis label updated to include units (%·days)
- **Toxicity — ET plot group controls** — group order, colour, and palette-preset controls now available for all ET mode plots (Therapeutic Window, Safety-Efficacy Scatter, Total Benefit Score, Weight-Corrected TGI)
- **Toxicity — Safety-Efficacy Scatter controls** — added centroid size slider, ellipse line-width slider, ellipse line-style selector (solid / dashed / dotted), and "Shape by group" checkbox

### Fixed
- **Toxicity — Safety-Efficacy Scatter only showing reference group** — `efficacy_toxicity_bivariate()` now uses composite `Treatment⁠ID` keys when iterating per-mouse efficacy, preventing shared IDs across groups from collapsing to a single group

### Removed
- **Toxicity — Correlation Heatmap** — removed from ET plot selector, UI, server, and download handler

## [v26.05.008] 2026 May 5th

### Added
- **Toxicity Analysis — Results tabs restructured** — the single scrollable "Results" tab has been replaced with separate, focused tabs for each analysis section:
  - *Body Weight mode*: **Mixed Model** (fixed effects + estimated marginal means), **AUC** (AUC summary, nadir, pairwise comparisons), **Weight Loss KM** (time-to-threshold event data)
  - *Efficacy-Toxicity mode*: **Therapeutic Window**, **Safety-Efficacy**, **Benefit Score**, **Corrected TGI**, **Combined Report** — dynamically inserted/removed when the Statistical Model selector is changed
- **Toxicity Analysis — Post-hoc comparison type** — new sidebar dropdown in Body Weight mode: "vs Reference Group" (BH-corrected Welch t-tests restricted to comparisons involving the control) or "All Pairwise (BH-corrected)" (all group pairs); the AUC Comparisons table, Comparisons title, and AUC Comparison Forest plot all react to this selector
- **Toxicity Analysis — BH p-value correction** — Benjamini-Hochberg correction is now actually applied to the AUC comparison p-values (previously the Info tab stated BH correction was used but the displayed p-values were unadjusted); a `P_Adjusted_BH` column is now shown alongside the raw `P_Value`
- **Toxicity Analysis — AUC Comparison Forest Plot** — new plot type in the BW Plots selector: effect-size forest plot of mean AUC differences with 95% CI, coloured by BH-corrected significance; respects the comparison type selector
- **Toxicity Analysis — Diagnostics tab** — new tab showing Q-Q plot of mixed-model residuals and residuals-vs-fitted scatter, available in Body Weight mode; displays a notice in ET mode
- **Toxicity Analysis — Download All tab** — new tab with a single button that exports all available result tables as a ZIP file of CSVs
- **Toxicity Analysis — Additional download buttons** — new per-table download buttons: EMMs CSV (Mixed Model tab), Nadir CSV and Comparisons CSV (AUC tab)
- **Toxicity Analysis — Reactive tab info** — AUC tab description and comparisons title now update dynamically to reflect the current comparison type

## [v26.05.007] 2026 May 5th

### Added
- **Survival Analysis — Endpoint source always visible** — the "Use Event/Censor column" / "Derive from tumour volume threshold" radio group is now shown unconditionally in the sidebar; previously it was hidden behind a static text message when no direct Volume column was mapped

### Fixed
- **Survival — volume-threshold endpoint with calculated volume** — when the Data Upload tab is configured to calculate tumour volume from length × width measurements (rather than mapping a direct Volume column), the volume-threshold event derivation now correctly uses the calculated `Volume` column produced by `ensure_preprocessed()`; previously `rv$settings$col_volume` was empty in this case, causing a validation error and blocking the analysis

## [v26.05.006] 2026 May 5th

### Added
- **Survival Analysis — Endpoint source always visible** — the "Use Event/Censor column" / "Derive from tumour volume threshold" radio group is now shown unconditionally in the Analysis Parameters sidebar; previously it was replaced by a static text note when no Volume column was configured in Data Upload. When "Derive from tumour volume threshold" is selected and no volume column is available, an inline amber warning directs the user to Data Upload instead of silently hiding the control.

## [v26.05.005] 2026 May 5th

### Added
- **Body Weight AUC plot** — new "Show box plot" checkbox (default: on) to toggle the boxplot layer, leaving only individual data points when unchecked; new "Point size" slider (0.5–6, default: 2); dot colours now follow the Groups palette when the boxplot is hidden
- **Toxicity Analysis — Info tab** — new first tab in the results panel with structured explanations of all analyses, metrics, and parameters, organised into four sections:
  - *Overview*: when to use Body Weight Analysis vs Efficacy-Toxicity Analysis
  - *Body Weight Analysis*: mixed-effects model formula and covariate effects, REML vs ML, AUC derivation, nadir analysis, KM weight-loss curve, and threshold definition
  - *Efficacy-Toxicity Analysis*: Final TGI, Tumor AUC, and Log Cell Kill formulas; Therapeutic Window Metric; Safety-Efficacy Scatter layers; Total Benefit Score formula and λ guidance; Weight-Corrected TGI rationale; Correlation Heatmap interpretation
  - *Shared Parameters*: reference group auto-detection, tumor weight adjustment, tumor density, and baseline day

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
