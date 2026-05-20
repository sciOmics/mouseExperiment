# Changelog

All notable changes to the mouseExperiment package will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.4.1] - 2026-05-20

### Added
- `bayesian_therapeutic_window()` — new exported function that propagates full posterior uncertainty through the Therapeutic Window Metric (TWM = TGI / |weight-loss %|) via draw-wise arithmetic on posterior predictive samples from `bayesian_tumor_growth()` and `bayesian_body_weight()` results. Key features:
  - Accepts `tg_result` (from `bayesian_tumor_growth(return_model = TRUE)`) and `bw_result` (from `bayesian_body_weight(return_model = TRUE)`); validates that both contain a `brmsfit` object.
  - Calls `brms::posterior_epred(..., re_formula = NA)` on each model for population-level predictions, back-transforms from the model scale using each model's `transform_used` field (`log` → `exp`, `sqrt` → `x^2`, otherwise identity).
  - TGI per draw: `1 - V_trt^(d) / V_ctrl^(d)` using the posterior predictive volume at the endpoint day. Weight-loss % per draw: `(W_last^(d) - W_first^(d)) / W_first^(d) × 100`.
  - `noise_floor` (default 1 %) clamps the TWM denominator so that near-zero weight-loss values do not inflate the metric to infinity: `TWM^(d) = TGI^(d) / max(|WL%^(d)|, noise_floor) / 100`.
  - Draws from the two independently fitted models are aligned draw-by-draw (using the minimum draw count when chains differ), yielding a proper posterior distribution of TWM with a direct probability interpretation.
  - `twm_table`: one row per treated group with `TWM_Median`, `TWM_Lower` (2.5 % CrI), `TWM_Upper` (97.5 % CrI), `TGI_Median`, `WL_Pct_Median`, `N_Draws`.
  - `tgi_summary` and `wl_summary`: per-group posterior medians and 95 % CrI for TGI and weight-loss % respectively, including the reference group.
  - `twm_plot`: forest plot of TWM ± 95 % CrI per treated group; vertical dashed line at TWM = 1 (break-even efficacy-safety boundary).
  - `tgi_wl_plot`: scatter of posterior-median TGI vs weight-loss % per group with the TWM = 1 isoline (dotted red) and the noise-floor vertical dashed reference.
  - Returns `model_type_used = "bayes_twm"`.
- `make_tg_for_twm()` test fixture added to `tests/testthat/test-bayesian_therapeutic_window.R`: 2 groups × 5 mice × 4 time-points (days 0, 7, 14, 21); Control log_slope = 0.12 (fast growth), TreatmentA log_slope = 0.03 (strongly inhibited), yielding positive TGI > 0 for TreatmentA. Combined with `make_bw_simple()` for weight-loss data.
- Tests for `bayesian_therapeutic_window()` in `tests/testthat/test-bayesian_therapeutic_window.R` (18 tests): return structure, `model_type_used`, `twm_table` columns and row count, exclusion of reference group, `TWM_Median` finiteness, positive `TGI_Median`, CrI brackets median, `N_Draws` positive integer, `tgi_summary` and `wl_summary` row count and columns, control TGI near zero, TreatmentA negative weight loss, summary metadata, `twm_plot`/`tgi_wl_plot` ggplot class, input-validation errors (NULL model, invalid reference group), and `plots = FALSE` returns NULL plots.

## [0.4.0] - 2026-05-19

### Added
- `bayesian_dose_response()` — new exported function for Bayesian Hill/Emax dose-response modelling via `brms`. Fits the nonlinear inhibition model TGI = Emax / (1 + (EC50/Dose)^Hill) to per-mouse Tumour Growth Inhibition values computed relative to the vehicle control mean at the endpoint day. Key features:
  - Log/logit reparameterisation ensures positivity constraints without hard bounds: Emax = inv_logit(logEmax) ∈ (0,1), EC50 = exp(logEC50) > 0, Hill = exp(logHill) > 0. The model is fitted on treated animals only (Dose > 0); the reference group defines TGI = 0 at Dose = 0 by construction.
  - Data-adaptive EC50 prior centred at the log geometric mean of the observed non-zero doses, so the prior is automatically reasonable across different dose scales.
  - Prior-strength presets: `"skeptical"` (default; narrow priors centred at Emax ≈ 0.82, Hill = 1, EC50 at geometric mean dose), `"weakly_informative"`, `"diffuse"`, and `"manual"` (supply `prior_emax`, `prior_ec50`, `prior_hill`, `prior_sigma` directly as brms prior strings on the log/logit scale).
  - `dr_parameters` table: EC50, Emax, and Hill back-transformed to interpretable scale with posterior median and 95 % CrI.
  - `dose_response_summary`: per-dose observed TGI mean/SD alongside posterior-predicted median and 95 % CrI.
  - `dose_response_curve_plot`: posterior median ± 95 % CrI ribbon over a dense dose grid with observed points; vertical red dashed line at the posterior median EC50 and horizontal reference at TGI = 0.5.
  - `prior_posterior_plot`, `pp_check_plot`, `mcmc_trace_plot` all returned when `plots = TRUE`.
  - Returns `model_type_used = "bayes_dr"`, `tgi_data` (analysis data frame), and `control_mean_volume`.
- Tests for `bayesian_dose_response()` in `tests/testthat/test-bayesian_dose_response.R` (17 tests): return structure, `model_type_used`, `brmsfit` class, EC50/Emax/Hill presence and biological validity (Emax ∈ (0,1), EC50 > 0, Hill > 0, CI brackets median), dose-response summary columns and row count, monotone observed TGI with dose, `tgi_data` content, `control_mean_volume` finiteness, MCMC diagnostics, summary metadata, input-validation errors, `return_model = FALSE`, and `dose_response_curve_plot` ggplot class.

## [0.3.9] - 2026-05-19

### Added
- `bayesian_body_weight()` — new exported function for Bayesian linear mixed-effects analysis of longitudinal body weight data using `brms`. Mirrors the `analyze_body_weight()` API with Bayesian-specific additions:
  - Default `transform = "none"` (body weight is approximately normal; `"log"` and `"sqrt"` also supported).
  - Optional tumor weight subtraction: `adjust_tumor_weight = TRUE` subtracts `Volume / 1000 × tumor_density` from gross body weight before modelling.
  - `weight_loss_summary` — per-group posterior percentage weight change from the first to last study day (median + 95 % CrI), computed via `brms::posterior_epred(..., re_formula = NA)` and back-transformed to the original weight scale.
  - `weight_trajectory_plot` — population-level posterior median ± 95 % CrI ribbon over the full study timeline with observed data points overlaid per treatment group.
  - Same prior-strength presets as `bayesian_tumor_growth()` (default `"weakly_informative"`), cage random-intercept support (`include_cage_effect`), `prior_posterior_plot`, `residuals_plot`, `pp_check_plot`, `posterior_dist_plot`, `credible_intervals_plot`, and `mcmc_trace_plot`.
  - Returns `model_type_used = "bayes_bw"` and `weight_data` (prepared analysis data frame, analogous to `analyze_body_weight()$weight_data`).
- `make_bw_simple()` test fixture added to `tests/testthat/helper-fixtures.R`: 2 groups × 5 mice × 5 time-points (days 0, 7, 14, 21, 28); Control stable at 22 g, TreatmentA losing −0.08 g/day (≈−10 % by day 28).
- Tests for `bayesian_body_weight()` in `tests/testthat/test-bayesian_body_weight.R` (15 tests): return structure, `model_type_used`, transform, `brmsfit` class, treatment-effects columns, reference-group labelling, MCMC diagnostics, `weight_loss_summary` columns and sign-check for the treated group, `weight_data` row count, summary metadata, input-validation errors, `return_model = FALSE`, and `weight_trajectory_plot` ggplot class.

## [0.3.8] - 2026-05-19

### Added
- `bayesian_survival()` — new exported function for Bayesian parametric survival / accelerated-failure-time (AFT) analysis using `brms`. Key features:
  - Families: `"weibull"` (default), `"lognormal"`, `"exponential"`, `"gamma"` — all fitted in AFT parameterisation.
  - Effect metric: Time Ratio (TR = exp(coef)); for Weibull and exponential, a Hazard Ratio (HR = TR^(−shape)) is also computed draw-wise from the joint posterior and summarised as a posterior median.
  - Frailty model: optional cage-level random intercept (`include_frailty = TRUE`) when `cage_column` is supplied, modelling inter-cage heterogeneity in baseline survival.
  - `treatment_effects` table schema is compatible with `survival_statistics()`: columns `Group`, `Time_Ratio`, `Lower_CL`, `Upper_CL`, `HR`, `Median_Survival`, `Events`, `Total`, `Event_Rate`, `Note`.
  - Censoring: automatically handles the brms censoring convention inversion (`event = 1` → `brms_cens = 0`; censored = 0 → `brms_cens = 1`).
  - Prior-strength presets (`"skeptical"` default through `"diffuse"`) with auxiliary-class priors (`"shape"` for Weibull/gamma, `"sigma"` for lognormal, none for exponential); manual prior specification via `prior_b`, `prior_intercept`, `prior_sd`, `prior_aux`.
  - `sample_prior = "yes"` enables `prior_posterior_plot` (grey prior vs blue posterior density overlay for each treatment coefficient).
  - Plots: `pp_check_plot`, `posterior_dist_plot`, `prior_posterior_plot`, `mcmc_trace_plot`, `survival_curve_plot` (200-draw posterior credible bands with KM overlay).
  - Returns `model_type_used = "bayes_survival"`, `family_used`, `frailty_used`, `survival_data` (standardised analysis data frame).
- Tests for `bayesian_survival()` in `tests/testthat/test-bayesian_survival.R` (13 tests plus 2 standalone): return structure, `model_type_used`, family, treatment-effects columns, reference TR = 1 and HR = 1, non-reference TR < 1 and Upper_CL < 1, event counts, positive-finite median survival, posterior_summary Rhat, MCMC diagnostics Converged logical, survival_data row count, missing-column error, invalid reference-group error, summary metadata, lognormal HR = NA, frailty_used = FALSE without cage column.

## [0.3.7] - 2026-05-19

### Added
- `bayesian_tumor_growth()` — cage random-intercept modelling via new `include_cage_effect` parameter (default `TRUE`). When a `cage_column` is supplied and `include_cage_effect = TRUE`, the random-effects structure is `(1|cage/ID)` for intercept-only or `(Day|ID) + (1|cage)` for random slopes, capturing between-cage heterogeneity. Set `include_cage_effect = FALSE` to use `cage_column` solely for composite mouse-key construction without modelling it.
- `bayesian_tumor_growth()` — `"skeptical"` prior preset added and set as the new default (`b ~ N(0, 0.25)`, `sd/σ ~ Exponential(2)`). Expresses the scientifically motivated belief that treatment effects are small and requires stronger evidence to support large posterior differences. Previous default `"weakly_informative"` remains available.
- `bayesian_tumor_growth()` — `"manual"` prior-strength option: supply `prior_b`, `prior_intercept`, `prior_sd`, and `prior_sigma` as brms prior strings for full control over the prior specification.
- `bayesian_tumor_growth()` — `prior_posterior_plot` in the return list: overlaid density plot of prior (grey) and posterior (blue) for each treatment-effect coefficient, with a dashed zero-effect reference line. Produced by fitting with `sample_prior = "yes"` and calling the new internal helper `bayes_prior_posterior_plot()`.
- `bayesian_tumor_growth()` — `residuals_plot` in the return list: posterior-mean residuals vs study day, faceted by treatment group with a loess smoother. Systematic curvature flags violations of the log-linear time assumption.
- Internal helper `bayes_prior_posterior_plot(model, treatment_column)` defined at the bottom of `R/bayesian_tumor_growth.R` (`@noRd`); shared by `bayesian_survival()` and `bayesian_body_weight()`.
- `summary$model_specification` now includes `random_effects` (the RE formula term) and `cage_effect_modelled` (logical).

## [0.3.6] - 2026-05-14

### Added
- `bayesian_tumor_growth()` — new exported function for Bayesian linear mixed-effects analysis of tumor growth data using `brms` (Stan backend). Accepts the same data format and column arguments as `tumor_growth_statistics()` and returns a compatible result structure (`treatment_effects`, `pairwise_comparisons`, `growth_rates`, `data_summary`) alongside Bayesian-specific outputs: `posterior_summary` (fixed-effect parameters with 95 % credible intervals, Rhat, ESS), `mcmc_diagnostics` (per-parameter Rhat, Bulk_ESS, Tail_ESS, Converged flag), `pp_check_plot` (posterior predictive density overlay), `posterior_dist_plot` (mcmc_areas for treatment parameters), `mcmc_trace_plot`, and `credible_intervals_plot` (forest plot of credible intervals). Three prior-strength presets: `"weakly_informative"` (default: b ~ N(0,1), σ ~ Exp(1)), `"informative"` (b ~ N(0,0.5)), `"diffuse"` (b ~ N(0,2.5)). Supports `log`, `sqrt`, and `none` volume transformations and both `intercept_only` and `slope` random-effects specifications. Requires `brms (>= 2.19)` and `bayesplot (>= 1.10)` (both in `Suggests`).
- `brms (>= 2.19)` and `bayesplot (>= 1.10)` added to `Suggests` in DESCRIPTION.
- Tests for `bayesian_tumor_growth()` in `tests/testthat/test-bayesian_tumor_growth.R`: return-structure coverage, column-name compatibility with the lme4 path, Rhat convergence check, `return_model = FALSE` path, and graceful error when `brms` is absent. All tests are skipped when `brms` is not installed.

## [0.3.5] - 2026-05-14

### Fixed
- `apriori_power_simulation()` — data are now generated on the **log scale** (`log(V) = log(baseline) + b0 + (rate + b1)*t + ε`; volumes exponentiated) and the LMM is fitted as `log(Volume) ~ Treatment * Day + (Day | ID)`, matching `tumor_growth_statistics()`. The previous linear data-generating model produced anti-conservative power estimates when the fitted model uses log-transformation. Function argument defaults updated to log-scale values: `control_growth_rate = 0.15`, `treatment_effect = 0.10`, `random_intercept_sd = 0.20`, `random_slope_sd = 0.05`, `residual_sd = 0.10`.
- `apriori_power_analysis()` — the Cohen's d → f conversion (`f = d / sqrt(2)`) used for ANOVA power (k ≥ 3 groups) is now fully documented: it assumes two extreme groups and overestimates f (and therefore power) when all treated groups are uniformly reduced vs. control. A `method_note` field is returned when `n_groups >= 3` to surface this at runtime.
- `tumor_auc_analysis()` — the `last_observation` (LOCF) AUC method was dimensionally incoherent: it returned `last_volume` (mm³) added to `extrapolated_value` (mm³·day). Fixed to compute the trapezoidal AUC for the observed period first, then add the LOCF extension (time gap × last volume) as a properly dimensioned mm³·day area.
- `survival_statistics()` — in the logrank fallback path (zero-event groups), the omnibus chi-squared p-value was assigned to every non-reference group. Replaced with per-pair log-rank tests (`survdiff()` restricted to each treatment vs. reference, df = 1) so each group receives its own p-value.
- `therapeutic_window_metric()` — TWM denominator changed from the single worst-case mouse in each group (`FUN = max`) to the group mean of per-mouse maximum weight loss (`FUN = mean`), making the metric less sensitive to outliers. Column renamed `Max_Pct_Weight_Loss` → `Mean_Pct_Weight_Loss`.
- `tumor_growth_statistics()` — the lme4 path now additionally computes estimated marginal means (EMMs) at five quintile study days (min, Q1, median, Q3, max) via `emmeans::emmeans()`, returned as `treatment_effects_over_time`. Marginalising at a single mean time point (the existing `treatment_effects`) discards the Treatment × Day interaction that is the primary analysis target.
- `dose_response_statistics()` — AIC and BIC from `stats::lm` and `drc::drm` are stored with a warning comment and an `aic_comparison_note` return field: the two likelihoods use different parameterisations and their information criteria are not directly comparable for model selection.
- `analyze_drug_synergy()`, `analyze_drug_synergy_over_time()` — added `@section Assumptions and Limitations:` documenting the Bliss ceiling effect (individual TGI > 50% makes synergy detection nearly impossible) and the Loewe single-dose linear approximation limitation.

## [0.3.4] - 2026-05-14

### Removed
- `post_power_analysis()` — eliminated entirely. Post-hoc power computed from observed effect sizes is uninformative: it is a 1-to-1 monotone function of the p-value and provides no information beyond it (Hoenig & Heisey, 2001). All power and sample-size questions should be answered prospectively with `apriori_power_analysis()` or `apriori_power_simulation()`. The `export(post_power_analysis)` entry has been removed from `NAMESPACE`.

### Fixed
- `body_weight_auc()`, `efficacy_toxicity_bivariate()`, `total_benefit_area()` — removed local `trap_auc()` definitions; all three now call the shared `calculate_auc()` from `utils_auc.R`. Identical implementations had drifted into three separate files.
- `apriori_power_simulation()` — likelihood-ratio test (LRT via `lme4`) is now the preferred method for extracting a p-value from the fitted LMM, with Satterthwaite (`lmerTest`) as a fallback. Previously the branches were inverted: because `lme4` does not populate p-values in the coefficient table, the LRT path was almost never taken and Satterthwaite was always used.
- `analyze_drug_synergy_over_time()` — removed internal Bliss-validation recalculation block that duplicated logic already computed by `analyze_drug_synergy()`; the Bliss Expected TGI is now read directly from the per-timepoint summary. Annotation x-coordinates in the synergy-trend plot now use additive range-based offsets instead of multiplicative offsets, which broke at `Time_Point = 0`. Row accumulation in the timepoint loop replaced with list-based pre-allocation.
- `survival_statistics()` — `print_results()` no longer re-derives median survival by refitting `survfit` internally; it now reads from the already-computed `results` data frame passed by the caller. The previous approach updated only a local copy, leaving `NA` values in the returned result list.
- `dose_response_statistics()` — `verbose` is now forwarded to `generate_summary_statistics()`, so message output is properly suppressed when `verbose = FALSE`. All `if (verbose)` guards standardised to `if (isTRUE(verbose))`.
- `survival_statistics()`, `tumor_growth_statistics()` — all `if (verbose)` guards standardised to `if (isTRUE(verbose))`.
- `tumor_auc_analysis()` — removed dead `exists("plot_auc", mode = "function")` guard and its unreachable fallback branch; `plot_auc()` is always available via the package namespace.

## [0.3.3] - 2026-05-12

### Added
- `apriori_power_analysis()` — analytic a priori power analysis for prospective sample-size planning. Supports two-sample t-test (two groups) and one-way ANOVA via `pwr::pwr.anova.test` (≥3 groups, with non-central F fallback). Accepts Cohen's d directly or a raw mean difference + pooled SD. Two modes: `"find_n"` (required N per group for a target power × alpha grid) and `"find_power"` (achieved power at a specified N). Optionally generates a variability-sensitivity table showing how required N changes as the assumed SD is perturbed by ±20% and ±40%.
- `apriori_power_simulation()` — LMM-based a priori power analysis via Monte Carlo. Simulates longitudinal tumour-growth data with per-mouse random intercepts and slopes, fits `Volume ~ Treatment * Day + (Day|ID)` via `lme4::lmer`, and uses a likelihood-ratio test to determine significance. Power is the proportion of simulations where the Treatment × Day interaction is significant. Supports a `progress_fn` callback for Shiny integration.

## [0.3.2] - 2026-05-05

### Removed
- `efficacy_toxicity_bivariate()` — Log Cell Kill (`log_cell_kill`) efficacy metric removed. The metric required exponential control-group tumor growth and reliable volume observations at a common target size, conditions that are rarely satisfied in practice. The `efficacy_metric` argument now accepts `"tgi"` or `"tumor_auc"` only; the control-group doubling-time computation and `ctrl_doubling_time` return-list element are also gone.

## [0.3.1] - 2026-05-04

### Fixed
- `body_weight_auc()` — now uses a composite `Treatment + ID` key when grouping observations per mouse. Previously, IDs shared across treatment groups (e.g. reused ear tags or cage labels) caused all mice to be attributed to the first treatment group encountered, so AUC results showed only the control group.
- `weight_loss_threshold()` — same composite-key fix applied to the per-mouse baseline computation, event detection, and KM strata. KM curves previously collapsed to a single stratum (control group) when IDs were not unique across treatments.

## [0.3.0] - 2026-04-20

### Added
- `analyze_body_weight()` — longitudinal mixed-effects modeling of body weight
  - Tumor weight adjustment (Volume × density / 1000)
  - Optional covariates: tumor volume, sex, initial mass
  - REML or ML estimation; auto-fallback from random slope+intercept to intercept-only on convergence failure
  - Returns fixed effects, random effects, emmeans, and model diagnostics
- `body_weight_auc()` — trapezoidal AUC of body weight change
  - Per-mouse AUC, % change from baseline, nadir weight analysis
  - Group summaries with pairwise comparisons (Welch t-tests)
- `weight_loss_threshold()` — time-to-threshold survival analysis
  - Kaplan-Meier + log-rank test for time to specified % weight loss
  - Optional Cox PH model; configurable baseline day and threshold
- `therapeutic_window_metric()` — TWM = TGI / MaxWeightLoss%
  - Noise floor: when weight loss ≤ threshold, safety score = TGI
  - Per-group ranking
- `efficacy_toxicity_bivariate()` — safety-efficacy data for bivariate plots
  - Per-mouse and per-group toxicity (max % weight loss) vs efficacy
  - Supports three efficacy metrics: Final TGI, Tumor AUC, Log-Cell Kill
  - Log-Cell Kill uses growth delay formula: LCK = (T − C) / (3.32 × Td)
- `total_benefit_area()` — integrated efficacy-toxicity benefit score
  - B = AUC_efficacy − λ × AUC_toxicity with adjustable λ
  - Per-group benefit ranking
- `weight_corrected_tgi()` — TGI excluding mice exceeding weight loss threshold
  - Compares corrected vs uncorrected TGI per group
  - Reports excluded mouse counts and identities

## [0.2.1] - 2026-04-17

### Changed
- Replaced simplified arithmetic mean synergy model with proper Loewe Additivity (Berenbaum, 1989)
  - Expected fractional effect: `min(FE_A + FE_B, 1.0)` (was `(FE_A + FE_B) / 2`)
  - Combination Index: `(FE_A + FE_B) / FE_combo` (was `(FE_A + FE_B) / (2 * FE_combo)`)
  - CI thresholds: < 0.85 synergistic, 0.85–1.15 additive, > 1.15 antagonistic
- Renamed `analyze_drug_synergy()` output field `$additive_model` to `$loewe_additivity`
- Renamed `analyze_drug_synergy_over_time()` columns `Additive_Mean_Expected_TGI` / `Additive_Mean_Difference` to `Loewe_Expected_TGI` / `Loewe_Difference`
- Updated `synergy_metrics` data frame labels to reference Loewe Additivity

### Fixed
- Fixed Combination Index formula: numerator now uses capped Loewe expected FE `min(FE_A + FE_B, 1.0) / FE_combo` — previously the uncapped sum `(FE_A + FE_B) / FE_combo` inflated CI and falsely indicated antagonism when individual drug effects sum to > 1.0
