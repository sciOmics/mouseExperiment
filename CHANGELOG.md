# Changelog

All notable changes to the mouseExperiment package will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
