# Changelog

All notable changes to the mouseExperiment package will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
