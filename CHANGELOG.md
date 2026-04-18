# Changelog

All notable changes to the mouseExperiment package will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.1] - 2026-04-17

### Changed
- Replaced simplified arithmetic mean synergy model with proper Loewe Additivity (Berenbaum, 1989)
  - Expected fractional effect: `min(FE_A + FE_B, 1.0)` (was `(FE_A + FE_B) / 2`)
  - Combination Index: `(FE_A + FE_B) / FE_combo` (was `(FE_A + FE_B) / (2 * FE_combo)`)
  - CI thresholds: < 0.85 synergistic, 0.85–1.15 additive, > 1.15 antagonistic
- Renamed `analyze_drug_synergy()` output field `$additive_model` to `$loewe_additivity`
- Renamed `analyze_drug_synergy_over_time()` columns `Additive_Mean_Expected_TGI` / `Additive_Mean_Difference` to `Loewe_Expected_TGI` / `Loewe_Difference`
- Updated `synergy_metrics` data frame labels to reference Loewe Additivity
