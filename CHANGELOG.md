# Changelog

All notable changes to the mouseExperiment package will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.4.7] - 2026-06-02

Closes CODE_REVIEW.md Round 2 items **D.2**, **D.3**, and **K.10** —
the long-deferred architecture / process items from the v0.3.6 Round 2
review.

### Added

- **`tg_priors()` and `tg_mcmc()` config helpers** (CODE_REVIEW.md D.2).
  New exported helpers that bundle the five prior-related and four
  MCMC-related arguments accepted by every `bayesian_*` entry point.
  Lets callers replace
  `bayesian_tumor_growth(df, prior_strength = "weakly_informative", n_chains = 2, n_iter = 800, ...)`
  with
  `bayesian_tumor_growth(df, priors = tg_priors(strength = "weakly_informative"), mcmc = tg_mcmc(chains = 2, iter = 800), ...)`.
  The change is **fully additive** — every individual `prior_*` /
  `n_chains` / `n_warmup` / `n_iter` / `seed` / `backend` argument
  continues to work exactly as before. When a config object is
  supplied, its fields take precedence over the legacy individual
  arguments. The two helpers live in `R/bayesian_config.R`.
- Wired the `priors` and `mcmc` arguments into
  `bayesian_tumor_growth()`, `bayesian_body_weight()`,
  `bayesian_survival()` (priors$sigma → prior_aux mapping documented),
  and `bayesian_synergy()`. `bayesian_dose_response()` accepts `mcmc`
  only because its Hill-model priors (`prior_emax`, `prior_ec50`,
  `prior_hill`) are model-specific and not bundled by `tg_priors()`.

- **Coverage measurement infrastructure** (CODE_REVIEW.md K.10).
  - `covr (>= 3.6.0)` added to `Suggests`.
  - `coverage.R` script in the package root: `Rscript coverage.R`
    produces a `covr::package_coverage()` summary; `--html` emits an
    HTML report. Bayesian functions are excluded by default (Stan
    compilation dominates the run).
  - Known caveat: the script's first run surfaces pre-existing test
    failures in `test-post_power_analysis.R` (function removed in
    v0.3.4) and `test-toxicity_functions.R` (`efficacy_metric =
    "log_cell_kill"` is no longer a valid choice). Both are stale-test
    issues tracked under CODE_REVIEW.md K.11 and block `covr` from
    producing a clean summary until they're cleaned up.

### Changed — documentation

- **Separating analysis from visualization** (CODE_REVIEW.md D.3).
  Added a new `@section Separating analysis from visualization` to
  `bayesian_tumor_growth()`'s roxygen documentation. Documents the
  pattern that every analysis function in this package already follows:
  data frames (`treatment_effects`, `posterior_summary`,
  `pairwise_comparisons`, …) are returned alongside pre-rendered plot
  objects (`pp_check_plot`, `credible_intervals_plot`, …), and callers
  can pass `plots = FALSE` to skip plot generation entirely and receive
  data-only results suitable for headless / CI pipelines. Verified that
  every `bayesian_*` function guards plot generation with
  `if (isTRUE(plots))`. The Shiny dashboard already uses this pattern.

## [0.4.6] - 2026-06-01

Closes 12 more CODE_REVIEW.md Round 2 items, including every open
"Enhancement" entry and three of the four open Bayesian-power
"Missing" entries. Only the two invasive API refactors (D.2, D.3)
remain deferred at this point.

### Added — Bayesian diagnostics & features
- **PSIS posterior predictive interval coverage** (`bayes_ppc_coverage()`,
  E.6). One-row data frame (`cov_50`, `cov_80`, `cov_95`, `n_obs`) for
  every Bayesian fit. Mis-calibration flags model mis-specification
  cheaply; complements the existing `loo_diagnostics` (out-of-sample).
- **Per-animal posterior growth rates from brms** (E.3). When
  `bayesian_tumor_growth(random_effects_specification = "slope")` is
  used, `growth_rates` is now derived from per-animal posterior random
  slopes (with credible intervals) instead of OLS on log-volumes. New
  internal helper `tg_brms_per_animal_growth_rates()`.
- **Bayesian power: multi-group, random-slope, null calibration**
  (J.5 / J.6 / J.7). `bayesian_power_analysis()` gains:
  - `n_groups`: total groups including control; per-arm
    `treatment_effect` accepted as scalar or vector of length
    `n_groups - 1`. Success criterion is "max P(β < −δ) > target_prob"
    across treated groups.
  - `random_effects_specification = c("intercept_only", "slope")` to
    match the downstream fitting pipeline; new `animal_slope_sd`
    parameter for the slope variance.
  - `null_calibration = TRUE` re-runs the simulation with
    `treatment_effect = 0` and reports the empirical Type-I rate
    (should be ≤ `1 - target_prob`). Power curve plot now overlays
    a dot-dash Type-I line.
- **Percentile bootstrap CIs for AUC pairwise differences** (E.4).
  `tumor_growth_statistics(model_type = "auc")` gains
  `auc_bootstrap_n` and `auc_bootstrap_seed`. When active,
  `pairwise_comparisons` gets `boot_ci_lower` / `boot_ci_upper`
  alongside the parametric Welch's t-CI. New helper
  `tgs_boot_diff_ci()`.
- **cmdstanr backend option** (E.1). All six Bayesian entry points and
  `bayesian_power_analysis()` accept `backend = c("rstan", "cmdstanr")`.
  cmdstanr is 3-10× faster to compile and produces better diagnostics;
  rstan stays the default for out-of-the-box use. New helper
  `resolve_brms_backend()` validates the choice and emits installation
  instructions on missing toolchain rather than silently falling back.

### Changed — Architecture & code organisation
- **`bs_fit_synergy_model()` helper extracts shared synergy fit code**
  (G.6). `bayesian_synergy()` and `bayesian_synergy_over_time()` each
  had ~180 LOC of bit-identical transform / cage-setup / RE-formula /
  priors / brms::brm / diagnostics extraction. Both call sites now
  delegate to the helper; any future synergy-fit change lands in one
  place. No behaviour change.
- **`tgs_path_auc()` extracted to R/tgs_path_auc.R** (D.1 partial).
  The 210-LOC AUC path inside `tumor_growth_statistics()` is now a
  separate file, dropping the main file from 1,389 → 1,197 LOC. Pure
  code-organisation refactor. The two other D-class items (D.2:
  config helpers for 20+ arg signatures, D.3: separating ggplot
  generation from statistical computation) remain open — both break
  the dashboard call sites and need paired upstream design.

### Changed — Documentation
- **`me_result` class scope corrected** (J.8). The docstring claimed
  every analysis function returns an `me_result`; in reality the main
  surface returns bespoke lists and only `repeated_measures_anova()`
  uses the class. Documentation rewritten to accurately describe the
  class as a small optional utility. No code change; no API change.

### Tests — discipline & regression
- **Defensive column-name masks tightened** (K.2). Tests that used
  `if (!is.na(col)) <real assertions>` would have silently passed if
  the column was renamed. `test-tumor_growth_statistics.R` (three
  blocks) now asserts `expect_false(is.na(col))` before the real
  expectations.
- **Test-side warning capture helpers** (K.3). New
  `capture_warnings()` + `expect_no_unexpected_warnings()` in
  `helper-fixtures.R`. Existing tests aren't switched in bulk (most
  use benign `lme4` warning suppression that's still appropriate);
  the helpers are in place for incremental adoption.
- **Parameter-sensitivity regression tests** (K.4). New
  `test-param_sensitivity.R` covers six assertions that the
  silent-ignore bug class catches: `random_effects_specification`,
  `transform`, `auc_bootstrap_n`, `ph_test` presence, `cage_column`
  formula entry, and TWM `TGI < 0` clamping. Locks in the post-Round-2
  behaviour and prevents the entire bug class from recurring.

### Open — still deferred
- D.2 (config-helper grouping for 20+ arg signatures): invasive API
  break across every dashboard call site; needs paired design with the
  dashboard module that consumes the configs.
- D.3 (separate ggplot from stat functions): same scope; affects every
  exported analysis function's return shape.

## [0.4.5] - 2026-05-31

Resolves 29 items from `CODE_REVIEW.md` Round 2 across statistical
validity, Bayesian diagnostics, composite-key consistency, and
documentation. See the individual commits for per-item context.

### Added — Bayesian diagnostics (every brms-based function)
- **NUTS sampler diagnostics** (`make_nuts_diagnostics()`, A.2):
  per-fit count of divergent transitions, max-treedepth hits, and
  the minimum per-chain E-BFMI. The Bayesian path no longer reports
  `Converged = TRUE` while hiding 500 divergent transitions.
- **PSIS-LOO + Pareto-k** (`bayes_loo()`, A.4): one-row summary of
  `elpd_loo / p_loo / looic / n_high_k` plus the per-observation
  Pareto-k vector. The Bayesian counterparts of AIC and Cook's
  distance — both previously absent.
- **bayes_R2 posterior summary** (`bayes_r2_summary()`, C.1): posterior
  median / 95 % CrI of variance explained.
- **ESS / N efficiency** (C.3): `make_mcmc_diagnostics()` now appends
  `ESS_per_draw = Bulk_ESS / total_draws`.
- **Directional posterior probability** for emmeans contrasts
  (`emm_p_direction()`, C.2): `P_direction` column on
  `pairwise_comparisons` for `bayesian_tumor_growth()` /
  `bayesian_body_weight()`. Replaces "does the 95 % CrI exclude zero?"
  with a quantitative posterior statement.
- **Standard diagnostic plots for synergy and TWM** (G.4, G.5): the
  four-plot battery (`pp_check_plot`, `prior_posterior_plot`,
  `mcmc_trace_plot`, plus the existing `posterior_dist_plot` /
  `synergy_plot`) now shipped by `bayesian_synergy()`; the TWM wrapper
  surfaces aliases from the two input models.

### Added — Survival diagnostics
- **Proportional-hazards check via `cox.zph()`** (A.1, J.17). Both
  `survival_statistics()` and `weight_loss_threshold()` now return a
  `ph_test` Schoenfeld-residual test object on the standard-Cox path
  and warn on the console when the global p < 0.05.
- **Concordance / C-index** in `survival_statistics()` (E.2).
- **Firth fallback** in `weight_loss_threshold()` (J.17), matching
  `survival_statistics()`'s coxphf path for complete-separation cases.

### Changed — Statistical validity
- **`tumor_growth_statistics()` `auc_method` argument removed** (A.3).
  It was matched but never propagated; users requesting LOCF silently
  received trapezoidal output. Same silent-ignore bug class as the
  Round 1 1.1 `handle_cage_effects` fix. For LOCF AUC, call
  `tumor_auc_analysis(method = "last_observation")` directly.
- **Cox CI uses `summary(coxph)$conf.int`** (B.2). The magic `1.96`
  literal is gone; CI bounds now come from the proper qnorm-derived
  `lower .95` / `upper .95` columns.
- **`analyze_polynomial_trends()` uses true dose scores** (G.7 — *Major,
  silent wrong answer*). `contr.poly()` is now called with the actual
  numeric dose values as `scores`, so linear / quadratic / cubic trend
  decompositions are on the dose axis rather than the (incorrect) index
  axis. Critical for the common unequally-spaced dose designs.
- **DR model labels separate shape from direction** (G.2): the model
  comparison now stores `dr_model_type ∈ {symmetric, asymmetric}` (the
  LL.4-vs-LL.5 axis the AIC actually selects on) plus a new
  `dr_model_direction ∈ {inhibition, stimulation}` derived from
  `sign(coef(model)["Slope:(Intercept)"])`. The previous labels
  (`inhibition`/`stimulation` on shape selection) were misleading.
- **`bayesian_dose_response()` warns on stimulatory data** (G.3): the
  Hill formula bounds Emax in [0, 1] and the Hill exponent > 0, so the
  curve is inhibition-only by construction. Documented prominently and
  surfaced via `warning()` at fit time when mean TGI is materially
  negative.
- **`analyze_body_weight()` actually uses `cage_column`** (J.11). The
  argument was previously attached to the working frame as `wd$Cage`
  but never appeared in the model formula. The LMM and the fallback
  random-intercept formulas now include `(1 | Cage)` when a cage column
  is supplied. Same bug class as Round 1 1.1.

### Changed — Bayesian fit hygiene
- **`bayesian_synergy()` / `bayesian_synergy_over_time()` no longer
  `suppressWarnings()` around `brms::brm()`** (G.1 — Critical). Stan's
  divergent-transition / max-treedepth / Rhat warnings are exactly the
  signals to surface; the prior session would silently report
  `Converged = TRUE` on broken fits. The new NUTS diagnostics provide
  programmatic detection (A.2) and Stan's warnings provide the
  user-visible second line.
- **`bayes_prior_params()` `stop()`s on unknown `prior_strength`** (B.4).

### Changed — Composite-key consistency
- **`tumor_doubling_time()` / `therapeutic_window_metric()` /
  `efficacy_toxicity_bivariate()` toxicity arm** (J.9 / J.14 / J.18).
  All three switched from `ID`-only or `ID + Treatment` aggregation to
  the canonical `make_mouse_key(ID, Treatment, Cage)` composite key so
  reused IDs across cages or treatments no longer collapse. Same fix
  class as Round 1 1.5 / 1.8.
- **`tgs_extrapolate()` reports the actual count of extrapolated
  subjects** (B.1). Previously
  `length(unique(df$Extrapolated[df$Extrapolated]))` always evaluated
  to 0 or 1.

### Changed — Zero-handling canonicalisation
- **`analyze_growth_rate()` and `repeated_measures_anova()` now use the
  canonical `log(x); x[x<=0] <- min(positive)/2` zero-handling pattern**
  (G.8, J.10), matching `tumor_growth_statistics()` and
  `bayesian_tumor_growth()`. The previous `log1p()` / `log(x + 1)`
  patterns introduced a fixed +1 mm^3 bias on the small early-day
  volumes that matter most.

### Changed — TWM / safety scoring
- **`therapeutic_window_metric()` clamps `TGI < 0` to 0** (J.13). A
  treatment that *accelerated* tumour growth (negative TGI) and caused
  weight loss previously earned a *positive* TWM via `abs(TGI)` — the
  worst case the metric is supposed to flag. Now scores 0.
- **`bayesian_therapeutic_window()` renders the exact TWM = 1 piecewise
  isoline** (J.2). The single-line `geom_abline()` approximation is
  replaced with a four-corner `geom_path()` matching the true
  horizontal-then-linear shape.
- **Independence-assumption caveat** documented in `@section` (J.1):
  the draw-pairing across the two independently-fitted TG/BW posteriors
  ignores within-animal correlation between efficacy and toxicity, so
  the TWM CrI can be too wide or too narrow depending on the sign of
  that correlation.

### Changed — Power / synergy / dose-response
- **`apriori_power_simulation()` `baseline_sd` is no longer a ghost
  parameter** (J.4). It's now converted to a log-scale SD via the
  first-order delta method and added as a per-mouse jitter on
  `log(baseline_volume)`. Setting `baseline_sd = 0` reproduces the
  previous deterministic behaviour.
- **`analyze_drug_synergy()` / `bayesian_synergy()` /
  `bayesian_synergy_over_time()` share `synergy_bliss_expected()` and
  `synergy_loewe_ci()` helpers** (I.1) in `R/utils_synergy.R`. All three
  formulas are scalar/vector polymorphic so a single implementation
  serves point-estimate and posterior-draw call sites.
- **`plot_combination_index()` annotation positions use additive
  offsets** (J.3). Round 1 3.9 fixed `plot_synergy_trend()` for studies
  that don't start at day 0; this finishes the pair.

### Changed — Documentation
- **`apriori_power_analysis()` clarifies effect-size scale** (E.5):
  Cohen's d is on the *modelling scale* (log-volume by default), not
  on raw volume. Adds a worked example of converting log-scale d to
  approximate fold-difference.
- **`bayesian_dose_response()` documents the data-aware EC50 prior**
  (G.9): the prior is centred on `log(median(non_zero_doses))` so it
  stays sensible regardless of dose units.
- **`total_benefit_area()` `lambda` parameter** documents the
  unit-mismatch between TGI-point-days and weight-loss-point-days, and
  the clinically arbitrary equivalence at the `lambda = 1` default
  (J.16).
- **`therapeutic_window_metric()` `noise_floor` default** of 1.0 % is
  now explicitly an experimentally pragmatic threshold for numerical
  stability with no formal clinical basis (I.2).

### Fixed — Utility hygiene
- **`calculate_volume()` `formula` argument is `match.arg`'d** (J.19).
  Unknown values (e.g. the typo `"ellipsoid_3axes"`) previously fell
  through a default-catch branch that silently computed the ellipsoid
  volume. Now fails loudly.
- **`calculate_dates()` `in_place = TRUE` is deprecated** (J.20). The
  `parent.frame()` + `assign()` mechanism silently failed when called
  from inside another function. Matches the calculate_volume()
  deprecation from Round 1 3.3.
- **Removed unused `my_data` example dataset** (J.21). Template
  scaffolding unrelated to the package's purpose; the two roxygen
  examples that referenced it now point at `master_synthetic_data`.
- **`total_benefit_area()` `is.finite()` guard** catches both NaN and
  ±Inf TGI values, not just NaN (J.15).

### Fixed — Tests
- **Stale Bayesian TG / BW / Survival tests updated to current API**
  (K.1): `model_type_used = "bayes_tg"` (not `"bayes"`), CrI columns
  use `Lower_CrI`/`Upper_CrI` (not `Lower_CL`). Tests previously would
  have failed whenever brms was installed; CI skips masked the
  staleness.
- **Bayesian convergence assertions tightened** (K.6): `Rhat <= 1.01`
  (Vehtari 2021) instead of the legacy `< 1.1`, plus new `Bulk_ESS`
  and `Tail_ESS` ≥ 400 assertions across all four Bayesian test files.

### Open — deferred to follow-up sessions
- G.6: extract `bs_fit_synergy_model()` helper to dedupe the
  ~300 LOC shared between `bayesian_synergy()` and
  `bayesian_synergy_over_time()`.
- J.5 / J.6 / J.7: multi-group support, random-slope option, and
  null-distribution / type-I calibration in `bayesian_power_analysis()`.
- J.8: decide whether to wire the analysis functions into the `me_result`
  S3 class (currently defined but unused by main analysis functions)
  or delete the class entirely.
- K.2 / K.3: replace defensive `if (!is.na(col))` masking + blanket
  `suppressWarnings()` patterns in the test suite.
- K.4: per-function parameter-sensitivity tests (the missing safety net
  that would have caught Round 1 1.1, A.3, J.11).
- D.1 / D.2 / D.3: file-size refactors (1000+ LOC `tumor_growth_statistics.R`,
  `bayesian_synergy.R`, `bayesian_body_weight.R`); param-group config
  helpers for 20+ arg signatures; separating ggplot generation from
  statistical computation.
- E.1: `cmdstanr` backend option across all Bayesian functions.
- E.3: per-animal posterior growth rates drawn from the brms model itself
  instead of via OLS on log-volumes.
- E.4: bootstrap CIs for AUC treatment effects.
- E.6: posterior predictive interval coverage check.

## [0.4.4] - 2026-05-30

### Added
- **GAM model option for tumor growth, body weight, and Bayesian tumor growth.**
  Replaces the previous mechanism of fitting higher-order polynomials in time:
  - `tumor_growth_statistics()` — new `model_type = "gam"` choice. Fits via
    `gamm4` (lme4 random effects + mgcv smoothers) with a group-specific
    smoother on Day: `y ~ Treatment + s(Day, by = Treatment, k)`. Smoother
    basis dimension is auto-chosen from the number of unique time points
    (clamped to `[3, 10]`). Returns the same result-list shape as the LME4
    path — `treatment_effects` (at mean day), `treatment_effects_over_time`
    (at the five study-day quantiles), `pairwise_comparisons` (smooth-vs-
    smooth differences at quantile days, pairwise vs the reference group),
    `anova` (mgcv smooth-term significance), and a `diagnostics` block
    including `k_check` (basis adequacy) and `deviance_explained`.
  - `bayesian_tumor_growth()` — new `model_type = c("lmm", "gam")` argument
    (default `"lmm"` preserves prior behaviour). When `"gam"`, brms fits
    `y ~ Treatment + s(Day, by = Treatment, k)` with the same random-effects
    spec. The returned `model_type_used` is `"bayes_tg_gam"`.
  - `analyze_body_weight()` — new `model_type = c("lmm", "gam")` argument.
    GAM path uses `gamm4` with `(1 | Cage) + (1 | ID)` random effects when a
    cage column is supplied.
- `R/tgs_gam.R` — new internal helpers shared across these paths:
  `tgs_fit_gamm4_model()`, `tgs_gam_treatment_effects()`, `tgs_gam_emm_time()`,
  `tgs_gam_pairwise()`, `tgs_gam_anova_table()`, `tgs_gam_diagnostics()`.

### Removed
- **`polynomial_degree` argument** removed from `tumor_growth_statistics()`
  entirely. Higher-order polynomial-in-time fits were a workaround for
  non-linear growth; `model_type = "gam"` is the principled replacement and
  extrapolates more honestly. Existing callers that previously passed
  `polynomial_degree = 1` are unaffected (linear time is still the default
  LME4 behaviour). Callers passing `polynomial_degree > 1` should migrate
  to `model_type = "gam"`.

### Dependencies
- `gamm4 (>= 0.2-6)` and `mgcv (>= 1.8-40)` added to `Suggests`. Required
  only when `model_type = "gam"` is requested.

## [0.4.3] - 2026-05-20

### Added
- `R/utils_bayes.R` — new internal module centralising shared Bayesian helpers
  (B2.1–B2.4, B4.4):
  - `bayes_backtransform(x, transform)` — single back-transform helper replacing
    three identical `bt()` closures in `bayesian_body_weight`, `bayesian_synergy`,
    and `bayesian_therapeutic_window`.
  - `make_mcmc_diagnostics(posterior_summary_df)` — standardised MCMC diagnostic
    data frame builder (Rhat, Bulk_ESS, Tail_ESS, Converged); replaces duplicated
    inline construction in five functions.
  - `bayes_prior_params(prior_strength)` — centralised prior-hyperparameter lookup
    (`b_sd`, `exp_rate`) replacing five independent `switch()` blocks.
  - `setup_cage_column(df, cage_column)` — cage placeholder setup (returns
    `list(df, cage_column, no_cage_mode)`); replaces duplicated inline logic in
    four functions.
  - `build_posterior_summary(model)` — standardised `summary(model)$fixed`
    extractor with `Lower_95_CrI`/`Upper_95_CrI` column renaming.
  - `bayes_prior_posterior_plot(model, treatment_column)` — moved from
    `bayesian_tumor_growth.R`; prior vs posterior density overlay shared by all
    Bayesian functions.
- `bayesian_power_analysis()` — simulation-based Bayesian a priori power analysis
  (B6.1). Estimates Bayesian power (fraction of simulated experiments where
  `P(β_trt:Day < -δ | data) > target_prob`) across a grid of sample sizes.
  Uses `brms::update()` to reuse compiled Stan programs across simulations.
  Returns `power_table`, `power_curve_data`, `params`, and `power_curve_plot`.
- `bayesian_twm_from_data()` — single-call convenience wrapper for the
  Therapeutic Window Metric (B6.2). Internally fits `bayesian_tumor_growth()`
  and `bayesian_body_weight()` then passes results to
  `bayesian_therapeutic_window()`. Returns the full TWM result list plus
  `tg_result` and `bw_result` sub-model outputs.
- `bayesian_synergy_over_time()` — time-series extension of `bayesian_synergy()`
  (B6.3). Fits the same `Treatment × Day` brms model then evaluates draw-wise
  Bliss excess and Loewe CI at every study day via a single
  `posterior_epred()` call on the full day-by-group grid. Returns
  `synergy_by_day`, `tgi_by_day`, `peak_bliss_day`, `peak_loewe_day`, and
  `synergy_time_plot` (faceted ribbon plot).

### Fixed
- **B1.1** `model_type_used` in `bayesian_tumor_growth()` changed from `"bayes"`
  to `"bayes_tg"` for consistency with `"bayes_bw"`, `"bayes_survival"`, etc.
- **B1.2** `prior_strength` default in `bayesian_body_weight()` aligned to
  `"skeptical"` (was `"weakly_informative"`).
- **B1.3** `"informative"` prior preset added to `bayesian_dose_response()`.
- **B1.4** `Lower_CL`/`Upper_CL` renamed to `Lower_CrI`/`Upper_CrI` throughout
  all Bayesian function outputs (`bayesian_tumor_growth`, `bayesian_body_weight`,
  `bayesian_survival`). "CL" (confidence limit) is frequentist terminology;
  "CrI" (credible interval) is correct for Bayesian posteriors.
- **B1.5** `include_frailty` parameter in `bayesian_survival()` renamed to
  `include_cage_effect` for consistency with all other Bayesian functions.
- **B1.6** Intercept prior width in `bayesian_synergy()` corrected from
  `b_sd * 2.0` to `b_sd * 2.5` to match all other Bayesian functions.
- **B1.7** `set_prior()` in `bayesian_synergy()` replaced with `prior_string()`
  (matching the brms API used everywhere else).
- **B3.1** Loewe floor in `bayesian_synergy()` changed from a fixed `1e-6` to a
  data-relative `max(max(FE_combo) * 1e-4, 1e-4)`. `Floor_Applied` flag added to
  `loewe_summary`.
- **B3.2** `bayesian_synergy()` MCMC diagnostics now include Bulk_ESS and
  Tail_ESS (previously missing; only Rhat was reported via `neff_ratio`).
- **B3.3** All-zero/negative volume guard added before log transform in
  `bayesian_tumor_growth()`, `bayesian_body_weight()`, and `bayesian_synergy()`.
  Stops with an informative error instead of silently producing `log(Inf)`.
- **B3.4** TWM group-name intersection in `bayesian_therapeutic_window()`
  normalised with `trimws()` to prevent mismatches from leading/trailing
  whitespace.
- **B3.5** `bliss_summary` rounding standardised to 3 dp throughout
  `bayesian_synergy()` (was inconsistently 3 vs 4 dp).
- **B4.1** `@param reference_group` in `bayesian_tumor_growth()` now documents
  the alphabetical auto-detection strategy.
- **B4.2** `@return` for `bayesian_body_weight()` now documents that weight-loss
  percentages use the earliest/latest day in the model data and are
  population-level predictions.
- **B4.3** `@section Assumptions and Limitations:` added to `bayesian_synergy()`
  documenting the Bliss ceiling effect and Loewe linear-dose approximation.
- **B5.1** `bayesian_body_weight()` endpoint weight-loss computation reduced from
  three `posterior_epred()` calls to two by batching first-day and last-day
  prediction grids into a single call.

## [0.4.2] - 2026-05-20

### Added
- `bayesian_synergy()` — new exported function for Bayesian drug combination synergy analysis with full posterior uncertainty quantification. Fits a single Bayesian linear mixed-effects model (Treatment × Day) on all four groups (control, drug A, drug B, combination) using `brms`, then propagates posterior uncertainty through both Bliss Independence and Loewe Combination Index metrics via draw-wise arithmetic on `posterior_epred` samples. Key features:
  - `tgi_summary`: per-group posterior median TGI with 95 % CrI at the endpoint day.
  - `bliss_summary`: posterior distribution of Bliss excess (Δ = observed FE_combo − Bliss expected FE); `P_Synergy = P(Δ > 0)` directly interpretable as the probability that the combination exceeds Bliss-independence additivity. Includes `Expected_FE_Median` and `Observed_FE_Median` for direct comparison.
  - `loewe_summary`: posterior distribution of Loewe CI (single-dose approximation: CI = min(FE_A + FE_B, 1) / max(FE_combo, ε)); `P_Synergy = P(CI < 1)`; `Interpretation` string ("Synergistic", "Additive", or "Antagonistic") based on the posterior-median CI.
  - `synergy_table`: combined six-row summary (four observed groups + Bliss-expected + Loewe-expected) with posterior median volumes, TGI %, and 95 % CrI.
  - `synergy_plot`: bar chart of observed TGI per group with 95 % CrI error bars; dashed blue line = Bliss expected TGI; dotted red line = Loewe expected TGI.
  - `posterior_dist_plot`: side-by-side density overlays of Bliss excess and Loewe CI draw distributions with reference lines at 0 / 1 respectively and P(synergy) annotation.
  - Same prior-strength presets, transform options, and cage random-intercept support as `bayesian_tumor_growth()`. Returns `model_type_used = "bayes_synergy"`, `transform_used`, `posterior_summary`, and `mcmc_diagnostics`.
- Tests for `bayesian_synergy()` in `tests/testthat/test-bayesian_synergy.R` (22 tests): return structure, `model_type_used`, `transform_used`, `brmsfit` class, `tgi_summary` row count and columns, control TGI near zero, treated groups positive TGI, `bliss_summary` fields, Bliss P_Synergy ∈ [0,1], CrI brackets median, `loewe_summary` fields, Loewe CI positive and finite, Loewe P_Synergy ∈ [0,1], Interpretation string, `synergy_table` row count and Type column, MCMC diagnostics, summary metadata, `synergy_plot`/`posterior_dist_plot` ggplot class, input-validation errors, `return_model = FALSE`, `plots = FALSE`.

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
