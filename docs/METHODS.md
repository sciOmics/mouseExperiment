# Statistical methods reference

For each analysis module: what the method does, when to choose it, what the result fields mean, and the assumptions you accept by using it.

This document complements the per-function roxygen — read it when you're trying to decide which model to use, not when you're looking up a specific argument. The roxygen is canonical for argument-level detail; this is canonical for the "why" and the "what's in the box".

For the Bayesian-specific diagnostics surface (Rhat, ESS, NUTS, LOO, Bayes R², PPC coverage, etc.), see [`BAYESIAN.md`](BAYESIAN.md).

---

## Tumor Growth

Function family: `tumor_growth_statistics()`, `bayesian_tumor_growth()`, `tumor_auc_analysis()`, `tumor_doubling_time()`.

The tumor growth pipeline supports four model types via `model_type`:

| Model | Function | Use when |
|---|---|---|
| `"lme4"` | `tumor_growth_statistics(..., model_type = "lme4")` | Default. Linear mixed-effects on log-volume; per-animal random intercept (and optionally slope). Comparable to common preclinical-oncology defaults. |
| `"auc"` | `tumor_growth_statistics(..., model_type = "auc")` | When trajectories are non-monotonic or the time-window matters more than the per-mouse slope. Per-mouse AUC under the growth curve; pairwise treatment comparisons via Welch's t-test. |
| `"gam"` | `tumor_growth_statistics(..., model_type = "gam")` | When growth is clearly non-linear (e.g., flattens then resumes after relapse). Penalized splines per group via `mgcv`. |
| Bayesian LMM | `bayesian_tumor_growth(...)` | When you want posterior probability statements, credible intervals with direct probability interpretation, or you have small N where the prior matters. Slower (3-12 min); use `cmdstanr` backend for faster compile. |

### Required columns

| Column | Notes |
|---|---|
| `id_column` | One identifier per animal (no duplicates across timepoints in the same group) |
| `time_column` | Numeric (day) or `Date` — date detection auto-converts |
| `volume_column` | Tumor volume (mm³). If you have `Length`/`Width`, run `calculate_volume()` first |
| `treatment_column` | Treatment group label |

Optional: `cage_column` (random intercept), `dose_column` (when crossed with treatment), `necrotic_column` (binary indicator → handled per `necrotic_handling`).

### Transforms (`transform` argument)

- `"log"` (default) — exponential growth is linear on log scale; standard.
- `"sqrt"` — variance-stabilizing for Poisson-like growth; rarely the right choice for tumors but available.
- `"none"` — for already-transformed data or when explicitly testing additive effects.

### LME4 path: what the result list contains

`tumor_growth_statistics(model_type = "lme4")` returns:

| Field | Description |
|---|---|
| `model` | The fitted `lmerMod` object (or `NULL` when `return_model = FALSE`) |
| `model_type_used` | `"lme4"` |
| `transform_used` | The transform that was applied |
| `treatment_effects` | One row per treatment group: `Adjusted_Mean`, `SE`, `Lower_CL`, `Upper_CL` (95% CI). Marginal means from `emmeans` at the mean study day |
| `pairwise_comparisons` | Treatment-pair contrasts. Columns: `contrast`, `estimate`, `SE`, `df`, `p.value`, `lower.CL`, `upper.CL` |
| `all_pairwise_comparisons` | All pairs, not Dunnett-filtered (used by the dashboard's forest plot) |
| `reference_comparisons_dunnett` | Only present when `reference_group` is specified; Dunnett-style multiplicity-adjusted contrasts vs reference |
| `anova` | Type II/III ANOVA from `car::Anova` |
| `growth_rates` | Per-animal exponential growth rate (= slope on log scale); one row per animal |
| `data_summary` | N per group, mean baseline volume, mean final volume |
| `tumor_growth_plot` | ggplot of trajectories with group means overlaid (only when `plots = TRUE`) |
| `diag_qq_plot`, `diag_resid_fitted_plot` | Diagnostic plots (only when `include_diagnostics = TRUE`) |

### AUC path

`model_type = "auc"` returns the same shape but with `auc_analysis` populated and `growth_rates` may be NULL. `auc_analysis$individual` is one row per animal with the integrated AUC; `auc_analysis$summary` is per-group mean ± SD. Pairwise comparisons use Welch's t on log-AUC.

### Assumptions you're accepting

- **LME4:** linearity on transformed scale; independent residuals within mouse (the random intercept absorbs mouse-level variation); homoscedasticity across treatments
- **AUC:** unequal variances allowed (Welch); independence between animals; the AUC is a meaningful summary of growth (often arguable on its own — pair with growth rates)
- **GAM:** smoothness penalty correctly trades fit vs wiggliness (the default `mgcv::s()` basis usually does fine; raise the basis dimension only if the smooth looks under-fitted)
- **Bayesian LMM:** prior choice is honest; see [`BAYESIAN.md`](BAYESIAN.md)

### When to pick what

- **Default:** LME4. Familiar; testable; fast.
- **Small N (< 4 per group)**: Bayesian LMM with a `weakly_informative` or `informative` prior. Frequentist CIs at small N are notoriously narrow / over-confident; the Bayesian path is more honest.
- **Non-exponential growth (flattening late):** GAM. The smooth captures the late-time shape.
- **Direct comparison of total tumor burden:** AUC. Use it as a confirmatory metric alongside the rate-based model, not as a replacement.

---

## Survival

Function family: `survival_statistics()`, `bayesian_survival()`.

### `survival_statistics()` — frequentist

Returns:

| Field | Description |
|---|---|
| `km_fit` | `survival::survfit` object (Kaplan-Meier) |
| `cox_model` | `survival::coxph` model |
| `log_rank_test` | `survival::survdiff` chi-squared test |
| `ph_test_table` | `cox.zph()` proportional-hazards check (one row per term + global). Reject the null (p < 0.05) means PH is violated for that term |
| `concordance` | C-index (Harrell's). 0.5 = random; 1.0 = perfect; report alongside HRs |
| `results` | Treatment-effect summary table |
| `survival_data` | The prepared data (Time, Event, Treatment) |

Note: when PH is violated, the Cox HR is no longer a constant time-ratio — it's an average over the time window. Switch to the Bayesian AFT path or stratify the Cox model in that case.

### `bayesian_survival()` — Bayesian AFT

Four parametric families:

| Family | Use when |
|---|---|
| `"weibull"` (default) | Flexible — handles increasing, decreasing, and constant hazards. Reportable as both AFT (time ratio) and PH (hazard ratio). |
| `"lognormal"` | Hazard rises then falls. Biologically plausible for treated tumors where late-time animals tend to survive longer once they've cleared an initial dose |
| `"exponential"` | Memoryless / constant hazard. Special case of Weibull with shape = 1. Use only when you have a strong reason to believe the hazard is constant |
| `"gamma"` | More flexible tail than Weibull. Good alternative when Weibull's diagnostics don't fit |

Returns the standard `treatment_effects` shape (Group, Time_Ratio, Lower_CrI, Upper_CrI, HR, Median_Survival, Events, Total, Event_Rate, Note) plus the Bayesian diagnostics block (Rhat, ESS, NUTS, LOO).

`include_cage_effect = TRUE` adds a cage-level frailty `(1 | cage)`. Worth using when you have ≥ 5 cages per group and reason to believe between-cage variation is non-trivial.

### Assumptions

- **Cox PH:** the hazard ratio is constant in time. Violated → switch to Bayesian AFT
- **Weibull AFT:** time on the log scale is linear in the covariates; the shape parameter is constant
- **All survival models:** non-informative censoring (animals censored for reasons unrelated to their hazard)

---

## Body Weight / Toxicity

Function family: `analyze_body_weight()`, `bayesian_body_weight()`, `body_weight_auc()`, `weight_corrected_tgi()`, `efficacy_toxicity_bivariate()`, `total_benefit_area()`, `therapeutic_window_metric()`, `bayesian_therapeutic_window()`, `bayesian_twm_from_data()`.

### `analyze_body_weight()` — frequentist LME4

Linear mixed-effects model with random intercept per animal (and optional cage random effect). Computes:

- Treatment-time interaction (does weight diverge over time?)
- Per-group adjusted means via `emmeans`
- Optionally `adjust_tumor_weight = TRUE`: subtracts estimated tumor weight (volume × `tumor_density`) before modelling

Returns the same `treatment_effects` shape as `tumor_growth_statistics()` plus a `weight_trajectory_plot`.

### `body_weight_auc()`

AUC of weight loss as % from baseline. Lower AUC = better tolerated.

### `weight_loss_threshold` (event derivation)

For survival-style analysis: time to ≥ X% weight loss. Wired into the dashboard's Survival module when "event source" is set to "weight loss".

### `weight_corrected_tgi()`

Tumor Growth Inhibition with weight loss correction:

```
TGI_uncorrected = 1 - (mean_treated / mean_control)   on volume scale
TGI_corrected   = TGI_uncorrected × penalty(weight_loss)
```

The penalty grows non-linearly with weight loss — small losses (≤ 5%) have minimal effect; losses above the safety threshold (default 20%) penalize the TGI heavily.

### `efficacy_toxicity_bivariate()`

Joint analysis of efficacy (TGI or tumor AUC) and toxicity (max weight loss %) per mouse. Output includes per-mouse, per-group means, and a scatter `safety_efficacy_scatter` plot.

### `total_benefit_area()`

Total benefit = ∫(efficacy(t) − λ × toxicity(t)) dt over the study window. The lambda is a clinician-set tradeoff weight (default 1.0).

### Therapeutic Window Metric

```
TWM = TGI(%) / mean_weight_loss(%)
```

Higher TWM = better therapeutic window. The Bayesian variants (`bayesian_therapeutic_window()`, `bayesian_twm_from_data()`) take posterior draws from the underlying TG + BW models and combine them draw-by-draw to propagate uncertainty.

### Assumptions

- **BW LME4:** linearity on natural weight scale; per-mouse random effect captures individual differences
- **Adjusted weight:** the tumor density assumption (default 1.0 g/cm³) — most tumors are within 10% of this; large discrepancies matter
- **TGI:** the comparison is on the chosen endpoint day; trajectories with different shapes can have the same TGI at one day and differ on another
- **TWM:** independence between TG and BW models (current limitation; documented in `bayesian_therapeutic_window()` `@details`)

---

## Drug Synergy

Function family: `analyze_drug_synergy()`, `analyze_drug_synergy_over_time()`, `bayesian_synergy()`, `bayesian_synergy_over_time()`.

Two synergy frameworks supported:

### Bliss Independence

Tests whether the combo effect exceeds what's expected from independent action of each monotherapy:

```
expected_combo = effect_A + effect_B - effect_A × effect_B   (on fractional scale)
synergy        = observed_combo - expected_combo
synergy > 0    = supra-additive (Bliss synergy)
synergy < 0    = sub-additive (Bliss antagonism)
synergy ≈ 0    = independent
```

Reported as `bliss_summary` with the synergy score + 95% CI (or CrI for the Bayesian version).

### Loewe Combination Index

```
CI = (concentration_A_in_combo / IC50_A) + (concentration_B_in_combo / IC50_B)
CI < 1 = synergy
CI ≈ 1 = additive
CI > 1 = antagonism
```

Reported as `loewe_summary`. Requires per-mouse concentration estimates — fragile when monotherapy is barely effective at the tested dose.

### `_over_time` variants

Same metric, computed at each timepoint, returning a `synergy_summary` table with one row per `(day, treatment_pair)`. Useful for detecting late-onset synergy or synergy that decays.

### Bayesian path

`bayesian_synergy()` fits a Bayesian LMM on the four-arm design (A / B / combo / control), then derives Bliss and Loewe scores from the posterior. Output includes `bliss_summary` and `loewe_summary` with median + 95% CrI.

Known caveat (CODE_REVIEW.md backend G.1): an earlier version wrapped `brms::brm()` in `suppressWarnings()`, hiding divergent-transition warnings. The current version exposes diagnostics via the standard `mcmc_diagnostics` + `nuts_diagnostics` fields.

### Assumptions

- **Bliss:** the two drugs act independently (no shared targets / pathways). When the drugs hit the same pathway, "Bliss synergy" can look high but is mechanistically expected
- **Loewe:** dose-response curves are well-approximated by Hill / Emax; the IC50 estimates are reliable

---

## Dose-Response

Function family: `dose_response_statistics()`, `bayesian_dose_response()`.

### Frequentist (`drc` backend)

Fits a 4-parameter logistic (Hill / Emax) via `drc::drm`:

```
volume(dose) = lower + (upper - lower) / (1 + (dose / EC50)^slope)
```

Returns: `linear_model`, `anova_model`, `statistics$ec50`, `hill_slope`, `lower_limit`, `upper_limit`, `growth_dose_p_value`, plus the analysis data.

### Bayesian

Same Hill / Emax form fit via `brms::brm()` with model-specific priors (`prior_emax`, `prior_ec50`, `prior_hill`, `prior_sigma`). Note: these priors are NOT bundled by `tg_priors()` because they're model-specific; only `mcmc = tg_mcmc()` applies here.

Returns posterior summaries for each Hill parameter + `dose_response_curve_plot` with credible bands.

### Assumptions

- Smooth, monotonic dose-response (no biphasic / hormetic effects)
- Independence of measurements across animals
- The Hill equation is the right functional form (it usually is for cytotoxic response; check the diagnostic fit plot)

---

## Power Analysis

Function family: `apriori_power_analysis()`, `apriori_power_simulation()`, `bayesian_power_analysis()`.

### Analytic (`apriori_power_analysis()`)

Standard t-test / one-way ANOVA power via `pwr::pwr.t.test` and `pwr::pwr.anova.test`. Fast (< 1 second). Use when:
- You have an effect-size estimate from a prior study
- The downstream analysis is a simple t-test or one-way ANOVA
- You're sketching the study design — first-pass numbers

### LMM simulation (`apriori_power_simulation()`)

Simulates data under an LME4 model with specified variance structure, fits, and counts the proportion of simulated runs where the treatment effect is significant. Slower (1-10 min depending on `n_sim`). Use when:
- Downstream analysis is the LMM path (you almost always want this for repeated-measures designs)
- You have an estimate of within-mouse variance from prior data
- You need to size for the right test (not the t-test approximation)

### Bayesian simulation (`bayesian_power_analysis()`)

Same idea as the LMM simulation but with Bayesian fitting at each simulated dataset. Reports posterior P(effect > 0) instead of frequentist power. Use when:
- The actual analysis will be Bayesian (consistency)
- You're sizing under genuine prior uncertainty

### Assumptions

- Effect size is on the same scale you'll analyze on (log volume for `transform = "log"`; raw volume for `"none"`)
- The simulated variance structure matches reality (the biggest assumption — sensitivity-analyze across plausible variance values)

---

## Utilities

| Function | Notes |
|---|---|
| `calculate_volume()` | Tumor volume from length × width. Default is the standard ellipsoid `V = (L × W²) × π / 6`. Other formulas: modified ellipsoid, cylinder, sphere. Dimensions are auto-corrected so the longer measurement is always L (handles "I measured them in the wrong order" data) |
| `calculate_dates()` | Date-to-day conversion using a reference date. Two methods: `"direct"` (study day already in data) and `"computed"` (compute from a date column) |
| `calculate_auc()` | Trapezoidal AUC over a time-volume vector |
| `repeated_measures_anova()` | Simple wrapper for rmANOVA via `car::Anova` |
| `export_diagnostics()` | Save diagnostic plots from an analysis result to disk |
| `tumor_doubling_time()` | Doubling time from exponential growth rate (returns `log(2) / rate`) |
| `tg_priors()`, `tg_mcmc()` | Config helpers (v0.4.7+) — bundle prior and MCMC arguments. See [`BAYESIAN.md`](BAYESIAN.md) |

---

## Choosing between frequentist and Bayesian

A short decision guide:

| Situation | Pick |
|---|---|
| N ≥ 6 per group, normal-ish residuals, simple comparison | Frequentist |
| N < 5 per group | Bayesian with `weakly_informative` or `informative` prior |
| You want "P(effect > 0)" or direct probability statements | Bayesian |
| You want to formally incorporate prior data | Bayesian |
| You need to publish in a venue that doesn't yet accept Bayesian methods | Frequentist (or both — Bayesian as a supplementary analysis) |
| You're short on compute (CI runs, exploratory analysis) | Frequentist |
| Frequentist CIs feel too narrow for the data you have | Bayesian (the discrepancy usually means you're under-N; the Bayesian prior makes uncertainty honest) |
| The frequentist model failed to converge | Bayesian (better-behaved when N is small) |

When in doubt, run both. The frequentist result tells you what a colleague reviewing the study will compute; the Bayesian result tells you what you actually believe given the data.

---

## See also

- [`BAYESIAN.md`](BAYESIAN.md) — diagnostics, priors, MCMC interpretation
- Roxygen `?<function_name>` — argument-level reference
- Vignettes in `vignettes/` — worked end-to-end examples
- Dashboard `docs/ARCHITECTURE.md` — how the UI exposes these methods
- `CODE_REVIEW.md` (repo root) — Round 2 audit; explains many design choices
