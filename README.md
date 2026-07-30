# mouseExperiment

An R package for statistical analysis of mouse tumor growth experiments. Covers the full preclinical-oncology analysis pipeline: tumor growth, survival, body weight toxicity, drug synergy, dose-response, power analysis, and Bayesian equivalents of every major model.

## Features

- **Tumor growth** — LME4 LMM, AUC, GAMM, Bayesian LMM (brms)
- **Survival** — Kaplan-Meier, Cox PH, log-rank, `cox.zph()` PH check, Bayesian AFT (Weibull / log-normal / exponential / gamma), C-index
- **Body weight / toxicity** — LME4 mixed model, Bayesian LMM, weight-loss threshold events, weight-corrected TGI
- **Drug synergy** — Bliss Independence, Loewe Combination Index, Bayesian synergy, Bayesian synergy over time
- **Dose-response** — Frequentist + Bayesian Hill / Emax curve fitting
- **Therapeutic window** — TWM from frequentist or Bayesian TG + BW models; single-call `bayesian_twm_from_data()` wrapper
- **Power analysis** — Analytic (t-test / ANOVA), LMM simulation, Bayesian simulation; multiplicity- and attrition-aware (`n_comparisons`, `dropout_rate`)
- **Randomisation tests** — `trajectory_permutation_test()` tests the treatment × time interaction without the denominator-df or normality approximations. You declare the unit of randomisation via `perm_spec(unit = "mouse" | "cage")`, because a permutation test is only valid if it mirrors how the study randomised; small designs are enumerated exhaustively for an exact p-value, and the design's resolution floor is reported so a null result cannot be over-read
- **Bayesian diagnostics** — Rhat, ESS, NUTS divergences / max_treedepth / E-BFMI, Bayes R², PPC coverage, PSIS-LOO with Pareto-k, posterior P(effect ≠ 0)
- **Comprehensive plots** — KM curves, growth trajectories, synergy bar charts, dose-response curves, forest plots, MCMC diagnostics
- **Config helpers** — `tg_priors()` and `tg_mcmc()` bundle prior + MCMC arguments so Bayesian entry-point signatures stay readable
- **Code-coverage** — `covr::package_coverage()` baseline; one-line `Rscript coverage.R` for HTML / stdout reports

## Current status

| Item | State |
|---|---|
| Version | 0.11.0 |
| `CODE_REVIEW.md` | Rounds 1–5 complete; Rounds 3–5 fully closed |
| Bayesian diagnostics surface | Rhat / ESS / NUTS / Bayes R² / PPC coverage / LOO / Pareto-k / posterior P direction |
| Test suite | testthat, 644 tests. **No test skips a required dependency** — the Bayesian and permutation paths always run (see below) |
| Coverage measurement | `Rscript coverage.R` |

### A note on dependencies

As of v0.10.0 the statistical packages (`brms`, `bayesplot`, `gamm4`, `mgcv`,
`pwr`, `coin`, `clinfun`, `ggpubr`, `posterior`) are **required**, not suggested.
Installing therefore needs a working C++ toolchain, because `brms` pulls the Stan
stack.

That cost is deliberate. While `brms` sat in `Suggests`, the Bayesian tests skipped
wholesale whenever it was absent — and two Critical defects survived five releases
behind that skip, including `bayesian_synergy()` being entirely non-functional from
v0.4.6 to v0.9.0 while this README advertised it. Optionality was not free; it was
the mechanism. See `CODE_REVIEW.md` §R3-L and §R4.

`cmdstanr` remains genuinely optional — it is a user-selected *alternative* Stan
backend and fails loudly with install instructions when chosen and missing.

## Installation

```r
# Development version from GitHub
devtools::install_github("sciOmics/mouseExperiment")

# For Bayesian analyses also install:
install.packages(c("brms", "bayesplot"))

# Optional: cmdstanr backend for 3-10x faster Stan compilation
install.packages("cmdstanr",
                 repos = c("https://mc-stan.org/r-packages/",
                           getOption("repos")))
cmdstanr::install_cmdstan()
```

## Usage

### Tumor Growth — frequentist + Bayesian

```r
library(mouseExperiment)

data(combo_treatment_synthetic_data)
df <- combo_treatment_synthetic_data

# ---- Frequentist LME4 ----
results <- tumor_growth_statistics(
  df               = df,
  time_column      = "Day",
  volume_column    = "Volume",
  id_column        = "ID",
  treatment_column = "Treatment",
  cage_column      = "Cage",
  model_type       = "lme4"
)
print(results$treatment_effects)
plot_tumor_growth(results$tumor_growth_plot)

# ---- Bayesian LMM (recommended config-helper style, v0.4.7+) ----
bayes_results <- bayesian_tumor_growth(
  df     = df,
  priors = tg_priors(strength = "weakly_informative"),
  mcmc   = tg_mcmc(chains = 4, warmup = 1000, iter = 500, seed = 42))

print(bayes_results$treatment_effects)     # group EMMs with 95 % CrI
print(bayes_results$mcmc_diagnostics)      # Rhat, ESS, Converged flag
print(bayes_results$bayes_R2)              # Bayes R² + 95 % CrI
print(bayes_results$loo_diagnostics)       # elpd_loo, n_high_k (Pareto-k > 0.7)
bayes_results$pp_check_plot                # posterior predictive density overlay
bayes_results$prior_posterior_plot         # how much the data updated each prior

# The individual `prior_strength`, `n_chains`, `n_iter`, etc. arguments
# still work as before — supply either the config helpers or the
# individual args, not both.
```

### Data-only path (analysis without ggplot generation)

Every analysis function accepts `plots = FALSE` to skip plot construction and return only the data frames. Useful for headless CI or when callers want to render their own plots from the data via the package's `plot_*()` helpers (or their own renderers). CODE_REVIEW.md D.3.

```r
data_only <- bayesian_tumor_growth(df = df, plots = FALSE,
                                   priors = tg_priors(),
                                   mcmc   = tg_mcmc(chains = 2, iter = 200))
data_only$treatment_effects            # populated
data_only$pp_check_plot                # NULL
data_only$prior_posterior_plot         # NULL
```

### Survival Analysis

```r
data(combo_treatment_synthetic_data)
df <- combo_treatment_synthetic_data

# Frequentist (KM + Cox + log-rank + cox.zph PH check)
surv_results <- survival_statistics(
  df               = df,
  time_column      = "Day",
  censor_column    = "Survival_Censor",
  treatment_column = "Treatment",
  id_column        = "ID"
)
print(surv_results$km_fit)
print(surv_results$ph_test_table)         # cox.zph() output
print(surv_results$concordance)           # C-index

# Bayesian AFT
bayes_surv <- bayesian_survival(
  df               = df,
  time_column      = "Day",
  event_column     = "Survival_Censor",
  treatment_column = "Treatment",
  id_column        = "ID",
  family           = "weibull",
  priors           = tg_priors(strength = "weakly_informative"),
  mcmc             = tg_mcmc(chains = 4, iter = 500))
print(bayes_surv$treatment_effects)       # Time Ratio + HR with 95 % CrI
```

### Body Weight / Toxicity

```r
# Frequentist LME4
bw_results <- analyze_body_weight(
  df               = weight_data,
  time_column      = "Day",
  weight_column    = "Weight",
  id_column        = "ID",
  treatment_column = "Treatment")

# Bayesian LMM
bayes_bw <- bayesian_body_weight(
  df               = weight_data,
  priors           = tg_priors(strength = "skeptical"),
  mcmc             = tg_mcmc())
print(bayes_bw$treatment_effects)
bayes_bw$weight_trajectory_plot
```

### Drug Synergy

```r
# Frequentist Bliss + Loewe
syn_results <- analyze_drug_synergy(
  df               = combo_treatment_synthetic_data,
  drug_a_name      = "DrugA",
  drug_b_name      = "DrugB",
  combo_name       = "Combo",
  control_name     = "Control")
plot_drug_synergy(syn_results)

# Bayesian combination synergy
bayes_syn <- bayesian_synergy(
  df               = combo_treatment_synthetic_data,
  drug_a_name      = "DrugA",
  drug_b_name      = "DrugB",
  combo_name       = "Combo",
  control_name     = "Control",
  priors           = tg_priors(),
  mcmc             = tg_mcmc())
print(bayes_syn$bliss_summary)
print(bayes_syn$loewe_summary)
```

### Dose-Response

```r
data(dose_levels_synthetic_data)
df <- dose_levels_synthetic_data

# Frequentist Hill / Emax (drc backend)
dr_results <- dose_response_statistics(
  df                 = df,
  dose_column        = "Dose",
  volume_column      = "Volume",
  treatment_column   = "Treatment",
  day_column         = "Day",
  id_column          = "ID")

# Bayesian Hill / Emax
# Note: tg_priors() doesn't bundle the Hill-model priors (prior_emax /
# prior_ec50 / prior_hill / prior_sigma) — those are model-specific.
# Use tg_mcmc() for the MCMC config; supply Hill priors individually.
bayes_dr <- bayesian_dose_response(
  df              = df,
  dose_column     = "Dose",
  volume_column   = "Volume",
  treatment_column = "Treatment",
  day_column      = "Day",
  id_column       = "ID",
  mcmc            = tg_mcmc())
print(bayes_dr$posterior_summary)
bayes_dr$pp_check_plot
```

### Therapeutic Window Metric

```r
# Frequentist TWM
twm <- therapeutic_window_metric(tg_results, bw_results)

# Bayesian TWM — single-call convenience wrapper
bayes_twm <- bayesian_twm_from_data(
  df_tg            = tumor_data,
  df_bw            = weight_data,
  priors           = tg_priors(strength = "skeptical"),
  mcmc             = tg_mcmc(chains = 4, iter = 500))
print(bayes_twm$twm_table)
bayes_twm$twm_plot

# Bayesian TWM from separately fitted models
bayes_twm2 <- bayesian_therapeutic_window(
  tg_result        = bayes_tg_result,   # bayesian_tumor_growth(return_model=TRUE)
  bw_result        = bayes_bw_result)   # bayesian_body_weight(return_model=TRUE)
```

### Power Analysis

```r
# Analytic a priori power
power <- apriori_power_analysis(
  effect_sizes = c(0.5, 0.8, 1.2),
  alpha        = 0.05)

# LMM simulation
sim_power <- apriori_power_simulation(
  n_per_group = seq(5, 20, by = 5),
  n_sim       = 500)

# Bayesian simulation
bayes_power <- bayesian_power_analysis(
  n_per_group = seq(5, 20, by = 5),
  target_prob = 0.95,
  n_sim       = 100,
  mcmc        = tg_mcmc(chains = 2, iter = 200))
```

## Functions

### Tumor Growth
| Function | Description |
|----------|-------------|
| `tumor_growth_statistics()` | LME4 / AUC / GAMM tumor growth analysis |
| `bayesian_tumor_growth()` | Bayesian LMM via brms |
| `tumor_auc_analysis()` | Area-under-the-curve analysis |
| `tumor_doubling_time()` | Doubling time estimation |
| `plot_tumor_growth()` | Growth trajectory plot |
| `plot_growth_rate()` | Growth-rate forest plot |

### Survival
| Function | Description |
|----------|-------------|
| `survival_statistics()` | KM, Cox PH, log-rank, `cox.zph()`, C-index |
| `bayesian_survival()` | Bayesian AFT model |

### Body Weight / Toxicity
| Function | Description |
|----------|-------------|
| `analyze_body_weight()` | LME4 mixed model for body weight |
| `bayesian_body_weight()` | Bayesian LMM for body weight |
| `body_weight_auc()` | AUC of body-weight trajectory |
| `therapeutic_window_metric()` | Frequentist TWM |
| `bayesian_therapeutic_window()` | Bayesian TWM from fitted models |
| `bayesian_twm_from_data()` | Single-call Bayesian TWM |
| `efficacy_toxicity_bivariate()` | Bivariate efficacy-toxicity analysis |
| `total_benefit_area()` | Total benefit area metric |
| `weight_corrected_tgi()` | Weight-corrected TGI |

### Drug Synergy
| Function | Description |
|----------|-------------|
| `analyze_drug_synergy()` | Bliss Independence + Loewe CI |
| `analyze_drug_synergy_over_time()` | Time-series synergy analysis |
| `bayesian_synergy()` | Bayesian combination synergy |
| `bayesian_synergy_over_time()` | Bayesian time-series synergy |
| `plot_drug_synergy()`, `plot_synergy_trend()`, `plot_combination_index()`, `plot_bliss()`, `plot_synergy_combined()` | Synergy visualisations |

### Dose-Response
| Function | Description |
|----------|-------------|
| `dose_response_statistics()` | Frequentist Hill / Emax |
| `bayesian_dose_response()` | Bayesian Hill / Emax |

### Power Analysis
| Function | Description |
|----------|-------------|
| `apriori_power_analysis()` | Analytic (t-test / ANOVA) |
| `apriori_power_simulation()` | LMM simulation-based |
| `bayesian_power_analysis()` | Bayesian simulation-based |

### Bayesian configuration helpers (v0.4.7+, CODE_REVIEW.md D.2)
| Function | Description |
|----------|-------------|
| `tg_priors(strength, b, intercept, sd, sigma)` | Prior-config object — bundles the five prior arguments |
| `tg_mcmc(chains, warmup, iter, seed, backend)` | MCMC-config object — bundles the four MCMC arguments |

### Utilities
| Function | Description |
|----------|-------------|
| `calculate_volume()` | Volume from length × width (multiple geometric formulas) |
| `calculate_dates()` | Date-to-day conversion |
| `calculate_auc()` | Trapezoidal AUC |
| `repeated_measures_anova()` | rmANOVA wrapper |
| `export_diagnostics()` | Export model diagnostic plots |
| `new_me_result()` | S3 constructor for typed analysis result lists |

## Datasets

| Dataset | Description |
|---------|-------------|
| `combo_treatment_synthetic_data` | Two-drug combination treatment with tumor + body-weight + survival columns |
| `dose_levels_synthetic_data` | Multiple dose levels of a single drug |
| `master_synthetic_data` | Master synthetic dataset with multiple treatment groups |
| `synthetic_data` | Single-treatment baseline data |
| `combo_treatment_schedule`, `dose_levels_treatment_schedule` | Companion dosing schedules |

## Coverage measurement

```bash
Rscript coverage.R           # per-file coverage summary to stdout
Rscript coverage.R --html    # also write covr_report.html
```

The script excludes Bayesian entry points by default because brms compilation dominates the run. See the K.10 entry in `CODE_REVIEW.md` for the known caveat about stale tests (`test-post_power_analysis.R` and one branch of `test-toxicity_functions.R`) blocking a clean baseline until they're cleaned up under K.11.

## Vignettes

| File | Description |
|------|-------------|
| `vignettes/mouseExperiment.Rmd` | Package overview vignette |
| `vignettes/mouseExperiment_combo_demo.qmd` | Worked combination-treatment example |
| `vignettes/mouseExperiment_dose_demo.qmd` | Worked dose-response example |

Build with `quarto render` (for `.qmd`) or `devtools::build_vignettes()` (for `.Rmd`).

## Development

```r
devtools::load_all()
devtools::test()                # 644 tests; nothing skips a required dependency
devtools::test(filter = "bayesian_tumor_growth")
devtools::document()            # regenerate NAMESPACE + .Rd files
devtools::check()               # R CMD check
Rscript coverage.R              # coverage baseline (excludes Bayesian fits by
                                # default -- a large exclusion now that the whole
                                # Bayesian surface is required and tested)
```

## License

MIT License — see [LICENSE](LICENSE) for details.

## Citation

If you use this package in your research, please cite:

```r
citation("mouseExperiment")
```
