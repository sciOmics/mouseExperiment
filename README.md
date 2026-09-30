# mouseExperiment

An R package for statistical analysis of mouse tumor growth experiments. Covers the preclinical-oncology analysis pipeline: tumor growth, survival, body weight toxicity, drug synergy, dose-response and power analysis, with Bayesian models for tumor growth and survival.

## Features

- **Tumor growth** — LME4 LMM, AUC, Bayesian LMM (brms)
- **Survival** — Kaplan-Meier, Cox PH, log-rank, `cox.zph()` PH check, Bayesian AFT (Weibull / log-normal), C-index
- **Body weight / toxicity** — LME4 mixed model, weight-loss threshold events
- **Drug synergy** — Bliss Independence with bootstrap intervals, plus an over-time variant
- **Dose-response** — Hill / Emax curve fitting (drc) for one agent's dose series, EC50 with a log-dose interval, two-sided Jonckheere–Terpstra trend test
- **Therapeutic window** — TGI against worst weight loss per arm, each with an interval, plus a tolerability flag at a declared threshold
- **Power analysis** — Analytic (t-test / ANOVA) and LMM simulation; multiplicity- and attrition-aware (`n_comparisons`, `dropout_rate`)
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

As of v0.10.0 the statistical packages (`brms`, `bayesplot`, `coin`,
`clinfun`, `ggpubr`, `posterior`) are **required**, not suggested.
Installing therefore needs a working C++ toolchain, because `brms` pulls the Stan
stack.

That cost is deliberate. While `brms` sat in `Suggests`, the Bayesian tests skipped
wholesale whenever it was absent — and two Critical defects survived five releases
behind that skip, including `bayesian_synergy()` (since removed) being entirely
non-functional from v0.4.6 to v0.9.0 while this README advertised it. Optionality
was not free; it was the mechanism. See `CODE_REVIEW.md` §R3-L and §R4.

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
```

### Drug Synergy

```r
# Frequentist Bliss independence
syn_results <- analyze_drug_synergy(
  df               = combo_treatment_synthetic_data,
  drug_a_name      = "DrugA",
  drug_b_name      = "DrugB",
  combo_name       = "Combo",
  control_name     = "Control")
plot_drug_synergy(syn_results)
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
```

### Therapeutic Window

```r
data(master_synthetic_data)
tw <- therapeutic_window_metric(
  df               = master_synthetic_data,
  reference_group  = "Vehicle",
  cage_column      = "Cage",
  volume_units     = "mm3",           # needed to subtract tumour mass
  tolerability_threshold = 20)        # percent of baseline weight
tw$window_table   # TGI and worst weight loss per arm, with intervals and a flag
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
```

## Functions

### Tumor Growth
| Function | Description |
|----------|-------------|
| `tumor_growth_statistics()` | LME4 / AUC tumor growth analysis |
| `bayesian_tumor_growth()` | Bayesian LMM via brms |
| `tumor_doubling_time()` | Doubling time estimation |
| `plot_tumor_growth()` | Growth trajectory plot |
| `plot_growth_rate()` | Growth-rate forest plot |

### Survival
| Function | Description |
|----------|-------------|
| `survival_statistics()` | KM, Cox PH, log-rank, `cox.zph()`, C-index |
| `bayesian_survival()` | Bayesian AFT model (Weibull / log-normal) |

### Body Weight / Toxicity
| Function | Description |
|----------|-------------|
| `analyze_body_weight()` | LME4 mixed model for body weight |
| `weight_loss_threshold()` | Time to a weight-loss threshold (KM / competing risks) |
| `therapeutic_window_metric()` | Therapeutic window: efficacy vs weight loss |

### Drug Synergy
| Function | Description |
|----------|-------------|
| `analyze_drug_synergy()` | Bliss Independence, with bootstrap intervals |
| `analyze_drug_synergy_over_time()` | Time-series synergy analysis |
| `plot_drug_synergy()`, `plot_synergy_trend()`, `plot_bliss()` | Synergy visualisations |

### Dose-Response
| Function | Description |
|----------|-------------|
| `dose_response_statistics()` | Hill / Emax dose-response |

### Power Analysis
| Function | Description |
|----------|-------------|
| `apriori_power_analysis()` | Analytic (t-test / ANOVA) |
| `apriori_power_simulation()` | LMM simulation-based |

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

## Datasets

| Dataset | Description |
|---------|-------------|
| `combo_treatment_synthetic_data` | Two-drug combination treatment with tumor + body-weight + survival columns |
| `dose_levels_synthetic_data` | Multiple dose levels of a single drug |
| `master_synthetic_data` | Master synthetic dataset with multiple treatment groups |
| `synthetic_data` | Single-treatment baseline data |

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
