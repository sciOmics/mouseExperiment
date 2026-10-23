# mouseExperiment

An R package for statistical analysis of mouse tumor growth experiments. Covers the preclinical-oncology analysis pipeline: tumor growth, survival, body weight toxicity, drug synergy, dose-response and power analysis, with Bayesian models for tumor growth and survival. Every analysis can also be run without code, from an interactive dashboard (`run_dashboard()`).

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
- **Dashboard** — `run_dashboard()` opens an interactive Shiny app: upload a CSV, map its columns and run each analysis from the browser, with downloadable tables, plots and an HTML report. The hosted version is at [mouseexperiment.com](https://mouseexperiment.com)
- **Code-coverage** — `covr::package_coverage()` baseline; one-line `Rscript coverage.R` for HTML / stdout reports

## Current status

| Item | State |
|---|---|
| Version | see `DESCRIPTION` |
| Review | `CODE_REVIEW.md`, cumulative rounds; Round 20 is the latest |
| Checks | `R CMD check` and the test suite run on every push and nightly (GitHub Actions) |
| README | every non-Bayesian code block below is run by `tests/testthat/test-readme.R` |

### A note on dependencies

The statistical packages (`brms`, `bayesplot`, `coin`, `clinfun`, `ggpubr`,
`posterior`) are **required**, not suggested, so installing needs a working C++
toolchain: `brms` pulls in the Stan stack.

That cost is deliberate. While `brms` sat in `Suggests`, the Bayesian tests were
skipped whenever it was absent, and two Critical defects survived five releases
behind that skip (`CODE_REVIEW.md` R3-L and R4).

`cmdstanr` remains optional. It is an alternative Stan backend, and asking for it
when it is missing fails with install instructions.

## Installation

```r
# From GitHub (brms and the Stan stack are installed as dependencies)
remotes::install_github("sciOmics/mouseExperiment")

# Optional: the cmdstanr backend, 3-10x faster Stan compilation
install.packages("cmdstanr",
                 repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
cmdstanr::install_cmdstan()
```

## Dashboard

```r
mouseExperiment::run_dashboard()
```

The dashboard opens in your browser and runs on your computer, so your data never
leaves it. Its packages (DT, shinyjs, future, promises and a few others) are
Suggests, so installing `mouseExperiment` for its functions does not pull them in;
`run_dashboard()` offers to install any that are missing. Usage analytics are off
unless a server enables them.

The developer documentation is in [`docs/dashboard/`](docs/dashboard/README.md).
Until v0.34.0 the dashboard was a separate package, `mouseExperimentDashboard`.

## Usage

Every example uses `master_synthetic_data`, which ships with the package and has
every column the functions use. The vignette, `vignette("mouseExperiment")`, walks
through the same analysis with explanations.

### Preparing data

```r
library(mouseExperiment)

# Raw laboratory data: calendar dates and caliper measurements
d <- calculate_dates(combo_treatment_synthetic_data, start_date = "03/24/2025")
d <- calculate_volume(d)     # ellipsoid, pi * L * W^2 / 6
head(d)

m <- master_synthetic_data  # already has Day and Volume
evaluable_days(m, volume_column = "Volume")$last_day
```

### Tumor growth

```r
library(mouseExperiment)
m <- master_synthetic_data

tg <- tumor_growth_statistics(m, reference_group = "Vehicle",
                              p_adjust_method = "holm")
tg$anova                   # Type III F-tests (Satterthwaite)
tg$pairwise_comparisons    # each arm vs Vehicle, Holm-adjusted
tg$random_effects          # structure requested and used
plot_tumor_growth(m)
plot_growth_rate(tg$growth_rates)
```

```r
# Bayesian LMM (compiles a Stan model on first use, about 30 s)
library(mouseExperiment)
m <- master_synthetic_data
fit <- bayesian_tumor_growth(
  subset(m, Treatment %in% c("Vehicle", "Drug_A High")),
  cage_column = "Cage",
  priors = tg_priors(strength = "weakly_informative"),
  mcmc   = tg_mcmc(chains = 4, warmup = 1000, iter = 500, seed = 42))
fit$treatment_effects      # group means with 95 % credible intervals
fit$mcmc_diagnostics       # Rhat, ESS, convergence flag
fit$loo_diagnostics        # elpd_loo, Pareto-k
```

Every analysis function takes `plots = FALSE` to skip its figures. The model diagnostics (residual and Q-Q plots, influence) come with `include_diagnostics = TRUE` where a function has it.

### Survival

```r
library(mouseExperiment)
m <- master_synthetic_data

# One row per animal: its last record carries the event indicator
mm <- m[order(m$ID, m$Day), ]
per_animal <- mm[!duplicated(mm$ID, fromLast = TRUE), ]

surv <- survival_statistics(per_animal, reference_group = "Vehicle",
                            verbose = FALSE)
surv$results               # hazard ratios, intervals, p-values, medians
surv$method_used           # "coxphf" (Firth) by default
```

```r
# Bayesian AFT model (Weibull or log-normal)
library(mouseExperiment)
m <- master_synthetic_data
mm <- m[order(m$ID, m$Day), ]
per_animal <- mm[!duplicated(mm$ID, fromLast = TRUE), ]
bayes_surv <- bayesian_survival(per_animal, reference_group = "Vehicle",
                                family = "weibull",
                                mcmc = tg_mcmc(chains = 4, iter = 500))
bayes_surv$treatment_effects   # time ratio and hazard ratio, with 95 % CrI
```

### Body weight and toxicity

```r
library(mouseExperiment)
m <- master_synthetic_data

bw <- analyze_body_weight(m, volume_column = "Volume", volume_units = "mm3",
                          cage_column = "Cage", reference_group = "Vehicle")
bw$pairwise_comparisons

wl <- weight_loss_threshold(m, cage_column = "Cage", threshold = 0.1,
                            reference_group = "Vehicle")
```

### Drug synergy

```r
library(mouseExperiment)
m <- master_synthetic_data

syn <- analyze_drug_synergy(m,
  drug_a_name = "Drug_A Mid", drug_b_name = "Drug_B",
  combo_name = "Drug_A Mid + Drug_B", control_name = "Vehicle",
  cage_column = "Cage", n_boot = 500, boot_seed = 1, verbose = FALSE)
syn$summary
plot_drug_synergy(syn)
```

### Dose-response

```r
library(mouseExperiment)
m <- master_synthetic_data

dr <- dose_response_statistics(m,
  treatments = c("Drug_A Low", "Drug_A Mid", "Drug_A High"),
  control_group_name = "Vehicle", cage_column = "Cage", verbose = FALSE)
dr$summary_table
dr$trend_test        # Jonckheere-Terpstra
```

### Therapeutic window

```r
library(mouseExperiment)
m <- master_synthetic_data

tw <- therapeutic_window_metric(m,
  reference_group = "Vehicle", cage_column = "Cage",
  volume_units = "mm3",               # needed to subtract tumour mass
  tolerability_threshold = 20,        # percent of baseline weight
  n_boot = 500, boot_seed = 1)
tw$window_table   # TGI and worst weight loss per arm, with intervals and a flag
```

### Power analysis

```r
library(mouseExperiment)

# Analytic: animals per group for a standardised effect of 1.2, three groups
apriori_power_analysis(effect_size = 1.2, n_groups = 3)

# Simulation with the same mixed model (use 500+ simulations in practice)
apriori_power_simulation(n_per_group = c(6L, 10L), n_simulations = 50L, seed = 1)
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
| `clear_compiled_model_cache()` | Empties the session's cache of compiled Stan models; a fit whose model has the same structure as an earlier one reuses its compiled model and skips the ~20 s compile (v0.28.0) |

### Utilities
| Function | Description |
|----------|-------------|
| `calculate_volume()` | Volume from length × width (multiple geometric formulas) |
| `calculate_dates()` | Date-to-day conversion |
| `calculate_auc()` | Trapezoidal AUC |

## Datasets

| Dataset | Description |
|---------|-------------|
| `master_synthetic_data` | Six arms with every column the package uses: Day, Length/Width/Height, Volume, Weight, Dose, Cage, Necrotic and Survival_Censor (the dashboard's "All Features" demo) |
| `combo_treatment_synthetic_data` | Two-drug combination as a laboratory records it: dates, ear-punch IDs, caliper length and width, survival events |
| `dose_levels_synthetic_data` | One drug at three doses and a vehicle control: dates, caliper length and width, survival events |
| `weight_synthetic_data` | Tumour volume with body weight; column names contain spaces |

## Coverage measurement

```bash
Rscript coverage.R           # per-file coverage summary to stdout
Rscript coverage.R --html    # also write covr_report.html
```

The script leaves out the Bayesian entry points by default, because brms
compilation dominates the run.

## Vignette

`vignette("mouseExperiment")`, a complete analysis of `master_synthetic_data`.

## Development

```r
devtools::load_all()
devtools::test()                # the full suite; nothing skips a required dependency
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
