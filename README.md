# mouseExperiment

An R package for statistical analysis of mouse tumor growth experiments. Covers the full preclinical oncology analysis pipeline: tumor growth, survival, body weight toxicity, drug synergy, dose-response, power analysis, and Bayesian equivalents of all major models.

## Features

- **Tumor growth** — LME4 linear mixed-effects model, AUC, Bayesian LMM (brms)
- **Survival** — Kaplan-Meier, Cox PH, log-rank, Bayesian AFT (Weibull/log-normal/exponential/gamma)
- **Body weight / toxicity** — LME4 mixed model, Bayesian LMM, weight-loss threshold events
- **Drug synergy** — Bliss Independence, Loewe Combination Index, Bayesian synergy, synergy over time
- **Dose-response** — Frequentist and Bayesian Hill/Emax curve fitting
- **Therapeutic window** — TWM from frequentist or Bayesian tumor growth and body weight models
- **Power analysis** — Analytic (t-test/ANOVA), LMM simulation, Bayesian simulation
- **Comprehensive plots** — KM curves, growth trajectories, synergy bar charts, dose-response curves, forest plots, MCMC diagnostics

## Installation

```r
# Install the development version from GitHub
devtools::install_github("sciOmics/mouseExperiment")

# For Bayesian analyses also install
install.packages(c("brms", "bayesplot"))
```

## Usage

### Tumor Growth Analysis

```r
library(mouseExperiment)

data(combo_treatment_synthetic_data)

# Frequentist LME4
results <- tumor_growth_statistics(
  data             = combo_treatment_synthetic_data,
  time_column      = "Day",
  volume_column    = "Volume",
  id_column        = "Mouse_ID",
  treatment_column = "Treatment",
  cage_column      = "Cage",
  model_type       = "lme4"
)
print(results$treatment_effects)
plot_tumor_growth(results)

# Bayesian LMM (runs asynchronously in interactive session)
bayes_results <- bayesian_tumor_growth(
  df               = combo_treatment_synthetic_data,
  prior_strength   = "skeptical",
  chains           = 4,
  iter             = 2000
)
print(bayes_results$treatment_effects)
print(bayes_results$mcmc_diagnostics)
bayes_results$pp_check_plot
bayes_results$prior_posterior_plot
```

### Survival Analysis

```r
data(survival_synthetic_data)

# Frequentist
surv_results <- survival_statistics(
  data             = survival_synthetic_data,
  time_column      = "Time",
  event_column     = "Event",
  treatment_column = "Treatment",
  id_column        = "Mouse_ID"
)
print(surv_results$km_summary)

# Bayesian AFT
bayes_surv <- bayesian_survival(
  df               = survival_synthetic_data,
  family           = "weibull",
  prior_strength   = "skeptical"
)
print(bayes_surv$treatment_effects)  # Time Ratio + HR with CrI
```

### Body Weight / Toxicity

```r
# Frequentist LME4
bw_results <- analyze_body_weight(
  data             = weight_data,
  time_column      = "Day",
  weight_column    = "Weight",
  id_column        = "ID",
  treatment_column = "Treatment"
)

# Bayesian LMM
bayes_bw <- bayesian_body_weight(
  df               = weight_data,
  prior_strength   = "skeptical"
)
print(bayes_bw$treatment_effects)
bayes_bw$weight_trajectory_plot
```

### Drug Synergy

```r
# Frequentist Bliss + Loewe
syn_results <- analyze_drug_synergy(
  data             = combo_treatment_synthetic_data,
  time_column      = "Day",
  volume_column    = "Volume",
  id_column        = "Mouse_ID",
  treatment_column = "Treatment"
)
plot_drug_synergy(syn_results)

# Bayesian combination synergy
bayes_syn <- bayesian_synergy(
  df               = combo_treatment_synthetic_data,
  prior_strength   = "skeptical"
)
print(bayes_syn$bliss_summary)
print(bayes_syn$loewe_summary)
```

### Dose-Response

```r
# Frequentist Hill/Emax
dr_results <- dose_response_statistics(
  data             = dose_data,
  dose_column      = "Dose",
  response_column  = "Volume",
  treatment_column = "Treatment"
)

# Bayesian Hill/Emax
bayes_dr <- bayesian_dose_response(
  df               = dose_data,
  prior_strength   = "skeptical"
)
print(bayes_dr$posterior_summary)
bayes_dr$pp_check_plot
bayes_dr$prior_posterior_plot
```

### Therapeutic Window Metric

```r
# Frequentist TWM
twm <- therapeutic_window_metric(tg_results, bw_results)

# Bayesian TWM (single-call convenience wrapper)
bayes_twm <- bayesian_twm_from_data(
  df_tg            = tumor_data,
  df_bw            = weight_data,
  prior_strength   = "skeptical"
)
print(bayes_twm$twm_table)
bayes_twm$twm_plot

# Bayesian TWM from separately fitted models
bayes_twm2 <- bayesian_therapeutic_window(
  tg_result        = bayes_tg_result,  # from bayesian_tumor_growth(return_model = TRUE)
  bw_result        = bayes_bw_result   # from bayesian_body_weight(return_model = TRUE)
)
```

### Power Analysis

```r
# Analytic a priori power
power <- apriori_power_analysis(
  effect_sizes = c(0.5, 0.8, 1.2),
  alpha        = 0.05
)

# LMM simulation
sim_power <- apriori_power_simulation(
  n_per_group = seq(5, 20, by = 5),
  n_sim       = 500
)

# Bayesian simulation
bayes_power <- bayesian_power_analysis(
  n_per_group    = seq(5, 20, by = 5),
  target_prob    = 0.95,
  n_sim          = 100
)
```

## Functions

### Tumor Growth
| Function | Description |
|----------|-------------|
| `tumor_growth_statistics()` | LME4 or AUC tumor growth analysis |
| `bayesian_tumor_growth()` | Bayesian LMM via brms |
| `tumor_auc_analysis()` | Area-under-the-curve analysis |
| `tumor_doubling_time()` | Doubling time estimation |
| `plot_tumor_growth()` | Growth trajectory plot |
| `plot_growth_rate()` | Growth rate forest plot |

### Survival
| Function | Description |
|----------|-------------|
| `survival_statistics()` | KM, Cox PH, log-rank |
| `bayesian_survival()` | Bayesian AFT model |

### Body Weight / Toxicity
| Function | Description |
|----------|-------------|
| `analyze_body_weight()` | LME4 mixed model for body weight |
| `bayesian_body_weight()` | Bayesian LMM for body weight |
| `body_weight_auc()` | AUC of body weight trajectory |
| `weight_loss_threshold()` | Time-to-threshold event derivation |
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
| `bayesian_synergy()` | Bayesian combination synergy |
| `analyze_drug_synergy_over_time()` | Time-series synergy analysis |
| `bayesian_synergy_over_time()` | Bayesian time-series synergy |
| `plot_drug_synergy()` | Synergy bar chart |
| `plot_synergy_trend()` | Synergy over time plot |

### Dose-Response
| Function | Description |
|----------|-------------|
| `dose_response_statistics()` | Frequentist Hill/Emax |
| `bayesian_dose_response()` | Bayesian Hill/Emax |

### Power Analysis
| Function | Description |
|----------|-------------|
| `apriori_power_analysis()` | Analytic power (t-test / ANOVA) |
| `apriori_power_simulation()` | LMM simulation-based power |
| `bayesian_power_analysis()` | Bayesian simulation-based power |

### Utilities
| Function | Description |
|----------|-------------|
| `calculate_volume()` | Volume from length × width |
| `calculate_dates()` | Date-to-day conversion |
| `calculate_auc()` | Trapezoidal AUC |
| `repeated_measures_anova()` | rmANOVA wrapper |
| `export_diagnostics()` | Export model diagnostic plots |

## Datasets

| Dataset | Description |
|---------|-------------|
| `combo_treatment_synthetic_data` | Combination treatment tumor growth data |
| `survival_synthetic_data` | Survival data with multiple treatment groups |
| `single_treatment_synthetic_data` | Single-treatment tumor growth data |

## License

MIT License — see [LICENSE](LICENSE) for details.

## Citation

If you use this package in your research, please cite:

```r
citation("mouseExperiment")
```
