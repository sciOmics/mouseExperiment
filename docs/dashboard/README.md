# mouseExperimentDashboard

An interactive Shiny web application for analyzing mouse tumor growth experiments. Provides a point-and-click interface to the `mouseExperiment` R package, covering the full analysis pipeline from raw data upload through publication-ready figures.

**Note**: This is a proprietary application for private use. See LICENSE for details.

## Overview

`mouseExperimentDashboard` wraps the statistical and visualization capabilities of `mouseExperiment`, exposing them through a Shiny web application. Users can:

- Upload tumor growth data (CSV or RDA format) and map columns interactively
- Analyze tumor growth with LME4 mixed-effects models, AUC, or Bayesian LMM (brms)
- Run survival analysis with Kaplan-Meier curves, Cox/AFT models, and Bayesian AFT (brms)
- Analyze drug synergy via Bliss Independence and Loewe CI, with a Bayesian combination synergy path
- Fit dose-response curves (frequentist Hill/Emax or Bayesian Hill/Emax)
- Assess toxicity via body weight mixed-effects or Bayesian LMM, with an Efficacy-Toxicity mode including Therapeutic Window Metric and Bayesian TWM
- Compute a priori power using analytic formulas, LMM simulation, or Bayesian simulation
- Customize and download all plots and tables

All long-running Bayesian fits run asynchronously via `future`/`promises` so the app remains responsive during Stan compilation and MCMC sampling.

## Local Installation & Development

### Prerequisites

- **R** ≥ 4.1
- **mouseExperiment** ≥ 0.4.3 (see [mouseExperiment repository](https://github.com/sciOmics/mouseExperiment))
- For Bayesian analysis: `brms` ≥ 2.19, `future` ≥ 1.33, `promises` ≥ 1.2, `bayesplot` ≥ 1.10

### Setup for Development

```r
# Navigate to local working directory
setwd("/path/to/mouseExperimentDashboard")

# Load dependencies
devtools::load_all()

# Run the app locally
runDashboard()
```

The app will open at `http://localhost:3838` (or another available port).

## Usage

### Basic Workflow

1. **Data Upload** — Select a CSV or RDA file containing data in long format (one row per observation)
2. **Column Mapping** — Map ID, time, volume, treatment, and optional columns
3. **Analysis Setup** — Choose analysis module, statistical method, and sidebar parameters
4. **Run Analysis** — Click "Run Analysis" to fit the model; Bayesian fits run in the background
5. **Explore Results** — Navigate result tabs: summary statistics, model output, plots, pairwise comparisons
6. **Download** — Export plots (PNG) and tables (CSV) from per-tab download buttons

### Analysis Modules

| Module | Methods |
|--------|---------|
| **Tumor Growth** | LME4 mixed model, AUC, Bayesian LMM (brms), Dose-Response (frequentist + Bayesian) |
| **Survival** | Kaplan-Meier, Cox PH, log-rank, Bayesian AFT (Weibull / log-normal / exponential / gamma) |
| **Toxicity** | Body weight LME4 or Bayesian LMM; Efficacy-Toxicity bivariate; Therapeutic Window Metric (frequentist + Bayesian) |
| **Drug Synergy** | Bliss Independence, Loewe CI (frequentist); Bayesian combination synergy |
| **Power Analysis** | A priori analytic (t-test / ANOVA), LMM simulation, Bayesian simulation |

## Data Format

The dashboard expects data in **long format**: one row per observation/measurement.

### Required Columns

- **ID / Mouse ID / Animal ID** — Unique identifier for each animal
- **Day / Time / Timepoint** — Study day or date
- **Tumor Measurement(s)** — Tumor volume, or tumor length and width
- **Treatment / Group** — Treatment group assignment

### Optional Columns

- **Cage / Cage_ID** — Cage identifier (random intercept in LME4/Bayesian models)
- **Survival_Censor / Event** — Event indicator (1 = event, 0 = censored) for survival analysis
- **Dose** — Dose level (for dose-response analysis)
- **Weight / Mass** — Animal weight/mass (for toxicity analysis)

### Example Data Structure

| ID | Day | Volume | Treatment | Cage |
|----|-----|--------|-----------|------|
| Mouse_001 | 0 | 45 | Control | A1 |
| Mouse_001 | 7 | 63 | Control | A1 |
| Mouse_002 | 0 | 41 | Drug_A | A2 |
| Mouse_002 | 7 | 45 | Drug_A | A2 |

### Tumor Volume Calculation from Length & Width

When volume is not provided directly, the dashboard calculates it from length and width:

| Formula | Equation |
|---------|----------|
| Ellipsoid (default) | V = (L × W² × π) / 6 |
| Modified Ellipsoid | V = (L × W²) / 2 |
| Cylinder | V = (π × W² × L) / 4 |
| Sphere | V = (π × W³) / 6 |

Dimensions are auto-corrected so that the longer measurement is always used as length.

## Requirements

| Package | Minimum Version | Purpose |
|---------|-----------------|---------|
| R | 4.1 | Language runtime |
| shiny | 1.7.0 | Web framework |
| mouseExperiment | 0.4.3 | Statistical backend |
| ggplot2 | 3.4.0 | Visualization |
| emmeans | 1.8.0 | Estimated marginal means |
| lme4 | 1.1-31 | Mixed-effects models |
| survival | 3.5 | Survival analysis |
| DT | 0.28 | Interactive tables |
| dplyr | 1.1.0 | Data manipulation |
| shinyjs | — | UI helpers, download button state |
| shinycssloaders | — | Spinner overlays |
| brms | 2.19 | Bayesian models (optional) |
| future | 1.33 | Async execution (optional) |
| promises | 1.2 | Async piping (optional) |
| bayesplot | 1.10 | MCMC diagnostic plots (optional) |

See `DESCRIPTION` for the complete dependency list.

## Key Files

| File | Purpose |
|------|---------|
| `R/app.R` | App entry point; UI layout, server dispatch, async plan setup |
| `R/mod_data_upload.R` | Data upload, column mapping, preprocessing cache |
| `R/mod_tumor_growth.R` | Tumor growth analysis UI + server (LME4, AUC, Bayesian LMM, Dose-Response) |
| `R/survival_modules.R` | Survival analysis UI + server (KM, Cox, Bayesian AFT) |
| `R/toxicity_modules.R` | Toxicity UI + server (body weight LMM, ET bivariate, TWM, Bayesian variants) |
| `R/synergy_modules.R` | Drug synergy UI + server (Bliss, Loewe, Bayesian synergy) |
| `R/dose_response_modules.R` | Dose-response UI + server (frequentist + Bayesian Hill/Emax) |
| `R/power_modules.R` | Power analysis UI + server (analytic, LMM sim, Bayesian sim) |
| `R/helpers.R` | Shared utilities: `translate_error()`, `toggle_downloads()`, `js_string()`, color helpers |
| `R/ui_modules.R` | Shared UI components and layout helpers |
| `www/css/` | Custom styles (light and dark modes) |
| `inst/sample_data/` | Example datasets |

## Development

```r
# Load dependencies
devtools::load_all()

# Run test suite
devtools::test()

# Check the package
devtools::check()
```

## License

**Proprietary** — All Rights Reserved. See [LICENSE](LICENSE) for details.

`mouseExperiment` (the underlying R package) is licensed under the MIT License.

## Citation

If you use the underlying `mouseExperiment` R package in your research, please cite:

```
mouseExperiment: An R package for statistical analysis of mouse tumor growth experiments.
https://github.com/sciOmics/mouseExperiment
```

`mouseExperimentDashboard` is a proprietary private application and is not for citation or redistribution.

## Support

For questions or issues, contact: david@sciomics.bio
