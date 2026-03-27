# mouseExperimentDashboard

An interactive web-based dashboard for analyzing mouse tumor growth experiments. Provides a user-friendly, point-and-click interface to the `mouseExperiment` R package — no programming required.

## Overview

`mouseExperimentDashboard` wraps the statistical and visualization capabilities of `mouseExperiment`, exposing them through a Shiny web application. Users can:

- Upload tumor growth data (CSV or RDA format)
- Perform statistical analyses using linear mixed-effects models (LME4) or area-under-the-curve (AUC) methods
- Generate customizable publication-quality plots
- Download results in multiple formats
- Analyze survival time-to-event data
- Evaluate drug synergy using Bliss independence and Loewe additivity models

## Installation

### Prerequisites

- **R** ≥ 4.0
- **mouseExperiment** package (see [mouseExperiment repository](https://github.com/sciOmics/mouseExperiment))

### Install the Dashboard Package

```r
# From GitHub (requires devtools)
devtools::install_github("sciOmics/mouseExperimentDashboard")
```

Or in development mode:

```r
devtools::load_all()
```

## Usage

Start the dashboard with:

```r
library(mouseExperimentDashboard)
runDashboard()
```

The app will open at `http://localhost:3838` (or another available port).

### Basic Workflow

1. **Data Upload** — Select a CSV or RDA file containing tumor growth data in long format (one row per observation)
2. **Column Mapping** — Ensure the dashboard identifies your ID, time, volume, treatment, and optional columns correctly
3. **Analysis Setup** — Choose statistical model (LME4 or AUC), data transformation, and polynomial degree
4. **Run Analysis** — Click "Run Analysis" to fit the model
5. **Explore Results** — View summary statistics, ANOVA, pairwise comparisons, growth rates, and plots
6. **Download** — Export plots and tables in PNG, PDF, or CSV format

## Features

### Implemented

- ✅ **Data Upload & Mapping**: CSV/RDA file import with automatic column detection
- ✅ **Tumor Growth Analysis**
  - LME4 (linear mixed-effects) models with random cage/animal effects
  - AUC (area under the curve) method for per-animal summary statistics
  - Data transformations: linear, log, or square-root scale
  - Polynomial terms for non-linear growth trajectories
- ✅ **Statistical Outputs**
  - Type III ANOVA (LME4) or one-way ANOVA (AUC)
  - Pairwise contrasts via emmeans (Dunnett's test or Tukey HSD)
  - Estimated marginal means (treatment effects)
  - Per-animal growth rates with confidence intervals
- ✅ **Visualization**
  - Tumor growth curves with optional group average overlay
  - Growth rates plot (spread of per-animal rates by group)
  - Waterfall plot (per-animal growth rates sorted by magnitude)
  - Effect size forest plot (pairwise comparisons with 95% CIs)
- ✅ **Survival Analysis**
  - Kaplan-Meier survival curves
  - Cox proportional hazards model with Hazard Ratios
  - Numbers-at-risk tables
- ✅ **Plot Customization**
  - Collapsible control panels (Tumor Growth, Survival, Synergy)
  - Group ordering and color scheme selection
  - Legend positioning (auto, none, top, bottom, left, right, in-plot)
  - Text size, title, subtitle, and export resolution (DPI)
  - Width/height and download buttons (PNG, PDF)
- ✅ **Model Diagnostics**
  - Q-Q plots of residuals for normality checking
  - Residuals vs. fitted values for heteroscedasticity assessment

### Planned (Disabled)

- 🔲 **Toxicity Analysis** — Assessment of organ toxicity endpoints; placeholder tab
- 🔲 **Dose-Response Analysis** — Relationship between dose level and tumor volume
- 🔲 **Drug Synergy Analysis** — Bliss independence and Combination Index (CI) evaluation; disabled pending UI/server completion
- 🔲 **Power Analysis** — Sample size and power calculations for future experiments

## Data Format

The dashboard expects data in **long format**: one row per observation/measurement.

### Required Columns

- **ID / Mouse ID / Animal ID / Subject** — Unique identifier for each animal
- **Day / Time / Timepoint** — Study day or timepoint (numeric)
- **Volume / Tumor_Volume** — Tumor volume measurement (numeric, mm³ or other units)
- **Treatment / Group** — Treatment group assignment (character/factor)

### Optional Columns

- **Cage / Cage_ID** — Cage identifier (for random effects in LME4)
- **Dose / DoseLevel** — Dose level (for dose-response analysis)
- **Survival_Censor / Censor / Event** — Event indicator (numeric: 1 = event, 0 = censored) for survival analysis
- **Weight / Mass** — Animal weight/mass (for toxicity analysis)

### Example Structure

```
ID,Day,Volume,Treatment,Cage
Mouse_001,0,45,Control,A1
Mouse_001,3,52,Control,A1
Mouse_001,7,63,Control,A1
Mouse_002,0,48,Control,A1
Mouse_002,3,55,Control,A1
Mouse_002,7,68,Control,A1
Mouse_003,0,41,Drug_A,A2
Mouse_003,3,42,Drug_A,A2
Mouse_003,7,45,Drug_A,A2
...
```

## Requirements

| Package | Minimum Version | Purpose |
|---------|-----------------|---------|
| R | 4.0 | Language runtime |
| shiny | 1.7.0 | Web framework |
| mouseExperiment | 0.3.0 | Statistical backend |
| ggplot2 | 3.4.0 | Visualization |
| emmeans | 1.8.0 | Estimated marginal means |
| lme4 | — | Mixed-effects models (via mouseExperiment) |
| DT | 0.28 | Interactive tables |
| dplyr | 1.1.0 | Data manipulation |
| tidyr | 1.3.0 | Data reshaping |

See `DESCRIPTION` for the complete dependency list.

## Development

### Setup

```r
# Clone the repository
git clone https://github.com/sciOmics/mouseExperimentDashboard.git
cd mouseExperimentDashboard

# Load dependencies
devtools::load_all()

# Run the app
runDashboard()
```

### Testing

```r
# Run test suite
devtools::test()

# Check the package
devtools::check()
```

### Key Files

- `R/app.R` — Main app UI and server layout
- `R/ui_modules.R` — Shared UI helper functions
- `R/server_modules.R` — Tumor Growth Analysis server logic
- `R/survival_modules.R` — Survival Analysis UI + server
- `R/synergy_modules.R` — Drug Synergy Analysis UI + server
- `www/css/` — Custom styles (light and dark modes)
- `inst/sample_data/` — Example datasets

## Troubleshooting

### File won't upload
- Ensure the file is a valid CSV or RDA
- Check file permissions

### Analysis tab is disabled/greyed out
- Upload a data file first
- Ensure all required columns are mapped in the Data Upload tab
- For Survival analysis, ensure you have a censor/event column
- For specific analyses (Dose Response, Synergy), ensure you have the required columns (Dose, multiple treatment groups)

### Columns not detected
- The dashboard looks for common column names. If your columns have non-standard names, use the Data Upload tab's column selector to map them manually
- Expected patterns: ID, Mouse_ID, Animal_ID, Subject; Day, Time, Timepoint; Volume, Tumor_Volume; Treatment, Group

### Analysis fails or shows errors
- Check that the data contains valid numeric values for volume and time
- Ensure you have at least 3 timepoints per animal for growth rate estimation
- For LME4 models, check that cage and animal levels are not perfectly confounded

## Citation

If you use `mouseExperimentDashboard` or `mouseExperiment` in your research, please cite:

```
mouseExperiment: An R package for statistical analysis of mouse tumor growth experiments.
https://github.com/sciOmics/mouseExperiment
```

## License

MIT License — See [LICENSE](LICENSE) for details.

## Support

For issues, feature requests, or questions:
- 📤 [Report an issue on GitHub](https://github.com/sciOmics/mouseExperiment/issues)
- 📖 [View the mouseExperiment documentation](https://github.com/sciOmics/mouseExperiment)
