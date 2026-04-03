# mouseExperimentDashboard

An interactive web-based dashboard for analyzing mouse tumor growth experiments. Provides a user-friendly, point-and-click interface to the `mouseExperiment` R package.

**Note**: This is a proprietary application for private use. See LICENSE for details.

## Overview

`mouseExperimentDashboard` wraps the statistical and visualization capabilities of `mouseExperiment`, exposing them through a Shiny web application. Users can:

- Upload tumor growth data (CSV or RDA format)
- Perform statistical analyses using linear mixed-effects models (LME4) or area-under-the-curve (AUC) methods
- Analyze survival time-to-event data
- Generate customizable plots
- Download results in multiple formats


## Local Installation & Development

### Prerequisites

- **R** ≥ 4.0
- **mouseExperiment** package (see [mouseExperiment repository](https://github.com/sciOmics/mouseExperiment))

### Setup for Development

```r
# Clone the repository (or navigate to local working directory)
cd /path/to/mouseExperimentDashboard

# Load dependencies
devtools::load_all()

# Run the app locally
runDashboard()
```

The app will open at `http://localhost:3838` (or another available port).

## Usage

### Basic Workflow — Select a CSV or RDA file containing tumor growth data in long format (one row per observation). 
2. **Column Mapping** — Map the ID, time, volume, treatment, and optional columns to columns in the data
3. **Tumor Growth Analysis Setup** — Choose statistical model (LME4 or AUC), data transformation, and additional parameters
4. **Run Analysis** — Click "Run Analysis" to fit the model
5. **Explore Results** — View summary statistics, ANOVA, pairwise comparisons, growth rates, and plots
6. **Download** — Customize and export plots and tables in PNG, PDF, or CSV format

## Data Format

The dashboard expects data in **long format**: one row per observation/measurement.

### Required Columns

- **ID / Mouse ID / Animal ID / Subject** — Unique identifier for each animal
- **Day / Time / Timepoint** — Study day or date 
- **Tumor Measurement(s)** — Tumor volume measurement or tumor length and width measurements
- **Treatment / Group** — Treatment group assignment 

### Optional Columns

- **Cage / Cage_ID** — Cage identifier (for random effects in LME4)
- **Survival_Censor / Censor / Event** — Event indicator (numeric: 1 = event, 0 = censored) for survival analysis
- **Dose / Dose Level** — Dose level (for dose-response analysis)
- **Weight / Mass** — Animal weight/mass (for toxicity analysis)

### Example Data Structure

| ID | Day | Volume | Treatment | Cage |  
|-------|--------|------|------|------|    
Mouse_001|0|45|Control|A1|  
Mouse_001|3|52|Control|A1|
Mouse_001|7|63|Control|A1|
Mouse_002|0|48|Control|A1|
Mouse_002|3|55|Control|A1|
Mouse_002|7|68|Control|A1|
Mouse_003|0|41|Drug_A|A2|
Mouse_003|3|42|Drug_A|A2|
Mouse_003|7|45|Drug_A|A2|


### Tumor Volume Calculation from Length & Width

The dashboard can calculate tumor volume from length and width measurements using several formulas.

#### Available Methods

The following volume formulas are supported (all assume perpendicular measurements):

| Formula | Equation | Notes |
|---------|----------|--------|
| **Ellipsoid** (default) | V = (L × W² × π) / 6 | Standard for 2D measurements | 
| **Modified Ellipsoid** | V = (L × W²) / 2 | Simplified formula used in older literature |
| **Cylinder** | V = (π × W² × L) / 4 | Used for more elongated/cylindrical tumors, narrow or needle-like |
| **Sphere** | V = (π × W³) / 6 | Spherical tumors; uses only diameter (width) |

#### Automatic Detection & Correction

The dashboard automatically:

1. **Detects** Length and Width columns by name pattern matching (e.g., "Length", "length", "L", "Width", "width", "W")
2. **Corrects dimensions**: Ensures the longer measurement is used as **Length** and the shorter as **Width**, regardless of how they were labeled in your data
   - If your data has Length=5 and Width=10, the system internally swaps them to L=10, W=5
   - This prevents incorrect volume calculations due to mis-labeled dimensions
   - Rows with swapped values are highlighted
3. **Selects an appropriate formula** based on available columns (defaults to Ellipsoid if 3D height unavailable)

#### Example

```
# Input data (width column larger than length column)
ID,Day,Length,Width,Treatment
M001,0,8,12,Control

# Dashboard correction
# → Internally uses: L=12, W=8

# Ellipsoid calculation:
# V = (12 × 8² × π) / 6 ≈ 201.1 mm³
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

## Planned Features for Future Releases
- **Bayesian Statistial Models**
- **Toxicity Analysis** — Assessment of toxicity using mouse weight
- **Handling of Necrotic Tumors** — Model tumor necrois in statistical analysis
- **Dose-Response Analysis** — Relationship between dose level and tumor volume
- **Drug Synergy Analysis** — Bliss independence and Combination Index (CI) evaluation; disabled pending UI/server completion
- **Power Analysis** — Sample size and power calculations for future experiments

# Citation

If you use the underlying `mouseExperiment` R package in your research, please cite:

```
mouseExperiment: An R package for statistical analysis of mouse tumor growth experiments.
https://github.com/sciOmics/mouseExperiment
```

**Note**: `mouseExperimentDashboard` is a proprietary private application and is not for citation or redistribution.

## License

**Proprietary** — All Rights Reserved. See [LICENSE](LICENSE) for details.

mouseExperiment (the underlying R package) is licensed under the MIT License.

## Support

For questions or issues, contact: david@sciomics.bio

