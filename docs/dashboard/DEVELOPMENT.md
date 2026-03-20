# Development Guide for mouseExperimentDashboard

## Quick Start

### 1. Prerequisites

First, ensure mouseExperiment is loaded from the local development folder:

```r
# Load mouseExperiment from local path (IMPORTANT - do this first!)
devtools::load_all("../mouseExperiment")
```

### 2. Install Dashboard Dependencies

```r
# Install required packages if you don't have them
install.packages(c("shiny", "shinycssloaders", "shinyFeedback", "shinyjs", "DT", 
                   "dplyr", "ggplot2", "scales", "tidyr"))
```

### 3. Load the Dashboard in Development

```r
# Navigate to dashboard directory
setwd("/home/david/Projects/mouseExperimentDashboard")

# Load the dashboard package
devtools::load_all()

# Run the dashboard
mouseExperimentDashboard::runDashboard()
```

This will:
1. Verify mouseExperiment is loaded (will error if not)
2. Load all functions from R/
3. Source global.R (loads libraries and initializes)
4. Launch the Shiny app in your default browser at http://localhost:3838

### Quick Reference: Complete Setup from Scratch

```r
# Step 1: Load mouseExperiment (parent package)
devtools::load_all("/home/david/Projects/mouseExperiment")

# Step 2: Navigate to dashboard
setwd("/home/david/Projects/mouseExperimentDashboard")

# Step 3: Install dependencies (one time only)
install.packages(c("shiny", "shinycssloaders", "shinyFeedback", "shinyjs", "DT", 
                   "dplyr", "ggplot2", "scales", "tidyr"))

# Step 4: Load dashboard
devtools::load_all()

# Step 5: Run the app
mouseExperimentDashboard::runDashboard()
```

### 3. Alternative: Direct Shiny Execution

If you already have both packages loaded in your R session:

```r
shiny::runApp("/home/david/Projects/mouseExperimentDashboard/app_standalone.R")
```

## Project Structure

```
mouseExperimentDashboard/
├── app_standalone.R              # Standalone development entry point
├── global.R                      # Global setup (libraries, constants)
├── mouseExperimentDashboard.Rproj# RStudio project file
├── DESCRIPTION                   # Package metadata
├── LICENSE                       # MIT License
├── NAMESPACE                     # Package exports
├── R/
│   ├── app.R                    # Main runDashboard() function
│   ├── ui_modules.R             # All UI functions (UI definitions)
│   ├── server_modules.R         # All server logic (reactive functions)
│   └── helpers.R                # Utility functions
├── www/
│   └── css/
│       └── custom.css           # Custom styling
└── inst/
    └── sample_data/             # Sample datasets for testing
```

## Key Features Implemented

### 1. Data Upload Module
- ✅ File upload (CSV, RDA)
- ✅ Automatic column detection
- ✅ Manual column mapping
- ✅ Data preview (first 20 rows)
- ✅ Sample dataset loading

**Files:** `R/ui_modules.R` (dataUploadUI), `R/server_modules.R` (data_upload_server)

### 2. Tumor Growth Analysis
- ✅ Parameter configuration UI
- ✅ Model type selection (LME4, AUC)
- ✅ Model options (transformation, polynomial degree, random effects)
- ✅ Results tabs (Summary, ANOVA, Comparisons, Effects, Growth Rates, Plots, Diagnostics)
- ✅ Download buttons for each result type
- ✅ Error handling

**Files:** `R/ui_modules.R` (tumorGrowthUI), `R/server_modules.R` (tumor_growth_server)

## Development Workflow

### Adding a New Feature

1. **Define UI** in `R/ui_modules.R`
   ```r
   newFeatureUI <- function() {
     fluidPage(
       # Your UI code here
     )
   }
   ```

2. **Add Server Logic** in `R/server_modules.R`
   ```r
   new_feature_server <- function(input, output, session, rv) {
     # Your reactive logic here
     return(rv)
   }
   ```

3. **Add Tab to Main App** in `R/app.R`
   ```r
   tabPanel("New Feature", newFeatureUI())
   ```

4. **Initialize Server** in `server` function
   ```r
   rv <- new_feature_server(input, output, session, rv)
   ```

5. **Test with sample data**
   ```r
   devtools::load_all()
   mouseExperimentDashboard::runDashboard()
   ```

### Common Debugging

**Issue: "could not find function" error**
- Run `devtools::load_all()` to reload packages
- Check that function is in R/ folder

**Issue: Shiny app won't start**
- Check console for error messages
- Verify all dependencies are installed: `library(shiny)`, etc.

**Issue: Data not appearing in dropdowns**
- Check that data was successfully loaded in Data Upload tab
- Verify column detection worked (check "File Information" section)

## Dependencies

### Core Shiny Packages
- **shiny** - Web application framework
- **shinycssloaders** - Loading spinners
- **shinyFeedback** - Validation feedback
- **shinyjs** - JavaScript interactions

### Data & Visualization
- **dplyr** - Data manipulation
- **ggplot2** - Plotting (via mouseExperiment)
- **DT** - Data tables
- **scales** - Scale helpers

### Analysis
- **mouseExperiment** (>= 0.3.0) - Core analysis functions

## Testing

### Manual Testing Checklist

- [ ] Data upload (CSV file)
- [ ] Column detection (automatic vs. manual)
- [ ] Sample data loading
- [ ] Tumor growth analysis (LME4)
- [ ] Tumor growth analysis (AUC)
- [ ] Results appear in tabs
- [ ] Download buttons work
- [ ] Error handling for invalid data

### Automated Testing (Future)

```r
# Once tests are written:
devtools::test()

# With coverage:
covr::package_coverage()
```

## Common Tasks

### Update Column Detection Patterns

Location: `R/helpers.R` - `detect_columns()` function

```r
detect_column <- function(col_names_lower, patterns) {
  # Add new patterns to detect additional column names
  patterns <- c("id", "mouse_id", "animal_id", "YOUR_PATTERN_HERE")
}
```

### Add New Sample Dataset

1. Add to SAMPLE_DATASETS in `global.R`
2. Implement loading logic in `data_upload_server()` function

### Customize Styling

Edit `www/css/custom.css` - changes apply immediately on page refresh

### Add New Analysis Module (e.g., Survival)

1. Create `survivalUI()` in `R/ui_modules.R`
2. Create `survival_server()` in `R/server_modules.R`
3. Add tab to main app in `R/app.R`
4. Add import to `global.R` if needed

## Performance Tips

1. **Server-side filtering**: Use reactiveValues() to cache data
2. **Avoid excessive re-rendering**: Use eventReactive() for computationally expensive operations
3. **Disable buttons during processing**: Use `shinyjs` to manage button states
4. **Progress bars**: Use `withProgress()` for long operations

## Next Steps

Current MVP (Complete):
- ✅ Data upload & configuration
- ✅ Tumor growth analysis
- ✅ Results visualization

Planned Features (Phase 2):
- [ ] Survival analysis module
- [ ] Synergy analysis module
- [ ] Power analysis module
- [ ] Report generation
- [ ] Session save/load
- [ ] Batch processing

## Resources

- [Shiny Documentation](https://shiny.rstudio.com/)
- [Shiny Modules](https://shiny.rstudio.com/articles/modules.html)
- [mouseExperiment Documentation](../mouseExperiment/README.md)
- [Shiny Tutorial](https://shiny.rstudio.com/tutorial/)

## Questions or Issues?

See the Help tab in the dashboard for user-facing documentation, or check console error messages for development debugging.
