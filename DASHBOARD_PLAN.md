## Plan: Interactive Dashboard for mouseExperiment

### TL;DR

Build a **Shiny web application** that provides researchers with a visual, point-and-click interface to:
- Upload and prepare tumor growth experimental data
- Execute statistical analyses (tumor growth, survival, synergy, power analysis)
- Interactively explore results with dynamic visualizations
- Generate downloadable reports

The dashboard would reduce the learning curve for R code while leveraging all existing package functionality. Estimated scope: 2-3 months for MVP, 4-6 months for production release. Core decision: Use R Shiny (vs. standalone web app) to minimize implementation complexity and maximize code reuse.

---

## 1. Architecture Overview

### 1.1 Technology Stack

| Layer | Technology | Rationale |
|-------|-----------|-----------|
| **Backend** | R + mouseExperiment pkg | Package already has all analysis logic |
| **Frontend Framework** | Shiny (R) | Native R, minimal learning curve, rapid development |
| **Styling** | Bootstrap 4 (built-in to Shiny) | Professional, responsive, no dependencies |
| **Data Viz** | ggplot2 (already in pkg) | Consistent with existing API |
| **Deployment** | Shiny Server / Posit Cloud | Simple R deployment, no infrastructure management |
| **Database** | SQLite (optional) | Light session management, data versioning |

### 1.2 System Architecture

```
┌─────────────────────────────────────────┐
│     USER INTERFACE (Shiny Browser)      │
│  ┌─────────────────────────────────────┐│
│  │ Step 1: Data Upload & Preview       ││
│  │ Step 2: Parameter Configuration     ││
│  │ Step 3: Run Analysis / Visualization││
│  │ Step 4: Explore & Export Results    ││
│  └─────────────────────────────────────┘│
└──────────────────├──────────────────────┘
                   │
           Shiny Server (R)
          (reactive functions)
                   │
┌──────────────────▼──────────────────────┐
│   mouseExperiment R Package              │
│  ┌─────────┬──────┬─────┬────┬────────┐ │
│  │Data     │Stats │Surv │Syn │Power   │ │
│  │Prep     │Tumor │Anal │Erg │Anal    │ │
│  └─────────┴──────┴─────┴────┴────────┘ │
└──────────────────┬──────────────────────┘
                   │
        ggplot2 Visualizations
        Statistical Results
        Tabular Data
```

---

## 2. Feature Scope (MVP - Minimum Viable Product)

### Phase 1: Core Analysis Workflows (MVP - 2-3 months)

**Dashboard Modules:**

#### **Module 1: Data Upload & Preparation**
- File upload (CSV, RDA)
- Data preview (first 20 rows with scrolling)
- Column mapping interface (select which columns represent ID, Day, Volume, Treatment, Cage, etc.)
- Automatic data validation with clear error messages
- Volume calculation tool (dropdown to select formula)
- Date parsing (auto-detect format or manual specification)
- Data summary statistics

**Key Features:**
```
├─ Upload file (CSV/RDA)
├─ Preview data with column detection
├─ Map columns to expected names
│  ├─ ID column
│  ├─ Day/Time column
│  ├─ Volume column OR Length+Width for calculation
│  ├─ Treatment column
│  ├─ Cage column (optional)
│  └─ Date column (optional)
├─ Configure volume formula (ellipsoid, modified, cylinder, sphere, box)
├─ Parse dates with format specification
├─ Display validation results
└─ "Ready to Analyze" confirmation
```

#### **Module 2: Tumor Growth Analysis**
- Primary analysis interface for tumor_growth_statistics()
- Parameter controls:
  - Model type selector (LME4 vs. AUC)
  - Transformation selector (log, sqrt, none)
  - Polynomial degree slider (1-3)
  - Random effects structure (intercept_only, slope, none)
  - Cage effect handling dropdown
  - Reference group selector
  - Optional: cage collinearity testing checkbox

**Output:**
```
├─ Model Summary (expandable tabsets)
│  ├─ ANOVA table
│  ├─ Parameter estimates
│  └─ Model diagnostics
├─ Pairwise Comparisons (table with significance stars)
├─ Treatment Effects (forest-style visualization)
├─ Growth Rates by Subject (interactive table)
└─ Interactive Plots (tabs)
   ├─ Tumor Growth Curves
   ├─ Growth Rate Box Plot
   ├─ AUC Comparison
   └─ Caterpillar Plot (random effects)
```

#### **Module 3: Survival Analysis**
- Interface for survival_statistics()
- Parameter controls:
  - Reference group selector
  - Firth correction toggle
  - Column mapping for time and event indicator

**Output:**
```
├─ Method Used (displays which approach: Cox/Coxph/LogRank)
├─ Results Table
│  ├─ Hazard Ratios with 95% CI
│  ├─ P-values
│  ├─ Median Survival Times
│  └─ Events/Total counts
├─ Kaplan-Meier Curves
└─ Forest Plot (hazard ratios)
```

#### **Module 4: Synergy Analysis**
- Interface for analyze_drug_synergy() + analyze_drug_synergy_over_time()
- Parameter controls:
  - Drug A name input
  - Drug B name input
  - Combination name input
  - Control group name input
  - Evaluation timepoint selector (dropdown or slider)

**Output:**
```
├─ Synergy Summary Table
│  ├─ TGI for each treatment
│  ├─ Bliss Expected vs. Observed
│  ├─ Combination Index (CI)
│  └─ Synergy determination (synergistic/additive/antagonistic)
├─ Bliss Plot (time-course)
├─ Combination Index Plot
└─ Synergy Trend Visualization
```

#### **Module 5: Power Analysis**
- Interface for post_power_analysis()
- Parameter controls:
  - Alpha levels checkboxes (0.05, 0.01, etc.)
  - Target power slider (0.70-0.99)
  - Method selector (parametric/simulation/auc)
  - N simulations slider (for simulation method)

**Output:**
```
├─ Effect Sizes Table (by treatment pair)
├─ Power Analysis Results
│  ├─ Power curves for different effect sizes
│  ├─ Power vs. sample size plots
│  └─ Sample size recommendations
└─ Download recommendations as CSV
```

#### **Module 6: Dose-Response Analysis** (Optional/Phase 2)
- Interface for dose_response_statistics()
- Interactive dose-response curves
- Trend testing visualization

---

### Phase 2: Extended Features (Post-MVP - months 4-6)

1. **Report Generation**
   - PDF/HTML report generation with standardized template
   - Include selected plots, tables, interpretation
   - Customizable title, author, date

2. **Batch Processing**
   - Upload multiple files, run same analysis on all
   - Comparative results across datasets
   - Export summary comparison tables

3. **Data Management**
   - Save analysis sessions (SQL backend)
   - Load previous analyses for modification
   - Compare multiple analyses side-by-side

4. **Advanced Visualization**
   - Interactive ggplotly plots (hover tooltips, zoom, pan)
   - 3D survival curves (plotly)
   - Heatmaps for multiple comparisons

5. **Help & Documentation**
   - Embedded vignettes per module
   - Tooltips for all parameters
   - Example dataset dropdown to tutorial analyses

6. **Export Options**
   - Download individual plots (PNG, PDF, SVG)
   - Download tables as CSV/Excel
   - Download R code to reproduce analysis
   - R Markdown script for reproducibility

---

## 3. Detailed Implementation Steps

### Step 1: Project Setup (Week 1)

**Create Shiny project structure:**
```
mouseExperiment-dashboard/
├── app.R                          # Main Shiny app
├── global.R                       # Load libraries, data, helpers
├── R/
│   ├── ui_modules.R              # UI functions
│   ├── server_modules.R          # Server logic
│   ├── helpers.R                 # Utility functions
│   └── validation.R              # Input validation
├── www/
│   ├── css/
│   │   └── custom.css           # Custom styling
│   └── images/
│       └── logo.png             # Package logo
├── inst/
│   ├── sample_data/             # Example datasets
│   └── templates/               # Report templates
├── tests/                        # Unit tests for dashboard
└── README.md                     # Dashboard documentation
```

**Decision:** Use modular Shiny design with separate UI/server modules to improve maintainability.

### Step 2: UI Framework (Weeks 1-2)

**Main Navigation Structure:**
```R
ui <- navbarPage(
  "mouseExperiment Dashboard",
  tabPanel("Home", homeUI()),  
  tabPanel("Data Upload", dataUploadUI()),
  tabPanel("Tumor Growth", tumorGrowthUI()),
  tabPanel("Survival", survivalUI()),
  tabPanel("Synergy", synergyUI()),
  tabPanel("Power Analysis", powerUI()),
  tabPanel("Help", helpUI())
)
```

**Key UI Components:**
- Sidebar for navigation and global settings
- Main content area with tabsets for outputs
- Status indicators (processing, error messages)
- Progress bars for long computations
- Modal dialogs for parameter configuration

### Step 3: Data Upload Module (Weeks 2-3)

**Server Logic:**
```r
# Reactive data frame with uploaded file
uploaded_data <- reactive({ ... })

# Validation and preview
data_preview <- reactive({ head(uploaded_data(), 20) })

# Column type detection
detected_columns <- reactive({
  # Auto-detect ID, Day, Volume, Treatment, Cage patterns
})

# Validation messages
validation_status <- reactive({
  # Check required columns, data types, value ranges
})
```

**UI Elements:**
- File input widget
- Data preview table
- Column mapping dropdowns
- Validation error/warning displays
- "Proceed to Analysis" button (disabled until valid)

### Step 4: Tumor Growth Module (Weeks 3-5)

**Server Logic:**
```r
# Parameter inputs
model_type <- reactive({ input$model_type })
transform <- reactive({ input$transform })
cage_effect <- reactive({ input$cage_effect })

# Execute analysis
tumor_results <- eventReactive(input$run_analysis, {
  tryCatch({
    tumor_growth_statistics(
      df = uploaded_data(),
      model_type = model_type(),
      ...
    )
  }, error = function(e) {
    return(list(error = as.character(e)))
  })
})

# Generate plots
plot_growth <- reactive({
  plot_growth_rate(tumor_results()$growth_rates)
})

plot_auc <- reactive({
  plot_auc(tumor_results()$auc_analysis$individual)
})
```

**UI Elements:**
- Expandable parameter panels (collapsible containers)
- Slider/checkbox controls for each parameter
- "Run Analysis" button with loading progress bar
- Tabset with multiple output visualizations
- Downloadable plot buttons
- Results summary table with export to CSV

### Step 5: Survival Module (Weeks 5-6)

**Server Logic:**
```r
survival_results <- eventReactive(input$run_survival, {
  tryCatch({
    survival_statistics(
      df = uploaded_data(),
      time_column = input$time_column,
      ...
    )
  }, error = function(e) list(error = as.character(e)))
})

# Generate KM curves and forest plots
plot_km <- reactive({
  survfit(Surv(...) ~ Treatment, data = uploaded_data())
})
```

**UI Elements:**
- Column selectors for time and event
- Reference group dropdown
- Results table (HR, CI, p-value, median survival)
- KM curves visualization
- Forest plot for hazard ratios

### Step 6: Synergy Module (Weeks 6-7)

**Server Logic:**
```r
synergy_results <- eventReactive(input$run_synergy, {
  analyze_drug_synergy(
    df = uploaded_data(),
    drug_a_name = input$drug_a,
    drug_b_name = input$drug_b,
    combo_name = input$combo
  )
})

synergy_overtime <- eventReactive(input$run_synergy_overtime, {
  analyze_drug_synergy_over_time(...)
})
```

**UI Elements:**
- Treatment name text inputs (with tab-completion from data)
- Timepoint selector
- Synergy summary table
- Bliss and CI visualizations
- Synergy interpretation text

### Step 7: Power Analysis Module (Weeks 7-8)

**Server Logic:**
```r
power_results <- eventReactive(input$run_power, {
  post_power_analysis(
    data = uploaded_data(),
    alpha = input$alpha_levels,
    power = input$target_power,
    method = input$power_method
  )
})
```

**UI Elements:**
- Checkboxes for alpha levels
- Power slider
- Method radio buttons
- Power curves (interactive with plotly)
- Sample size recommendations table

### Step 8: Report Generation (Weeks 8-10)

**Implementation:**
- Create R Markdown template with placeholder sections
- Dynamically populate with:
  - Data summary
  - Selected analysis results
  - Generated plots
  - Interpretation text
- Button to "Download PDF Report"

**Using:**
```r
# rmarkdown::render() with params
rmarkdown::render(
  "template.Rmd",
  params = list(
    title = input$report_title,
    analysis_results = results,
    plots = list(...)
  ),
  output_file = "report.pdf"
)
```

### Step 9: Testing & QA (Weeks 10-11)

- Unit tests for R helper functions
- Integration tests with sample datasets
- UI/UX testing with domain experts
- Performance testing with large datasets
- Error handling and edge cases

### Step 10: Deployment & Documentation (Week 12)

- Deploy to Shiny Server OR Posit Cloud
- Create user documentation/manual
- Record tutorial videos
- Deploy as R package `mouseExperimentDashboard` on GitHub

---

## 4. Technical Details & Decisions

### 4.1 Reactive Programming Strategy

**Reactive Pattern:**
```r
# Minimal reactives - only what's needed
uploaded_data <- reactive({ ... })
result_analysis <- eventReactive(input$run, { ... })
plot_output <- reactive({ generate_plot(result_analysis()) })

# Avoid over-reactivity that causes unnecessary recalculation
```

**Error Handling:**
```r
safeAnalysis <- function(df, ...) {
  tryCatch({
    tumor_growth_statistics(df, ...)
  }, error = function(e) {
    list(error = TRUE, message = conditionMessage(e))
  })
}
```

### 4.2 Performance Optimization

**For Fast Response:**
- Cache data frames (avoid re-reading uploads)
- Memoize heavy computations
- Use `future` package for non-blocking simulations (power analysis)

**Handling Long Computations:**
- Progress bars for slow analyses (power = "simulation")
- Background processing option (Shiny background jobs)
- Disable run button during execution

### 4.3 Data Security (For Future)

**Considerations:**
- File upload sanitization (no path traversal attacks)
- Session isolation (each user gets private session)
- Optional: Authenticate users (not necessary for MVP)
- Optional: Data encryption at rest (for clinical data)

### 4.4 Code Reuse from Package

**Strategy:**
- No code duplication - call package functions directly
- Package maintains analysis logic
- Dashboard adds UI/UX layer only
- Easy to keep dashboard updated as package evolves

---

## 5. UI/UX Design Principles

### Layout Philosophy

```
┌────────────────────────────────────────┐
│  Header: Logo | Title | Help           │
├─────────────────┬──────────────────────┤
│  SIDEBAR        │  MAIN CONTENT        │
│  ┌────────────┐ │  ┌────────────────┐  │
│  │ Navigation │ │  │ Analysis       │  │
│  │ Controls   │ │  │ Results &      │  │
│  │ Parameters │ │  │ Visualizations │  │
│  │            │ │  │                │  │
│  └────────────┘ │  └────────────────┘  │
├─────────────────┴──────────────────────┤
│  Footer: Version | Feedback            │
└────────────────────────────────────────┘
```

### Key UX Decisions

**1. Wizard Pattern (Optional)**
- Could use a multi-step wizard for first-time users
- Step 1: Upload data
- Step 2: Configure columns
- Step 3: Select analysis
- Step 4: Set parameters
- Step 5: Review results

**2. Contextual Help**
- Tooltips on hover for all parameters
- "Learn more" links to documentation
- Embedded example analyses with sample data

**3. Input Validation Feedback**
- Real-time validation (column detection)
- Clear error messages with suggestions
- Highlight invalid fields in red
- Provide "Start Over" option after errors

**4. Results Organization**
- Always show summary first (most important finding)
- Detailed tables below (expand/collapse)
- Plots in separate tabset
- Export options for each component

**5. Color & Accessibility**
- Use colorblind-safe palettes
- High contrast for readability
- Clear button labels (not just icons)
- Responsive design for tablets

---

## 6. Integration with Existing Package

### 6.1 Dependency Requirements

```r
# In DESCRIPTION for mouseExperimentDashboard package
Imports:
    mouseExperiment (>= 0.3.0),
    shiny (>= 1.7.0),
    shinycssloaders (>= 1.0.0),
    shinyFeedback (>= 0.4.0),
    tidyverse,
    rmarkdown,
    knitr,
    plotly,
    DT
```

### 6.2 Package as Dependency

**Strategy:**
- Create separate R package: `mouseExperimentDashboard`
- Depends on `mouseExperiment` package
- Installation: `devtools::install_github("user/mouseExperimentDashboard")`
- Avoids bloating main package with Shiny dependency

**Alternatively:**
- Keep dashboard in `/inst/shiny/app.R` within mouseExperiment
- Use `mouseExperiment::runDashboard()` function to launch

### 6.3 Version Compatibility

- Dashboard should version-lock to mouseExperiment API
- Add version checks at startup
- Display warning if package/dashboard versions mismatch

---

## 7. Deployment Options

### Option A: Shiny Server (Local/On-premises)
- Pros: Full control, no monthly cost, can handle internal data
- Cons: Infrastructure management required
- **Cost:** Server + admin time
- **Setup:** 1-2 hours

### Option B: Posit Cloud (formerly shinyapps.io)
- Pros: No infrastructure, simple deployment, SSL included
- Cons: Monthly subscription ($9-999/month depending on usage)
- **Cost:** Data-dependent
- **Setup:** 30 minutes (GitHub-based deployment)

### Option C: Docker Container
- Pros: Reproducible environments, scalable
- Cons: Requires Docker infrastructure
- **Cost:** Cloud provider dependent ($10-50/month minimum)
- **Setup:** 2-3 hours

**Recommendation:** Start with Posit Cloud free tier (25 active hours/month), upgrade as user base grows.

---

## 8. Implementation Timeline

```
MONTHS 1-2 (MVP):
├─ Week 1-2: Project setup, UI framework, data upload
├─ Week 3-5: Tumor growth analysis module  
├─ Week 5-6: Survival analysis module
├─ Week 6-7: Synergy analysis module
├─ Week 7-8: Power analysis module
├─ Week 9-10: Testing, QA, bug fixes
├─ Week 11-12: Deployment preparation
└─ Deliverable: Functional MVP deployed to Posit Cloud

MONTHS 3-4 (Polish & Testing):
├─ Report generation
├─ User testing with domain experts
├─ Documentation & tutorial videos
├─ Performance optimization
└─ Deliverable: Production-ready v1.0

MONTHS 4-6 (Extended Features - Optional):
├─ Batch processing
├─ Session management
├─ Interactive plots (plotly)
├─ Advanced visualization
└─ Deliverable: v2.0 with enhanced features
```

---

## 9. Success Criteria

### MVP Success (v1.0)
- ✅ Upload CSV/RDA files with validation
- ✅ Run tumor growth analysis from UI
- ✅ Run survival analysis from UI  
- ✅ Run synergy analysis from UI
- ✅ Run power analysis from UI
- ✅ Download plots/tables from all analyses
- ✅ Error handling with helpful messages
- ✅ Deployed and accessible
- ✅ Basic documentation

### Production Success (v2.0)
- ✅ PDF report generation
- ✅ Session save/load
- ✅ Batch processing
- ✅ Interactive visualizations
- ✅ >90% code coverage
- ✅ Zero critical bugs
- ✅ <2 second load time
- ✅ <10 second analysis runtime (or progress bar shown)
- ✅ Comprehensive user documentation

---

## 10. Risk Analysis & Mitigation

| Risk | Probability | Impact | Mitigation |
|------|-------------|--------|-----------|
| Package API changes break dashboard | Medium | High | Version lock dependencies, add API tests |
| Large dataset performance issues | Medium | High | Test with large datasets early, optimize reactives |
| Users confused by statistical parameters | High | Medium | Add extensive tooltips, example workflows, tutorials |
| File upload security issues | Low | High | Validate uploads, run in sandbox, security audit |
| Complex UI overwhelms users | Medium | Medium | Start simple, add complexity incrementally, user testing |
| Deployment downtime | Low | Medium | Use Posit Cloud (99.9% SLA), maintain version control |

---

## 11. Future Enhancements (Phase 3+)

1. **Collaborative Features**
   - Share analyses with team members
   - Comments on results/visualizations
   - Version history with diffs

2. **Advanced Modeling**
   - Custom model specification
   - Model comparison (AIC/BIC)
   - Interaction terms UI

3. **Machine Learning**
   - Predictive modeling on tumor growth
   - Clustering of similar response profiles
   - Variable importance visualization

4. **Real Data Integration**
   - Direct database connections
   - API endpoints for programmatic access
   - Export to clinical trial systems

5. **Mobile App**
   - React Native wrapper
   - Offline analysis capability
   - Push notifications for long jobs

---

## 12. Recommendations & Next Steps

### Immediate Actions (This Week)

1. **Clarify Scope with Stakeholders**
   - Which analyses are highest priority?
   - Who are the primary users?
   - What's the timeline/budget?
   - Deployment environment (public/private)?

2. **Prototype Decision**
   - Build quick Shiny proof-of-concept with just data upload + tumor growth analysis
   - Get feedback from potential users
   - Validate technical approach

3. **Resource Planning**
   - Assign 1 FTE developer for MVP (2-3 months)
   - Designate domain expert reviewer (10-15 hrs/month)
   - User testing group (5-10 researchers)

### For Development

**Technology Stack Confirmed:**
- ✅ R Shiny (best fit for R package ecosystem)
- ✅ mouseExperiment as dependency (no duplication)
- ✅ Posit Cloud deployment (simplest for startup)
- ✅ Modular UI/server design (maintainability)

**Development Best Practices:**
- Use Git version control from day 1
- Write tests as you build (not after)
- Get frequent user feedback
- Performance test with real data early
- Document parameters and decisions

---

## Conclusion

A **Shiny dashboard for mouseExperiment is highly feasible** and would dramatically improve accessibility for non-R-proficient researchers. The package's excellent modular design and comprehensive documentation make it an ideal candidate for dashboard-ification.

**Key Advantages:**
- Leverages all existing package functionality without duplication
- Rapid development (Shiny reduces time-to-MVP)
- Easy to keep updated as package evolves
- Familiar to R developers (not a black box)
- Can be deployed in various environments
- Natural upgrade path for R users

**Recommend proceeding to Phase 1 with:**
1. Detailed user persona research
2. Shiny prototype (2-week spike)
3. MVP development plan (12-week timeline)
4. Deployment to Posit Cloud
