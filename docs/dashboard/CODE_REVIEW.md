# Code Review: mouseExperimentDashboard

**Date:** March 25, 2026  
**Package:** mouseExperimentDashboard v0.1.0  
**Scope:** Deep review of all 9 R source files + global.R (~6,700 lines)

---

## Executive Summary

mouseExperimentDashboard is a Shiny dashboard that provides a GUI for the mouseExperiment analysis package, covering tumor growth, survival, dose-response, drug synergy, and post-hoc power analysis. The dashboard is feature-rich with extensive user controls, dynamic plot customization, and downloadable outputs. However, it has significant **architectural debt** (no Shiny module namespacing, 1,850-line god function), widespread **code duplication** (date/volume preprocessing repeated 4×), **stale download handlers**, and several **functional bugs**.

**Overall assessment:** Usable but fragile. Prioritize DRY refactoring and download handler fixes.

---

## 1. Bugs

### D1: QQ Plot Stored as "Caterpillar Plot" — LOW

**File:** server_modules.R  
**Impact:** Misleading variable name; downstream code expecting a caterpillar plot gets QQ instead

```r
result$caterpillar_plot <- qq_p  # QQ plot assigned to caterpillar_plot slot
```

The caterpillar plot from the package is not stored in the result; the QQ plot is stored under its name. If UI code references `caterpillar_plot` expecting the actual caterpillar plot, users see a QQ diagnostic plot instead.

---

### D2: Forest Plot CI Uses Hardcoded 1.96 × SE — MEDIUM

**File:** server_modules.R (effect forest plot rendering)  
**Impact:** Confidence intervals ignore actual degrees of freedom

```r
ci_lower = Estimate - 1.96 * SE,
ci_upper = Estimate + 1.96 * SE
```

For LME4 models with small sample sizes, the appropriate multiplier is `qt(0.975, df)` where `df` comes from the model's Satterthwaite or Kenward-Roger approximation. Using 1.96 (the z-value) underestimates CI width for small samples.

---

### D3: Download Handler Opens File Before Null Check — LOW

**File:** server_modules.R (waterfall plot download handler)  
**Impact:** Creates blank PNG files when no plot exists

```r
png(file, width = 10, height = 6, units = "in", res = 300)
# ... then checks if plot is NULL
```

Should check for NULL before calling `png()`.

---

### D4: Synergy Download Handlers Save Stale Plots — MEDIUM

**File:** synergy_modules.R  
**Impact:** Downloaded plots don't reflect user's color/legend/layout customizations

Download handlers save eagerly-built plots from `rv$syn_results` rather than the reactively-rendered versions that incorporate user's dynamic UI customizations (colors, legend position, etc.). The downloaded file shows the default plot, not what the user sees on screen.

---

### D5: `rv$syn_results$plot_trend` Never Set — LOW

**File:** synergy_modules.R  
**Impact:** "Download Trend Plot" always produces empty file

The download handler references `rv$syn_results$plot_trend`, but the synergy analysis code never assigns this slot. The download silently produces a blank file.

---

### D6: Survival `input$surv_palette` Referenced But Never Defined — LOW

**File:** survival_modules.R  
**Impact:** Warnings in reactive expressions; palette defaults to NULL

`input$surv_palette` is used in the KM plot rendering code but no corresponding UI input element exists in ui_modules.R.

---

### D7: Power Curves Hardcode `sig.level = 0.05` — MEDIUM

**File:** power_modules.R  
**Impact:** Power curve visualization ignores user-selected alpha level

The `renderPlot` for power curves always uses `sig.level = 0.05` regardless of what the user selected in the UI controls. The stored results use the correct alpha, but the curve regeneration does not.

---

### D8: Survival Forest Plot Double Computation — LOW

**File:** survival_modules.R  
**Impact:** Unnecessary duplicate computation (performance, not correctness)

`renderPlot` builds `plot_df`, reference row, and axis limits, then calls `build_surv_forest_plot()` which internally recomputes the same reference row and axis limits.

---

### D9: Hardcoded Legend Position in Synergy — LOW

**File:** synergy_modules.R  
**Impact:** User's legend position selection (`syn_opt_legend_pos`) is ignored

The Bliss/CI plots use a hardcoded legend position instead of reading from the user's dropdown selection.

---

## 2. Code Quality Issues

### 2.1 God Function: `tumor_growth_server()`

server_modules.R contains `tumor_growth_server()` at ~1,850 lines. This single function includes:
- Data upload and validation
- Date/volume preprocessing
- Treatment schedule merging
- LME4/AUC analysis dispatch
- 25+ `output$*` render definitions
- 10+ download handlers
- Plot customization reactives
- Info box rendering
- Growth rate, AUC, forest, waterfall, caterpillar, diagnostics rendering

**Recommendation:** Break into Shiny modules using `moduleServer()` + `NS()`:
- `mod_upload_server()` — file upload and validation
- `mod_tg_analysis_server()` — analysis dispatch
- `mod_tg_plots_server()` — all TG plot rendering
- `mod_tg_downloads_server()` — download handlers

---

### 2.2 No Shiny Module Namespacing

The entire app uses flat input/output IDs with no `NS()` namespacing. All IDs must be globally unique across all tabs. This makes:
- Name collisions likely as the app grows
- Refactoring into separate files risky
- Testing individual components impossible

Examples of collision-prone IDs: `input$ref_group`, `input$analysis_type`, `input$legend_pos` (used in TG), vs `input$surv_ref_group`, `input$surv_legend_pos` (used in survival) — manual prefix convention instead of proper NS().

---

### 2.3 Date/Volume Preprocessing Duplicated 4×

The same date-format detection + volume calculation + Day column derivation logic appears in:
1. **server_modules.R** — tumor growth tab
2. **survival_modules.R** — survival tab  
3. **synergy_modules.R** — synergy tab
4. **dose_response_modules.R** — dose-response tab

Each copy has ~30–50 lines of nearly identical code doing:
```r
# Detect date format (try multiple formats)
# Calculate Volume if L/W columns present
# Calculate Day from dates
# Handle treatment schedule merge
```

**Recommendation:** Create `preprocess_uploaded_data()` in helpers.R that all modules call. Accept the raw uploaded data and return a standardized data frame with Volume, Day, and Treatment columns.

---

### 2.4 Debug `cat()` Statements in Production

~15+ `cat("DEBUG: ...")` statements scattered across production code:

| File | Approximate Count |
|------|-------------------|
| server_modules.R | ~10 |
| survival_modules.R | ~5 |
| synergy_modules.R | ~3 |
| dose_response_modules.R | ~2 |

These print to the R console during production use — invisible to end users but noisy for anyone monitoring the server log.

**Recommendation:** Remove all debug `cat()` calls, or replace with a proper logging framework (e.g., `futile.logger` or `logger` package) that supports log levels.

---

### 2.5 Download Handlers Don't Match Rendered Plots

A systematic issue across multiple modules: download handlers rebuild plots from scratch (or use stale cached versions) instead of saving the plot that the user actually sees. Affected:

| Module | Issue |
|--------|-------|
| synergy_modules.R | Downloads save `rv$syn_results$plot_*` (initial), not the `renderPlot` version with user colors/legend |
| power_modules.R | Downloads save `rv$pw_results$plots$power_curves` but `renderPlot` builds its own curves with different parameters |
| survival_modules.R | KM plot download handler reconstructs the plot instead of using the rendered version |
| dose_response_modules.R | Download handlers reconstruct plots from scratch |

**Recommendation:** Use `reactiveVal` to store the current plot object, then have both `renderPlot` and the download handler reference the same reactive.

---

### 2.6 Triple-Colon Access to Unexported Functions

dose_response_modules.R uses `mouseExperiment:::prepare_dose_data()` — the triple-colon operator accesses unexported (internal) functions. This:
- Breaks if the internal function is renamed/removed
- Generates R CMD check warnings
- Is not a stable API

**Recommendation:** Export `prepare_dose_data()` from mouseExperiment (add `@export` roxygen tag), or duplicate the needed logic locally.

---

### 2.7 Deprecated `aes_string()` in dose_response_modules.R

`aes_string()` is used throughout the dose-response plotting code. Deprecated since ggplot2 3.0.0 (2018). Should migrate to `aes(.data[[col]])`.

---

### 2.8 Redundant `library()` Calls

global.R contains multiple `library()` calls:
```r
library(shiny)
library(shinyjs)
library(colourpicker)
# etc.
```

For a package-based Shiny app, these should be handled through the DESCRIPTION Imports field and `@import`/`@importFrom` directives. The current approach works but creates implicit dependencies not tracked by the package system.

---

### 2.9 `%||%` Operator Redefined

dose_response_modules.R defines:
```r
`%||%` <- function(a, b) if (!is.null(a)) a else b
```

This is already available from `rlang::%||%` (which is in mouseExperiment's Imports). Redefining it risks subtle differences and masks the standard implementation.

---

### 2.10 `shinyjs::useShinyjs()` Called Twice

Both global.R and app.R UI include `shinyjs::useShinyjs()`. Only one call is needed in the UI.

---

## 3. Architecture Issues

### 3.1 File Organization

Current structure mixes concerns:

| File | Lines | Responsibilities |
|------|-------|-----------------|
| server_modules.R | 2,396 | TG upload, TG analysis, TG 25+ renders, TG 10+ downloads |
| survival_modules.R | 1,007 | Surv upload, analysis, renders, downloads |
| ui_modules.R | 817 | All UI definitions for all tabs |
| synergy_modules.R | 683 | Synergy upload, analysis, renders, downloads |
| helpers.R | 644 | Shared utilities |
| dose_response_modules.R | 603 | DR upload, analysis, renders, downloads |
| power_modules.R | 409 | Power upload, analysis, renders, downloads |
| app.R | 260 | Main app shell, JS injection |
| zzz.R | 19 | `.onLoad` hook |

**Proposed restructuring:**

```
R/
  app.R                    # Main app (keep)
  helpers.R                # Shared utilities (keep)
  preprocessing.R          # NEW: Consolidated date/volume/day preprocessing
  mod_tumor_growth.R       # TG server module (from server_modules.R)
  mod_tumor_growth_ui.R    # TG UI module (from ui_modules.R)
  mod_survival.R           # Survival module
  mod_survival_ui.R        # Survival UI
  mod_synergy.R            # Synergy module  
  mod_synergy_ui.R         # Synergy UI
  mod_dose_response.R      # DR module
  mod_dose_response_ui.R   # DR UI
  mod_power.R              # Power module
  mod_power_ui.R           # Power UI
  zzz.R                    # Keep
```

### 3.2 Reactive Value Sprawl

Each module stores results in a shared `rv` (reactiveValues) object with deeply nested paths:
- `rv$tg_results$model`, `rv$tg_results$posthoc`, `rv$tg_results$data`, ...
- `rv$surv_results$fit`, `rv$surv_results$pairwise`, ...
- `rv$syn_results$bliss`, `rv$syn_results$loewe`, ...
- `rv$pw_results$plots`, `rv$pw_results$power_table`, ...
- `rv$dr_results$volume_stats`, `rv$dr_results$growth_stats`, ...

With proper Shiny modules, each module would manage its own reactives internally and expose only necessary values via return values.

### 3.3 JavaScript Injection Pattern

app.R injects ~50 lines of inline JavaScript via `tags$script()` for:
- Restructuring the navbar (moving disclaimer link)
- Disabling/greying out the diagnostics tab
- Tab switching handlers

This works but is fragile — it depends on specific Shiny HTML structure (CSS class names, element ordering). Consider using `shinyjs` functions or dedicated JS files in `www/js/` for maintainability.

### 3.4 No Input Validation Framework

Each module independently validates uploaded data with ad-hoc checks. There's no consistent validation framework that:
- Checks required columns exist
- Validates data types (numeric volumes, date formats)
- Provides user-friendly error messages via `validate(need(...))` consistently
- Handles edge cases (empty files, single-group data)

helpers.R has `validate_uploaded_file()` and `detect_columns()`, but modules don't all use them consistently.

### 3.5 Error Handling Strategy

Modules use a mix of:
- `tryCatch()` with `showNotification()`
- `validate(need(...))` in render functions
- `req()` for NULL checks
- `try()` with silent failure

There's no consistent pattern for which strategy to use where, leading to some errors being silently swallowed and others shown as full stack traces.

---

## 4. Enhancement Opportunities

### 4.1 Additional Functionality

| Feature | Description | Priority |
|---------|-------------|----------|
| **Body weight tab** | Tab exists but is disabled — implement toxicity/weight tracking | High |
| **Data export (CSV)** | Export analysis results as CSV, not just plots | High |
| **Report generation** | R Markdown report combining all analyses | Medium |
| **Session save/restore** | Save analysis state and reload later | Medium |
| **Batch analysis** | Upload multiple datasets and compare | Low |
| **Interactive plots** | plotly integration for hover/zoom on key plots | Medium |
| **Undo/redo** | Allow undoing parameter changes | Low |
| **Analysis log** | Show what analyses were run and with what parameters | Medium |

### 4.2 UI/UX Improvements

1. **Loading indicators:** Long analyses (LME4, power) have no progress indicator. Add `withProgress()` or `waiter` package spinners.

2. **Input persistence:** When switching tabs, re-uploading data clears all settings. Consider caching last-used settings.

3. **Help tooltips:** Add `bsTooltip()` or inline help for statistical parameters (e.g., "What is Bonferroni correction?").

4. **Responsive design:** Current fixed-width layout doesn't adapt well to different screen sizes. Consider `fluidPage()` responsive breakpoints.

5. **Color accessibility:** No built-in colorblind-friendly palette option. Add a "Colorblind safe" preset with viridis or similar palette.

6. **Validation feedback:** Show green checkmarks when data passes validation, red X when it doesn't, instead of only showing errors after analysis attempt.

### 4.3 Testing

No automated tests exist. Add:
1. **Unit tests** for helpers.R functions (pure functions, easy to test)
2. **Shiny testServer** tests for each module's reactive logic
3. **shinytest2** end-to-end tests with demo data uploads
4. **Snapshot tests** for plot outputs (vdiffr package)

---

## 5. File-by-File Summary

| File | Lines | Status | Key Issues |
|------|-------|--------|------------|
| server_modules.R | 2,396 | Warn | 1,850-line god function; D1–D3; 10+ debug cat(); hardcoded CI |
| survival_modules.R | 1,007 | Warn | Duplicated preprocessing; D6, D8; debug cat(); double-computation |
| ui_modules.R | 817 | OK | Well-structured with helpers; no NS() |
| synergy_modules.R | 683 | Bug | D4, D5, D9; stale downloads; duplicated preprocessing |
| helpers.R | 644 | Good | Clean utilities; well-documented |
| dose_response_modules.R | 603 | Warn | Triple-colon access; deprecated aes_string(); redefined %||% |
| power_modules.R | 409 | Warn | D7 (hardcoded sig.level); stale download; limited config |
| app.R | 260 | OK | Clean entry point; inline JS could be external |
| global.R | 62 | OK | Redundant library() calls; double useShinyjs() |
| zzz.R | 19 | Good | Simple .onLoad hook |

---

## 6. Priority Action Items

### Critical (Fix Before Release)
1. **D4:** Fix synergy download handlers to save current (customized) plots
2. **D7:** Use user-selected alpha in power curve rendering
3. **D2:** Use `qt()` with proper degrees of freedom for forest plot CIs

### High Priority
4. Consolidate date/volume preprocessing into a single shared function
5. Remove all debug `cat()` statements
6. Fix download-render mismatches across all modules (D4, power, survival, dose-response)
7. **D5:** Either set `rv$syn_results$plot_trend` or remove the download button
8. **D6:** Add `surv_palette` UI input or remove reference

### Medium Priority
9. Refactor `tumor_growth_server()` into Shiny modules with NS()
10. Replace deprecated `aes_string()` in dose_response_modules.R
11. Export `prepare_dose_data()` from mouseExperiment (remove triple-colon)
12. Add loading indicators for long-running analyses
13. Add unit tests for helpers.R
14. Move inline JS to external file in www/js/

### Low Priority
15. Remove redundant `%||%` definition
16. Remove redundant `useShinyjs()` call
17. Add colorblind-friendly palette option
18. Add interactive (plotly) plots for key visualizations
19. Add R Markdown report generation
20. Implement body weight analysis tab

---

## 7. Cross-Package Issues

Issues that span both mouseExperiment and mouseExperimentDashboard:

1. **AUC return structure inconsistency:** The package returns AUC results in a different format than LME4 results, forcing the dashboard to have if/else branches in comparisons, forest plots, etc. Standardize at the package level.

2. **No progress callbacks:** Long-running package functions (power analysis, large dataset LME4) block the Shiny UI with no feedback. Consider adding `progress` parameter support.

3. **Growth rate display:** The double-log bug (B2 in mouseExperiment) means the dashboard displays incorrect growth rates. Fixing the package bug will automatically fix the dashboard display.

4. **Download fidelity:** The dashboard builds its own ggplot layers on top of package-generated plots (custom colors, shapes, legend positions). Downloads should capture these additions, not the bare package output.

---

*Review based on complete source code examination of all dashboard files and cross-reference with mouseExperiment package findings.*
