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

### 2.1 God Function: `tumor_growth_server()` — ✅ COMPLETED

Originally, server_modules.R contained `tumor_growth_server()` at ~1,850 lines with all TG concerns mixed together.

**Status:** Refactored via commit `51a4c58`
- Split into `mod_data_upload.R` (530 lines, data upload/validation/settings sync)
- Split into `mod_tumor_growth.R` (1,838 lines, analysis/renders/downloads)
- Both now use `moduleServer()` pattern with proper namespacing
- Each module manages its own inputs/outputs within its namespace
- rv$settings bridge enables clean inter-module communication for column mapping

---

### 2.2 No Shiny Module Namespacing — ✅ COMPLETED

Originally, the app used flat input/output IDs with no `NS()` namespacing, creating collision risk and refactoring brittleness.

**Status:** Fully modularized via commit `51a4c58`
- All 6 server functions converted to `moduleServer()` pattern
- All 6 UI functions now use `NS(id)` wrapping
- All static input IDs wrapped with `ns()` in UI functions
- All dynamic renderUI inputs wrapped with `session$ns()` in server functions
- All conditionalPanel calls include `ns = ns` parameter
- Helper functions (`group_order_and_color_ui`, `column_mapping_ui`, etc.) updated to accept `ns` parameter
- Result: All IDs properly namespaced (e.g., `data-data_file`, `tg-tg_model`, `surv-surv_fit`)

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

### 3.1 File Organization — ✅ COMPLETED

Originally, structure mixed concerns with server_modules.R containing a 1,850-line god function and no clear module boundaries.

**Current structure (post-commit `51a4c58`):**

| File | Lines | Responsibilities |
|------|-------|-----------------|
| app.R | 189 | Main app shell with module wiring and rv observers |
| mod_data_upload.R | 530 | Data upload, column mapping, preview, rv$settings sync |
| mod_tumor_growth.R | 1,838 | TG analysis dispatch, renders (25+), downloads (10+) |
| ui_modules.R | 817 | UI definitions for all tabs (with NS wrapping) |
| survival_modules.R | 921 | Survival analysis server (moduleServer pattern) |
| synergy_modules.R | 708 | Synergy analysis server (moduleServer pattern) |
| dose_response_modules.R | 607 | DR analysis server (moduleServer pattern) |
| power_modules.R | 427 | Power analysis server (moduleServer pattern) |
| helpers.R | 754 | Shared utilities (updated with ns parameters) |
| zzz.R | 19 | `.onLoad` hook |

**Completed:**
- ✅ server_modules.R split into mod_data_upload.R + mod_tumor_growth.R
- ✅ All modules follow moduleServer()/NS() pattern
- ✅ Clear separation of concerns per file
- ⏳ Future: Consolidate date/volume preprocessing (still duplicated 4×, deferred to next phase)

### 3.2 Reactive Value Sprawl — ✅ PARTIALLY COMPLETED

Originally, modules stored results in shared `rv` with deeply nested paths: `rv$tg_results$model`, `rv$surv_results$fit`, etc.

**Status:** Improved via commit `51a4c58`
- ✅ Each module now uses `moduleServer()` with internal namespacing
- ✅ Reduced cross-module input collision via NS() — each module's inputs isolated
- ✅ Introduced `rv$settings` bridge pattern for column-mapping shared state (14 column inputs synced once in data_upload_server, read by all downstream modules)
- ✅ Each module still maintains results in `rv$*_results`, but this is appropriate for caching across re-renders
- ⏳ Future: Further isolation could expose modules' results via return values instead of rv paths (requires changes to app.R observer structure)

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
| app.R | 189 | Good | Module wiring complete; inline JS still present (3.3) |
| mod_data_upload.R | 530 | Good | Data upload, column mapping, rv$settings sync ✅ |
| mod_tumor_growth.R | 1,838 | Good | TG analysis, renders, downloads ✅; D1, D2, D3 remain |
| ui_modules.R | 817 | Good | Well-structured with NS() wrapping; modularized ✅ |
| survival_modules.R | 921 | Warn | D6, D8; duplicated preprocessing; debug cat(); double-computation |
| synergy_modules.R | 708 | Bug | D4, D5, D9; stale downloads; duplicated preprocessing |
| dose_response_modules.R | 607 | Warn | Triple-colon access; deprecated aes_string(); redefined %||% |
| power_modules.R | 427 | Warn | D7 (hardcoded sig.level); stale download; escaped quote fix ✅ |
| helpers.R | 754 | Good | Clean utilities; NS parameters added ✅ |
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
9. ✅ **COMPLETED:** Refactor `tumor_growth_server()` into Shiny modules with NS() (commit `51a4c58`)
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
