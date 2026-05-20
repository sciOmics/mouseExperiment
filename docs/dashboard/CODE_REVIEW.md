# mouseExperimentDashboard Code Review

**Reviewer:** Claude Code
**Date:** 2026-05-20
**Version reviewed:** v26.05.017
**Scope:** All Shiny module R files, server logic, and UI definitions

---

## Summary

The dashboard covers the full analysis surface of the `mouseExperiment` package and is broadly well-structured around Shiny modules. The main themes are: (1) a handful of real bugs — a hardcoded cage fallback, silent plot download failures, and a missing volume-column guard; (2) UI/UX inconsistencies between otherwise similar analysis modules; (3) some structural improvements that would make all Bayesian modules consistent with each other; and (4) missing features in Dose-Response and Power Analysis relative to the other modules.

Issues are categorised **Bug** (wrong behaviour or crashes), **Inconsistency** (same thing done differently in comparable modules), **Performance**, **Missing**, and **Enhancement** (a worthwhile addition).

---

## 1. Bugs and Logic Errors

### 1.1 Hardcoded cage column fallback in Toxicity module
**File:** `toxicity_modules.R` (Mixed Model analysis call, ~line 596)

```r
cage_column = if (!is.null(rv$settings$col_cage) && nzchar(rv$settings$col_cage))
                rv$settings$col_cage
              else "Cage"   # ← hardcoded fallback
```

When the user has not mapped a cage column, the code passes `cage_column = "Cage"` instead of `NULL`. `analyze_body_weight()` (and the Bayesian equivalent) will then try to read a column named "Cage" that does not exist in the uploaded data. This causes a downstream error rather than simply omitting the cage effect.

**Fix:** Change the `else` branch to `else NULL`.

---

### 1.2 Plot download handlers do not validate the plot object
**Files:** `survival_modules.R:1023-1103`, `synergy_modules.R:779-948`, `toxicity_modules.R` (multiple download handlers)

Download handlers call `req(result_object)` to gate on analysis completion, but do not verify that the specific plot within the result is non-NULL. If `plots = FALSE` was used or a plot-rendering step threw a silent error, the handler will call `ggsave()` (or equivalent) with `plot = NULL`, producing an empty or corrupt PNG file with no user-facing message.

**Fix:**
```r
downloadHandler(
  filename = ...,
  content = function(file) {
    p <- rv$surv_results$km_plot
    if (is.null(p)) {
      showNotification("Plot is not available.", type = "warning")
      return(invisible(NULL))
    }
    ggsave(file, plot = p, ...)
  }
)
```

---

### 1.3 Event-from-volume branch does not guard against missing volume column
**File:** `survival_modules.R` (~line 556)

In the `event_src == "volume"` code path, the code reads `full_data[[vol_col_surv]]` where `vol_col_surv` is derived from `rv$settings$col_volume`. If volume preprocessing failed silently or the user hasn't mapped a volume column, `vol_col_surv` is NULL and the `[[NULL]]` access throws a cryptic error inside `safe_analysis()`.

**Fix:** Add an early check before entering the volume-threshold path:
```r
if (is.null(rv$settings$col_volume) || !nzchar(rv$settings$col_volume)) {
  showNotification(
    "Volume column not mapped. Cannot derive survival events from volume threshold.",
    type = "error"
  )
  return(NULL)
}
```

---

### 1.4 Date format detection failure is not caught before analysis
**File:** `helpers.R` (~line 44-76, `detect_date_format()`)

`detect_date_format()` returns `NULL` when no pattern matches. `calculate_dates()` is then called with `date_format = NULL`, which silently guesses the format. If the guess is wrong, all computed `Day` values are incorrect and no error or warning is surfaced to the user. Analysis runs to completion on corrupted time data.

**Fix:** After calling `detect_date_format()`, check the result:
```r
fmt <- detect_date_format(rv$uploaded_data[[rv$settings$col_date]])
if (is.null(fmt)) {
  showNotification(
    "Could not auto-detect date format. Please verify your date column.",
    type = "warning"
  )
}
```

---

### 1.5 `ensure_preprocessed()` result ignored in survival loop
**File:** `survival_modules.R` (~line 461 vs 545)

`ensure_preprocessed()` is called and its result stored in `pp`, but the subsequent loop that computes per-group statistics still reads from `rv$uploaded_data` in places — using the raw (unpreprocessed) data rather than `pp$processed`. Any preprocessing transformation (volume calculation, day calculation) is bypassed in that path.

**Fix:** Audit all reads inside the survival analysis loop to ensure `pp$processed` is used consistently.

---

## 2. UI/UX Inconsistencies

### 2.1 Tab names differ across analysis modules for the same concept
**Files:** `survival_modules.R`, `toxicity_modules.R`, `synergy_modules.R`, `dose_response_modules.R`

| Module | Model output tab | Pairwise comparison tab |
|--------|-----------------|------------------------|
| Survival | "Cox Model" | "Pairwise" |
| Toxicity | "Mixed Model" | "Pairwise Comparisons" |
| Synergy | "Synergy Metrics" | *(none)* |
| Dose-Response | "Statistics" | "Comparisons" |

Users switching between modules must re-learn tab names each time. The result summary should always be called "Summary"; the model output "Model Results"; pairwise comparisons "Pairwise Comparisons".

---

### 2.2 Download button placement is not consistent
**Files:** All modules

- Survival: download buttons inline in tab headers
- Toxicity: download buttons at the bottom of some tabs, in the header of others
- Synergy: download buttons appear only after a plot selector
- Dose-Response: download buttons only for tables (no plot export — see §4.1)

A user who learns where to find download buttons in Toxicity will not find them in the same location in Synergy. All download buttons should be at the top of the relevant tab panel, to the right of the tab title.

---

### 2.3 Reference group selector label is inconsistent
**Files:** `survival_modules.R:358`, `toxicity_modules.R:710`, `synergy_modules.R:303`

| Module | Label |
|--------|-------|
| Survival | "Reference Group" |
| Toxicity | "Control / Reference Group" |
| Synergy | "Control Group" |

Standardise to "Reference Group" across all modules.

---

### 2.4 Plot customization controls vary substantially between modules
**Files:** Plot sidebar sections in all module files

- **Survival:** confidence intervals toggle, risk table toggle, X-axis start, legend position, title, subtitle, text size, plot dimensions
- **Toxicity:** trajectory type, center bar, error bars, point size, jitter, group order/colour, text size, line width, plot dimensions
- **Synergy:** colours, text size, legend position, shading toggle, CI region, title, subtitle, plot dimensions
- **Dose-Response:** minimal or no customization controls

Text size and plot dimensions should be universally available. Group colour pickers and legend position should be available wherever groups are plotted. Confidence interval / error bar toggles should be available for all quantitative plots.

---

### 2.5 Progress feedback granularity differs across modules
**Files:** All modules' `observeEvent(input$run_*)` handlers

- Toxicity: six progress steps (0.05, 0.2, 0.5, 0.75, 0.85, 1.0)
- Survival: two progress steps (0.2, 0.9)
- Synergy: one progress step (0.2, then implicitly 1.0)

Survival and Synergy feel "stuck" at 20% for most of the run. All modules should have at least three progress checkpoints: data validation (~10%), model fitting (~60%), post-processing/plots (~95%).

---

## 3. Performance

### 3.1 Data preprocessing is potentially re-run in every module
**File:** `helpers.R` (`ensure_preprocessed()`), all module files

`ensure_preprocessed()` is designed as a cache check, but the implementation re-runs `calculate_volume()` and `calculate_dates()` if `rv$processed_data` is NULL. Since `rv$processed_data` is reset on new data upload, the first call from each module re-runs preprocessing from scratch. For large datasets this is a multi-second stall in the middle of analysis.

**Fix:** Cache the preprocessed data once at upload time (in the upload `observeEvent`) rather than lazily on first analysis call.

---

### 3.2 `get_ordered_groups()` parses input IDs on every render
**Files:** `survival_modules.R:762`, `toxicity_modules.R:1194`, `synergy_modules.R:872`

`get_ordered_groups()` iterates over all Shiny input IDs to reconstruct the current group order. This runs on every `renderPlot()` call, including minor reactive invalidations unrelated to group order. Wrap the result in `reactive()` and memoize it.

---

### 3.3 Full model objects stored in session reactive values
**Files:** `survival_modules.R:321`, `toxicity_modules.R:539`, `synergy_modules.R:267`

`rv$surv_results`, `rv$tox_results`, `rv$syn_results` store the complete analysis result including `brmsfit` objects (when Bayesian mode is used). A single `brmsfit` can be 50–200 MB for realistic models. Storing multiple of these in session memory limits concurrent users and can cause OOM crashes.

**Fix:** Store only the display-ready artefacts (data frames and plot objects) in reactive values. If the `brmsfit` is needed (e.g., for `bayesian_therapeutic_window`), store it in a separate `rv$bayes_models` list and use `gc()` after extracting the required summaries.

---

## 4. Missing Features

### 4.1 Dose-Response module has no plot export
**File:** `dose_response_modules.R`

Every other analysis module exposes PNG download handlers for its plots (KM curve, forest plot, trajectory, synergy bar chart). The Dose-Response module provides CSV downloads only; the dose-response curve, the model fit plot, and the statistics plot have no download button. This prevents publishing the dose-response figure without leaving the app.

**Fix:** Add `downloadButton()` / `downloadHandler()` for `dose_response_curve_plot` and any other plots returned by `bayesian_dose_response()` or `dose_response_statistics()`.

---

### 4.2 Comparison type selector absent from Dose-Response
**File:** `dose_response_modules.R`

Toxicity has "vs. Reference Group" / "All Pairwise" selector for comparisons. Dose-Response has no equivalent — ANOVA always produces Tukey HSD pairwise results. Adding a "vs. Control only" option would reduce output clutter for the common case.

---

### 4.3 Bayesian MCMC diagnostic plots not exposed for all Bayesian modules
**Files:** Survival, Toxicity (Bayesian paths), Power Analysis

Tumor Growth Bayesian mode surfaces `pp_check_plot` and `mcmc_trace_plot`. The Survival Bayesian mode (`bayesian_survival`) returns these same plots in the result list but they are not surfaced in a dedicated tab. Power Analysis (LMM) has no diagnostic output at all.

**Fix:** Add "Posterior Check" and "MCMC Diagnostics" tabs to Survival and Toxicity Bayesian result panels, populated from `res$pp_check_plot` and `res$mcmc_trace_plot` respectively.

---

### 4.4 Prior vs Posterior plot not exposed for Bayesian Body Weight
**File:** `toxicity_modules.R` (Bayesian BW section)

`bayesian_body_weight()` returns `prior_posterior_plot` (same as `bayesian_tumor_growth`), but the Toxicity module does not display it. This means users cannot inspect prior sensitivity for body weight models.

**Fix:** Add "Prior vs Posterior" tab to the Bayesian Body Weight result panel.

---

### 4.5 No UI for `bayesian_therapeutic_window()`
**File:** `toxicity_modules.R`

`bayesian_therapeutic_window()` (v0.4.1) is implemented in the package but not yet exposed in the dashboard. The function requires both a TG and BW Bayesian result to be available — the dashboard should detect when both are present and enable a "Bayesian TWM" run button.

---

### 4.6 No UI for `bayesian_synergy()`
**File:** `synergy_modules.R`

`bayesian_synergy()` (v0.4.2) is implemented in the package but not yet exposed in the dashboard. The synergy module should offer a "Bayesian (brms)" model type alongside the existing frequentist mode, following the same pattern as Tumor Growth.

---

### 4.7 No UI for `bayesian_dose_response()`
**File:** `dose_response_modules.R`

`bayesian_dose_response()` (v0.4.0) is implemented in the package but not yet exposed. The dose-response module should offer a "Bayesian Hill/Emax (brms)" model type.

---

## 5. Error Handling

### 5.1 Raw R error messages surface to users
**Files:** All modules' `safe_analysis()` wrappers

`safe_analysis()` catches errors and stores `conditionMessage(e)` directly in `result$message`. These raw messages (e.g., `"object 'TreatmentGroup' not found"`, `"subscript out of bounds"`) are displayed in the UI without translation. Non-technical users cannot act on these messages.

**Fix:** Build a table of known error patterns → user-friendly messages:
```r
translate_error <- function(msg) {
  if (grepl("not found in", msg, fixed = TRUE))
    return("A required column was not found. Check your column mappings in Data Upload.")
  if (grepl("singular", msg, fixed = TRUE))
    return("Model fitting failed (singular fit). Try increasing sample size or simplifying the random effects structure.")
  msg  # fallback: show original
}
```

---

### 5.2 Silent failures in download handlers
**Files:** Multiple download handlers (see §1.2)

When `req()` returns FALSE (required reactive is NULL), the download handler silently cancels. The user sees the download button disappear briefly with no feedback. This is indistinguishable from a network error.

**Fix:** Use `shinyjs::disable()` on download buttons when their data is not ready, with a tooltip explaining what needs to run first.

---

### 5.3 Analytics JS escaping is incomplete
**Files:** `survival_modules.R:668`, `synergy_modules.R:473`

Error messages are escaped for JavaScript with `gsub("'", "\\\\'", result$message)` but backticks, newlines, and Unicode characters are not escaped. A multi-line R traceback in `result$message` will break the JS call and swallow the analytics event.

**Fix:** Use `jsonlite::toJSON(result$message, auto_unbox = TRUE)` to produce a properly escaped JSON string, then embed it in the JS call.

---

## 6. Code Structure

### 6.1 Server functions are too long for effective maintenance
**Files:** `survival_modules.R:318-1164` (~850 lines), `toxicity_modules.R:536-1341+` (~800 lines), `synergy_modules.R:264-950` (~690 lines)

A single `moduleServer()` function of 800+ lines is difficult to navigate, test, or refactor. Typical Shiny best practice targets <200 lines per server function by extracting plot-building, table-building, and download logic into named functions.

**Recommendation:** Extract each plot rendering block into a standalone function (e.g., `render_km_plot(results, settings)`) that takes reactive values as arguments. This makes each piece independently testable.

---

### 6.2 Plot-building code is duplicated between render and download
**Files:** `survival_modules.R`, `toxicity_modules.R`

The plot rendering for `renderPlot()` and the corresponding `downloadHandler()` contain near-identical code, differing only in output target. A code change to the KM plot appearance must be made in two places.

**Fix:** Compute the plot once in a reactive, then reference it in both `renderPlot()` and `downloadHandler()`:
```r
km_plot_r <- reactive({
  req(rv$surv_results)
  build_km_plot(rv$surv_results, input)
})
output$plot_km <- renderPlot(km_plot_r())
output$download_km <- downloadHandler(
  content = function(file) ggsave(file, plot = km_plot_r(), ...)
)
```

---

### 6.3 Synergy module bypasses `build_color_vector()` in some code paths
**File:** `synergy_modules.R:619-625`

Most modules use the shared `build_color_vector()` helper for consistent colour assignment. The synergy module has an inline colour construction block in at least one plot path that does not call `build_color_vector()`. This means user-selected colours in the colour picker are not applied to that plot.

**Fix:** Audit all `ggplot` calls in `synergy_modules.R` and replace inline colour construction with `build_color_vector(groups, input)`.

---

## 7. Harmonization — Making All Bayesian Modules Consistent

### 7.1 Pattern for Bayesian analysis is not yet implemented for Synergy, Dose-Response, and TWM
The Tumor Growth module established the pattern for Bayesian analysis integration: model type selector → async brms fit → result tabs including Summary, Posterior Summary, MCMC Diagnostics, Posterior Check, MCMC Trace, Prior vs Posterior. Survival and Body Weight partially follow this pattern. Synergy, Dose-Response, and TWM have no Bayesian path at all yet (per §4.5–4.7).

**Template to follow:** For each new Bayesian module:
1. Add `"Bayesian (brms)"` to the model type `selectInput` 
2. Add prior-strength, chains, and iterations sidebar controls under `conditionalPanel(condition = "input.model_type == 'bayesian'")`
3. Run the fit asynchronously using `future_promise()`
4. Display results in standardised tabs: **Summary** (treatment effects table), **Posterior Summary** (posterior_summary DF), **MCMC Diagnostics** (mcmc_diagnostics DF with Rhat / ESS / Converged), **Model Plots** (pp_check, trace), **Prior vs Posterior** (prior_posterior_plot)

---

### 7.2 MCMC Diagnostics tab should display identically in all Bayesian modules
**Files:** `mod_tumor_growth.R`, `survival_modules.R`, `toxicity_modules.R`

Tumor Growth's MCMC Diagnostics tab shows: a DT table of Rhat/ESS/Converged, a colour rule highlighting Rhat > 1.01 in red, and the trace plot below. Survival shows only the table. Body Weight shows only the table. They should all show the table + trace plot in the same layout.

---

### 7.3 Download filenames do not include study context
**Files:** All `downloadHandler()` calls

Files are named generically: `"kaplan_meier_plot.png"`, `"treatment_effects.csv"`. When users download from multiple analyses in a session, files overwrite each other or require manual renaming. A date stamp and reference group name would make files self-documenting.

**Fix:**
```r
filename = function() {
  ref <- input$surv_reference_group
  paste0("km_plot_", ref, "_", format(Sys.Date(), "%Y%m%d"), ".png")
}
```

---

## 8. Quick-Reference Checklist

| # | Issue | Type | File | Status |
|---|-------|------|------|--------|
| 1.1 | Hardcoded `"Cage"` fallback instead of `NULL` | Bug | `toxicity_modules.R:596` | Open |
| 1.2 | Plot download handlers do not validate plot object | Bug | Multiple | Open |
| 1.3 | Volume column not guarded in event-from-volume path | Bug | `survival_modules.R:556` | Open |
| 1.4 | `detect_date_format()` NULL return not caught | Bug | `helpers.R:44` | Open |
| 1.5 | Survival loop reads raw data instead of `pp$processed` | Bug | `survival_modules.R:545` | Open |
| 2.1 | Tab names differ across modules for the same concept | Inconsistency | Multiple | Open |
| 2.2 | Download button placement not consistent | Inconsistency | Multiple | Open |
| 2.3 | Reference group label: "Reference Group" vs "Control Group" vs "Control / Reference Group" | Inconsistency | Multiple | Open |
| 2.4 | Plot customization controls differ substantially between modules | Inconsistency | Multiple | Open |
| 2.5 | Progress feedback: Toxicity has 6 steps, Synergy has 1 | Inconsistency | Multiple | Open |
| 3.1 | `ensure_preprocessed()` may re-run preprocessing for each module | Performance | `helpers.R` | Open |
| 3.2 | `get_ordered_groups()` runs on every render; not memoized | Performance | Multiple | Open |
| 3.3 | Full `brmsfit` objects stored in session reactive values | Performance | Multiple | Open |
| 4.1 | No plot export in Dose-Response module | Missing | `dose_response_modules.R` | Open |
| 4.2 | No comparison type selector in Dose-Response | Missing | `dose_response_modules.R` | Open |
| 4.3 | Bayesian MCMC diagnostic plots not in Survival/Toxicity Bayesian panels | Missing | Multiple | Open |
| 4.4 | Prior vs Posterior plot not shown for Bayesian Body Weight | Missing | `toxicity_modules.R` | Open |
| 4.5 | No UI for `bayesian_therapeutic_window()` | Missing | `toxicity_modules.R` | Open |
| 4.6 | No UI for `bayesian_synergy()` | Missing | `synergy_modules.R` | Open |
| 4.7 | No UI for `bayesian_dose_response()` | Missing | `dose_response_modules.R` | Open |
| 5.1 | Raw R error messages surface to users | Error handling | Multiple | Open |
| 5.2 | Silent failures in download handlers (no user feedback) | Error handling | Multiple | Open |
| 5.3 | Analytics JS escaping incomplete for multiline/unicode errors | Error handling | Multiple | Open |
| 6.1 | Server functions 800+ lines; hard to maintain | Structure | Multiple | Open |
| 6.2 | Plot code duplicated in renderPlot + downloadHandler | Structure | Multiple | Open |
| 6.3 | Synergy inline colour construction bypasses `build_color_vector()` | Structure / Bug | `synergy_modules.R:619` | Open |
| 7.1 | Bayesian analysis pattern not yet implemented for Synergy / DR / TWM | Harmonization | Multiple | Open |
| 7.2 | MCMC Diagnostics tab layout differs: TG shows trace plot, others don't | Harmonization | Multiple | Open |
| 7.3 | Download filenames are generic; no date/reference context | Enhancement | Multiple | Open |
