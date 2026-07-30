# mouseExperimentDashboard Code Review

**Reviewer:** Claude Code
**Date:** 2026-05-20
**Version reviewed:** v26.05.019
**Scope:** All Shiny module R files, server logic, and UI definitions

**Resolution status:** ✅ All 25 items from the v26.05.019 review are resolved. The last open item (**2.4** Dose-Response plot customization) was closed in v26.06.005 by adding per-dose centre / error overlays + a jitter-width slider to the DR Plots tab.

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

**✅ Resolved across v26.06.002 + v26.06.005.** v26.06.002 added the
universal controls (legend position, axis / title / legend label text
inputs). v26.06.005 closed the remaining gap with per-dose centre
overlay (mean / median), error-bar overlay (SEM / SD / 95 % CI), and a
jitter-width slider on both the Tumor Volume vs Dose scatter and the
TGI subplot.

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
| 1.1 | Hardcoded `"Cage"` fallback instead of `NULL` | Bug | `toxicity_modules.R:596` | Done |
| 1.2 | Plot download handlers do not validate plot object | Bug | Multiple | Done |
| 1.3 | Volume column not guarded in event-from-volume path | Bug | `survival_modules.R:556` | Done |
| 1.4 | `detect_date_format()` NULL return not caught | Bug | `helpers.R:44` | Done |
| 1.5 | Survival loop reads raw data instead of `pp$processed` | Bug | `survival_modules.R:545` | Done |
| 2.1 | Tab names differ across modules for the same concept | Inconsistency | Multiple | Done |
| 2.2 | Download button placement not consistent | Inconsistency | Multiple | Done |
| 2.3 | Reference group label: "Reference Group" vs "Control Group" vs "Control / Reference Group" | Inconsistency | Multiple | Done |
| 2.4 | Plot customization controls differ substantially between modules | Inconsistency | Multiple | ✅ Fixed v26.06.002 (DR caught up to other modules) |
| 2.5 | Progress feedback: Toxicity has 6 steps, Synergy has 1 | Inconsistency | Multiple | Done |
| 3.1 | `ensure_preprocessed()` may re-run preprocessing for each module | Performance | `helpers.R` | Done |
| 3.2 | `get_ordered_groups()` runs on every render; not memoized | Performance | Multiple | Done |
| 3.3 | Full `brmsfit` objects stored in session reactive values | Performance | Multiple | Done |
| 4.1 | No plot export in Dose-Response module | Missing | `dose_response_modules.R` | Done |
| 4.2 | No comparison type selector in Dose-Response | Missing | `dose_response_modules.R` | Done |
| 4.3 | Bayesian MCMC diagnostic plots not in Survival/Toxicity Bayesian panels | Missing | Multiple | Done |
| 4.4 | Prior vs Posterior plot not shown for Bayesian Body Weight | Missing | `toxicity_modules.R` | Done |
| 4.5 | No UI for `bayesian_therapeutic_window()` | Missing | `toxicity_modules.R` | Done |
| 4.6 | No UI for `bayesian_synergy()` | Missing | `synergy_modules.R` | Done |
| 4.7 | No UI for `bayesian_dose_response()` | Missing | `dose_response_modules.R` | Done |
| 5.1 | Raw R error messages surface to users | Error handling | Multiple | Done |
| 5.2 | Silent failures in download handlers (no user feedback) | Error handling | Multiple | Done |
| 5.3 | Analytics JS escaping incomplete for multiline/unicode errors | Error handling | Multiple | Done |
| 6.1 | Server functions 800+ lines; hard to maintain | Structure | Multiple | Done (v26.06.003 + v26.06.004) — all 5 mega-files split by concern; 7 new helper files |
| 6.2 | Plot code duplicated in renderPlot + downloadHandler | Structure | Multiple | Done |
| 6.3 | Synergy inline colour construction bypasses `build_color_vector()` | Structure / Bug | `synergy_modules.R:619` | Done |
| 7.1 | Bayesian analysis pattern not yet implemented for Synergy / DR / TWM | Harmonization | Multiple | Done |
| 7.2 | MCMC Diagnostics tab layout differs: TG shows trace plot, others don't | Harmonization | Multiple | Done |
| 7.3 | Download filenames are generic; no date/reference context | Enhancement | Multiple | Done |

---

# Review Round 2 (v26.05.022 — 2026-05-30)

**Reviewer:** Claude Code
**Scope:** Items introduced or remaining after the v26.05.019 review, focusing on the modal/progress migration, bugs surfaced during it, deployment correctness, and statistical-diagnostic surface area inherited from the backend.

## Round-1 status

Of the 25 items in the v26.05.019 review, **all 25 are resolved**. Item **6.1 (server function length)** was closed in v26.06.004 — all five mega-files have been split by concern (see C.1 below). Item **2.4 (Dose-Response plot customization)** was closed in v26.06.005 by adding per-dose centre / error overlays and a jitter-width slider that bring the DR Plots tab in line with the Tumor Growth and Toxicity plot tabs.

The Bayesian module pattern (item 7.1) has been fully harmonised: all five analysis modules (Tumor Growth, Survival, Synergy, Dose Response, Toxicity TWM/BW) now share a common `show_analysis_modal()` running indicator, `Posterior Check` / `MCMC Diagnostics` / `Prior vs Posterior` tab structure, and posterior-summary table layout.

---

## A. Bugs found and fixed during this review cycle

The v26.05.020–022 progression surfaced four interacting issues. All four are now fixed but are documented here because they collectively explain a long stretch of confusing behaviour.

### A.1 Orphaned `progress$set()` references after Progress → modal migration — **FIXED (56bfaf5)**
**Files:** all `*_modules.R`, `mod_tumor_growth.R`

When `Progress$new()` was replaced with `show_analysis_modal()` in v26.05.020, the `progress <-` definition was removed but 13 trailing `progress$set(...)` calls were left scattered through the analysis observers. The TG-LME4 / TG-AUC versions sat inside `safe_analysis()` and were translated by `helpers.R::translate_error()` (which maps `"object .* not found"` to the user-facing "A required column was not found") — making the bug look like a column-mapping problem. The Synergy/Survival/Power/Dose-Response/Toxicity versions sat *after* `safe_analysis()`, so when reached they crashed the observer; Shiny does not flush the queued `removeModal()` on an observer error, so the modal stayed visible while the analysis silently never produced output.

**Lesson:** `translate_error()` should distinguish "missing column in data" from "undefined R object" — the regex `"object .* not found"` is too broad.

### A.2 Toxicity Bayesian BW modal torn down by stray `on.exit()` — **FIXED (56bfaf5)**
**File:** `toxicity_modules.R`

The frequentist `on.exit(shiny::removeModal())` was registered before the BW Bayesian branch (lines 1090 / 1095). When that branch showed its own `show_analysis_modal("Bayesian body weight fit running…")` and `return()`-ed to start its future, the still-pending `on.exit` fired and removed the just-shown modal. The future ran headless and only displayed once it finished. Fixed by calling `on.exit()` (clear) inside the Bayesian branch.

### A.3 Persistent "overlay on top of everything" was Shiny's `#shiny-disconnected-overlay` — **NOT a UI bug**
**Root cause:** NPM `proxy_read_timeout` (60 s default) was killing the websocket during long-running Bayesian fits, after which the Shiny client renders its built-in `#shiny-disconnected-overlay` (z-index 99998, position fixed, grey 50% opacity). This survived every CSS/JS backdrop fix because it's emitted by Shiny itself, only after the connection is already dead, and only for analyses long enough to exceed the timeout (Bayesian only — frequentist fits finished under 60 s).

The dashboard-side cache-busting fix in `5e63303` is necessary but **not sufficient**; the root cure is **applying the 1 h `proxy_read_timeout` to the live NPM proxy hosts**. See VPS-repo commit `043419e` and the live-update curl snippet in this conversation. Until the live NPM config is updated, Bayesian fits on staging will continue to drop their connection.

### A.4 Stale browser cache of CSS/JS — **FIXED (5e63303)**
Every prior fix to `overlay.js` and `custom.css` was deployed correctly but never reached clients running the dashboard, because browsers were holding pre-fix cached assets indefinitely. Cache-busting via a `?v=<package-version>` query string on every local CSS/JS include forces refetch on each release.

---

## B. New deployment-correctness items

These are net-new to this review and live partly in this repo, partly in `mouseExperimentVPS`. The dashboard side is documented here so it's recoverable from this file alone.

### B.1 Silent rebuild failures
**File (companion repo):** `scripts/webhook-handlers/rebuild-shiny.sh`

`remotes::install_github()` treats a missing CRAN dependency as a warning and exits 0; `docker build` reports success; the rebuild script declared "Build successful" and recreated the container with an image whose `mouseExperimentDashboard` package never installed → R failed at startup → 502 with no signal. **Fixed in VPS commit `d2034c5`** which inserts a verification step (`docker run --rm <image> R -e "library(mouseExperimentDashboard)"`) before swapping the container.

### B.2 NPM proxy host config drift
**File (companion repo):** `scripts/06-configure-npm.sh`

The canonical provisioner created the Shiny hosts with neither `allow_websocket_upgrade: true` nor any `advanced_config`. A from-scratch reprovision would have produced a non-functional dashboard. **Fixed in VPS commit `043419e`**. Live NPM must still be updated via the API snippet in this conversation.

---

## C. Architecture

### C.1 Server function length (6.1 from v26.05.019 — ✅ Resolved across v26.06.003 + v26.06.004)
| File | LOC at v26.05.022 | LOC at v26.06.004 |
|---|---|---|
| `mod_tumor_growth.R` | **2,367** | **1,727** (split into 3 files, v26.06.003) |
| `mod_toxicity.R` | **2,069** | **1,595** (split into 3 files, v26.06.004) |
| `mod_survival.R` | **1,405** | **1,467** (downloads extracted, v26.06.004) |
| `mod_synergy.R` | **1,205** | **1,171** (downloads extracted, v26.06.004) |
| `helpers.R` | **1,083** | **1,153** (bayes + report extracted, v26.06.004) |

Note: `mod_survival.R` / `helpers.R` line counts grew between the v26.05.022
review snapshot and v26.06.002 because of feature additions (Round 2 D / F
items). The v26.06.004 numbers show post-extraction state.

**Resolution pattern** — splitting concerns into helpers in the same package
namespace, rather than the originally-proposed per-tab subdirectory layout
(R/ is flattened on install, so subdirs aren't usable). The split is by
concern, not by Bayesian/frequentist:

| Original file | Helper file(s) created | What lives there |
|---|---|---|
| `mod_tumor_growth.R` | `_downloads.R`, `_plots.R` | 16 downloads + necrosis DT / 5 plot reactives + 7 renderPlots + diagnostics + group-order UI |
| `mod_toxicity.R` | `_downloads.R`, `_plots.R` | 17 downloads / 10 freq plots + diagnostics + group-order UI |
| `mod_survival.R` | `_downloads.R` | 6 downloads (2 Bayesian + 4 frequentist) |
| `mod_synergy.R` | `_downloads.R` | 6 downloads |
| `helpers.R` | `helpers_bayes.R`, `helpers_report.R` | Bayesian diagnostics panel + LOO influential obs / HTML report builder |

Plot reactives are passed explicitly to the downloads helper rather than
captured by closure, so the cross-file contract is visible at the call
site. No subdirectory layout used.

### C.2 Module file naming inconsistency
Only `mod_tumor_growth.R` uses the `mod_` prefix; the others (`survival_modules.R`, `synergy_modules.R`, …) use the `_modules.R` suffix. Pick one convention.

### C.3 `output$` rendering is concentrated in a few mega-files
`survival_modules.R`, `mod_tumor_growth.R`, and `mod_data_upload.R` each contain dozens of `output$X <- render*()` blocks. Extracting per-tab `*_render_helpers.R` files would let test cases hit them individually.

**✅ Resolved across v26.06.003 + v26.06.004** for the four module mega-files — see C.1 above. Download handlers were extracted in every case; plot renderers were additionally extracted for tumor-growth and toxicity (the modules with the largest plot blocks). Table renderers and Bayesian renderers remain in the main files where they're interleaved with lifecycle code — the marginal gain from further splitting was small. `mod_data_upload.R` is 656 LOC and doesn't need the same treatment.

---

## D. Inherited statistical-diagnostic surface area

The dashboard should surface diagnostic outputs from the backend once they're added (see the companion `mouseExperiment` review, Round 2, sections A–C). Required dashboard work for each:

| Backend addition | Dashboard work |
|---|---|
| `cox.zph()` PH check (backend A.1) | New "PH Assumption" tab in Survival results showing the test table + Schoenfeld plot |
| NUTS divergences / max_treedepth / E-BFMI (backend A.2) | Inline warning banner in MCMC Diagnostics tab when any of these fail |
| LOO / Pareto-k (backend A.4) | New "Influential Observations" panel listing mice with Pareto-k > 0.7 |
| `bayes_R2()` (backend C.1) | New row in the Summary table for Bayesian paths |
| Posterior P(effect > 0) (backend C.2) | Add `P(effect ≠ 0)` column to Bayesian treatment-effects tables |
| Concordance / C-index (backend E.2) | Add to Survival Summary |

These should land alongside the backend changes — no value in shipping diagnostics that can't be inspected.

---

## E. UX / behavior items

### E.1 Power Analysis has no Summary tab → excluded from auto-switch (v26.05.022)
The auto-switch-to-Summary behaviour added in `8ee551d` skips Power Analysis because its first tab is "Info" — every other module's first tab is "Summary". Either:
- Add a Summary tab to Power Analysis (preferred — consistency), or
- Pick the first existing tab (e.g. "Scenario Table") as Power's "completion landing"

### E.2 CHANGELOG gap: no `v26.05.020` entry
`DESCRIPTION` was bumped to 26.05.020 (commit `de78be3` — Log tab removal, Help rewrite, Umami repoint) without a CHANGELOG entry. The v26.05.022 review folded forward without back-filling 020; a one-paragraph 020 entry would make the changelog complete.

### E.3 No startup health check for NPM timeout
Given the disconnect-overlay surprise (A.3), it would be defensible for `runDashboard()` to log a startup warning if `Sys.getenv("SHINY_PORT") != ""` (running under shiny-server) and a known-good proxy doesn't expose long-timeout config. Hard to detect from inside the container; alternatively document the requirement prominently in the deploy README.

---

## F. Recommended Additional Functionality

### F.1 Reproducible analysis report export
Generate a parameterised Quarto/R Markdown report for any completed analysis (PDF + HTML). All inputs, the modal-side analysis settings, the result tables, plots, and the model diagnostic outputs would be captured. Useful for inclusion in manuscripts and as an audit artifact for regulatory submissions.

### F.2 Save / load analysis presets
A user who has carefully tuned priors, transforms, random-effects spec, and reference group for one dataset has no way to apply that to a re-run with new data. A "Save Settings" button writing a small JSON blob, plus an "Import Settings" upload, would shorten the iteration loop.

### F.3 Side-by-side frequentist vs Bayesian comparison
For any analysis with both paths (TG, Survival, Synergy, DR), run both, present the treatment-effects tables side-by-side and the credible-interval forest plot overlaid with the confidence-interval forest plot. Useful pedagogically and for cross-validation.

### F.4 Per-mouse drill-down
A modal or new tab showing a single mouse's volume trajectory, weight trajectory, censoring, with the model's per-animal predicted curve overlaid. Currently the closest is the caterpillar plot of random effects.

### F.5 Re-run with subset
Quick filter ("exclude this mouse / cage / treatment group") that re-runs the analysis. Currently requires editing the source data and re-uploading.

### F.6 Cmdstanr backend toggle
Once the backend exposes it (`backend = c("rstan", "cmdstanr")`), add a checkbox in the Advanced section of every Bayesian sidebar. Will substantially reduce fit time on the VPS.

### F.7 Session resume
If the user accidentally closes the tab or the connection drops mid-Bayesian-fit, the in-flight future is lost. Storing the future-promise handle keyed by a session-ID cookie would let the user reconnect to a still-running fit.

**✅ Resolved v26.06.006.** Process-level job registry keyed by a
256-bit CSPRNG cookie sid; futures live outside the per-session
reactive scope so a fit survives tab disconnect / NPM proxy timeout
/ browser sleep. Resume banner surfaces outstanding jobs on
reconnect (only when the sid has any — no UI noise on fresh
sessions). All six Bayesian fits (TG, BW, Survival, Synergy, DR,
TWM) are routed through the registry; each has a per-kind resume
handler that re-attaches the poller on Resume click and switches
to the matching module tab. Cancel button per job. TTL sweeper
evicts completed jobs after 30 min and abandoned running jobs
after 4 h. Multi-user isolation enforced at four independent
layers (registry keying, server-side cookie reading, sid leak
audit, worker non-eval). 210/210 tests pass including 75
F.7-specific assertions.

---

## G. Status Summary — Round 2

| ID | Issue | Severity | File / Area | Status |
|---|---|---|---|---|
| A.1 | Orphaned `progress$set()` after modal migration | Critical | All modules | ✅ Fixed (56bfaf5) |
| A.2 | Toxicity BW `on.exit` removed Bayesian modal | Critical | `toxicity_modules.R` | ✅ Fixed (56bfaf5) |
| A.3 | Bayesian "overlay" = Shiny disconnect overlay from NPM timeout | Critical | Deployment | ⚠ Partial — script fixed; live NPM still needs PUT |
| A.4 | Stale browser cache of CSS/JS | Major | `app.R` | ✅ Fixed (5e63303) |
| B.1 | Silent rebuild failures (CRAN miss → broken image swap) | Critical | VPS `rebuild-shiny.sh` | ✅ Fixed (VPS `d2034c5`) |
| B.2 | NPM hosts missing websocket + advanced_config | Critical | VPS `06-configure-npm.sh` | ✅ Fixed (VPS `043419e`) |
| C.1 | Server function length still growing | Architecture | Multiple | ✅ Fixed v26.06.003 + v26.06.004 (all 5 mega-files split by concern; 7 new helper files; ~2,100 LOC moved out of main module files) |
| C.2 | Module file naming inconsistency | Minor | Module files | ✅ Fixed v26.06.002 (mod_*.R everywhere) |
| C.3 | Mega-files concentrate render handlers | Architecture | Multiple | ✅ Fixed v26.06.003 + v26.06.004 (downloads + plots extracted for the four module files; helpers.R split by category) |
| D.x | Diagnostic surfacing | Major | All Bayesian / Survival tabs | ✅ Fixed v26.06.001 (bayes_R2, ppc_coverage, NUTS banner, LOO panel, cox.zph tab, C-index, P_direction column) |
| E.1 | Power Analysis missing Summary tab | Minor | `power_modules.R` | ✅ Fixed v26.06.001 |
| E.2 | CHANGELOG missing v26.05.020 entry | Minor | `CHANGELOG.md` | ✅ Fixed v26.06.001 |
| E.3 | No startup warning for missing NPM timeout | Minor | `app.R` | ✅ Fixed v26.06.002 |
| F.1 | Reproducible report export | Enhancement | New | ✅ Fixed v26.06.002 (HTML report download) |
| F.2 | Save/load analysis presets | Enhancement | New | ✅ Fixed v26.06.002 |
| F.3 | Side-by-side frequentist vs Bayesian | Enhancement | New | ✅ Fixed v26.06.002 (consecutive-run compare) |
| F.4 | Per-mouse drill-down | Enhancement | New | ✅ Fixed v26.06.002 |
| F.5 | Re-run with subset | Enhancement | New | ✅ Fixed v26.06.002 (Subset Filter) |
| F.6 | Cmdstanr backend toggle | Enhancement | All Bayesian sidebars | ✅ Fixed v26.06.001 |
| F.7 | Session resume | Enhancement | `app.R` + 11 new files | ✅ Fixed v26.06.006 (job registry + cookie identity + resume banner; all 6 Bayesian fits resumable; 75 F.7-specific test assertions) |

---

## H. Cross-reference to backend deep audit (2026-05-30)

The backend `mouseExperiment` CODE_REVIEW.md Round 2 includes a Section G "Deep Audit Findings" (added after this review). Three findings there directly affect dashboard behavior:

- **Backend G.1** — `bayesian_synergy()` and `bayesian_synergy_over_time()` wrap `brms::brm()` in `suppressWarnings()`. The dashboard's Synergy "MCMC Diagnostics" tab consequently never displays divergent-transition warnings even when the underlying fit has them. Resolving the dashboard side is blocked on backend G.1.
- **Backend G.4** — `bayesian_synergy()` doesn't return `pp_check_plot`, `prior_posterior_plot`, or `mcmc_trace_plot`. The dashboard's Synergy "Posterior Check" / "MCMC Diagnostics" / "Prior vs Posterior" tabs are wired up but **empty** because the backend never populates these outputs.
- **Backend G.5** — `bayesian_therapeutic_window()` returns **none** of the standard Bayesian diagnostic plots. The Toxicity TWM Bayesian results show no diagnostic tabs by design.

When the corresponding backend gaps are filled, the dashboard side requires no code change — the tab wiring is in place and will start populating automatically.

The misleading dose-response model-type label (backend G.2) and the inhibition-only constraint of `bayesian_dose_response()` (backend G.3) propagate directly to the Dose Response module's "Summary" tab — `model_type = "inhibition"` or `"stimulation"` is shown to the user as a categorical result there. These will need a coordinated backend+dashboard fix.

---

## I. Test suite review (2026-05-30)

**Scope:** all 6 files under `tests/testthat/` (895 LOC total, 89 `expect_` calls). Compared to the backend's 4,498 / 458 — the dashboard's test surface is much thinner per LOC of source code.

### I.0 Strengths

1. **`test-shiny_app.R` uses `shinytest2`** — proper headless-browser integration testing with graceful skips when `chromote` / Chrome are unavailable. The infrastructure is in place even if only the Data tab is exercised today.
2. **`test-helpers.R` exercises internal helpers** through namespace access (`getNamespace("mouseExperimentDashboard")` + `get()`), so private helpers can be tested without exporting them.
3. **Test fixtures load real demo data** (`combo_treatment_synthetic_data`) — exercises realistic shapes rather than minimal hand-built data frames.

---

### I.1 `expect_true(TRUE)` placeholders pad test counts without testing anything — **Major**
**File:** `test-module_servers.R` lines 102–142

```r
test_that("dose_response_server initialises without error", {
  shiny::testServer(dose_response_server_fn, args = list(rv = make_test_rv()), {
    expect_true(TRUE)   # ← no assertion about the module's behavior
  })
})
```
Five module servers (`dose_response`, `power`, `survival`, `synergy`, plus a stub for analysis_log) are "tested" with `expect_true(TRUE)` inside `testServer()`. These pass as long as the module initialises without throwing. They do **not** verify a single reactive value, output, or observer. They inflate the green-check count without checking anything.

**Fix:** for each module, drive `session$setInputs(...)` to simulate the "Run analysis" click and assert `rv$X_results` is populated. The fixtures already include the data needed for this.

### I.2 `tumor_growth_server`, `toxicity_server`, `data_upload_server` not tested — **Major**
The three largest server modules in the repo (2,367 / 2,069 / 653 LOC) have **zero** `testServer()` entries in `test-module_servers.R`. These contain the bulk of the dashboard's reactive logic, including all of this session's modal / overlay / orphaned-progress bug surface. The auto-switch-to-Summary behaviour added in `8ee551d` could have been broken at any point and the test suite wouldn't notice.

### I.3 No tests for this session's new features — Missing
- `show_analysis_modal()` ([helpers.R:16](R/helpers.R#L16)) — modal generation, message/detail wiring
- `updateTabsetPanel(..., selected = "Summary")` auto-switch — assert `session$getReturned()$X_results_tabs == "Summary"` after analysis completes
- Cache-busting `?v=<version>` query string in static asset includes — assert UI HTML contains the version query
- MutationObserver in `overlay.js` — untestable from R, but the integration test could assert the rendered DOM has `.modal-backdrop` with `background-color: transparent` after a modal opens

### I.4 `analysis_log_server` tests still present despite Log tab removal — Minor (verify-and-prune)
**File:** `test-module_servers.R` lines 76–97

`analysis_log_server` is still defined in `R/log_session_modules.R` so the tests still pass — but the Log tab itself was removed in v26.05.020 ("remove Log tab; rewrite Help tab"). Either restore the tab, or delete the unused server module + its tests.

### I.5 `shinytest2` coverage is single-tab — Minor
`test-shiny_app.R` builds an `AppDriver` and asserts the Data tab renders. It doesn't drive a workflow (upload → configure → run an analysis → inspect results). A second test exercising "upload combo data → click Run on Tumor Growth → verify results panel populates" would catch the entire integration-level bug class that this session has been chasing (modal stuck, no results, etc.).

### I.6 No test asserts the dashboard's NPM-timeout / cache-bust deployment contract — Architecture
Section A.3 / A.4 of this review document deployment-side bugs that the test suite cannot catch (NPM proxy timeout, browser-cached assets). The companion VPS repo has rebuild and verification scripts; a smoke test that pings `staging.mouseexperiment.com` from CI after a webhook fires and asserts a Bayesian analysis can complete would close that loop. Out of scope for the dashboard test suite proper, but worth tracking.

### I.7 `make_test_rv()` doesn't pre-populate `rv$X_results` — Minor
The module-server fixture sets `tg_results = NULL`, `surv_results = NULL`, etc. So any test of the **post-completion** UI behaviour (auto-switch to Summary, status alert content, button enable/disable on results-ready) has no way to set up the precondition. A `make_test_rv(with_results = c("tg", "surv"))` variant would unlock those tests.

---

## I summary

| ID | Issue | Severity |
|---|---|---|
| I.1 | `expect_true(TRUE)` placeholders in module-server tests | **Major** | ✅ Fixed v26.06.001 |
| I.2 | `tumor_growth_server` / `toxicity_server` / `data_upload_server` untested | **Major** | ✅ Fixed v26.06.002 (init tests) |
| I.3 | No tests for `show_analysis_modal()`, auto-switch, cache-busting | Missing | ✅ Fixed v26.06.002 |
| I.4 | `analysis_log_server` tests outlive the UI tab | Minor (process) | ✅ Fixed v26.06.001 (function + tests removed) |
| I.5 | shinytest2 covers only the Data tab | Minor | ✅ Fixed v26.06.002 (end-to-end workflow added) |
| I.6 | No CI smoke test against deployed staging | Architecture | ❌ Open (VPS-side) |
| I.7 | `make_test_rv()` doesn't support `with_results` | Minor | ✅ Fixed v26.06.001 |

---

# Review Round 3 (v26.06.013 — 2026-07-29)

**Reviewer:** Claude Code
**Scope:** Statistical correctness of what the dashboard *displays*, as distinct from what the backend computes. Round 1 and Round 2 covered Shiny mechanics, deployment, and module harmonisation thoroughly; this round looks specifically at the places where the dashboard re-derives or re-labels a statistic rather than passing the backend's through, since those are the points where the UI can disagree with the package.

Round-1 and Round-2 items spot-check as resolved. All findings below are net-new. Backend findings referenced as `R3.n` live in `mouseExperiment/CODE_REVIEW.md` Round 3.

---

## R3.D1 The dashboard silently discards the backend's pairwise comparisons and recomputes its own — **Major**
**File:** `R/mod_tumor_growth.R:633-653`

For every non-AUC model the dashboard ignores `result$pairwise_comparisons` and rebuilds contrasts from `result$model`:

```r
emm <- emmeans::emmeans(result$model, as.formula(paste("~", treatment_col_name)))
result$all_pairwise_comparisons <- as.data.frame(
  summary(pairs(emm, adjust = "tukey"), infer = c(TRUE, TRUE)))
result$reference_comparisons_dunnett <- as.data.frame(
  summary(emmeans::contrast(emm, method = "trt.vs.ctrl", ref = ref_idx, adjust = "dunnett"),
          infer = c(TRUE, TRUE)))
```

The point estimates agree with the backend's (both marginalise `Day` at its mean — the backend via an explicit `at =`, the dashboard via emmeans' default reference grid), but the **p-values do not**: the backend applies no adjustment at all (backend R3.1), the dashboard applies Tukey or Dunnett. The same contrast on the same data therefore yields one p-value in the Shiny table and a different one from `results$pairwise_comparisons` in an R session or a downstream script.

This is not primarily a dashboard defect — the dashboard's behaviour is the *more* defensible of the two — but a product that reports two different p-values for one comparison cannot be used for a regulatory submission or a manuscript. Whichever adjustment is decided on (see backend R3-F question 1), one implementation should produce it and the other should consume it.

**Fix:** have the backend apply `p_adjust_method` (backend R3.1) and have the dashboard render `result$pairwise_comparisons` directly, deleting the recomputation. That also removes the dependency on `return_model = TRUE` that currently forces full model objects into the reactive values Round 2 §3.3 tried to keep out.

---

## R3.D2 `adjust = "dunnett"` is silently downgraded to the `dunnettx` approximation — **Minor**
**File:** `R/mod_tumor_growth.R:650`; help text at `:1098`

`"dunnett"` is not one of emmeans' recognised `adjust` values. emmeans falls back to `"dunnettx"` — its fast approximation to the Dunnett procedure — without erroring. Verified by execution; the returned summary footnote reads:

```
Conf-level adjustment: dunnettx method for 2 estimates
P value adjustment: dunnettx method for 2 tests
```

The help text at `:1098` tells the user "adjusted p-value (Dunnett correction for 'vs Reference')". `dunnettx` is documented in emmeans as an approximation that is accurate for two-sided tests with equal correlations and degrades otherwise; the exact procedure is `adjust = "mvt"`, which is affordable at the group counts these studies use.

**Fix:** use `adjust = "mvt"` and label it "Dunnett (exact, multivariate-t)", or keep `dunnettx` and say "Dunnett (approximate)" in the help text. Either is fine; claiming exact Dunnett while computing the approximation is not.

---

## R3.D3 `p_adjust_method` is not exposed anywhere in the UI — **Minor**
`tumor_growth_statistics()` documents four multiplicity options. No `selectInput`, `radioButtons`, or any other binding for it exists across `R/` — grep returns nothing. Dashboard users get whatever each code path hardcodes: Tukey or `dunnettx` for lme4/GAM (D1, D2), Bonferroni over all pairs for AUC. The "Post-hoc Comparisons" dropdown chooses the *family* (reference vs all-pairs) but not the correction.

Worth deciding alongside backend R3.1 — if the backend gains a working `p_adjust_method`, the UI control is a two-line addition and the comparison type dropdown becomes a natural place for it.

---

## R3.D4 The Tumor Growth "All Comparisons" tab is permanently empty for GAM models — **Major**
**File:** `R/mod_tumor_growth.R:633`, table renderer at `:1224-1231`

For `model_type = "gam"`, `result$model` is a `gamm4` object — a bare `list(mer =, gam =)` with class `c("gamm4", "list")`. `emmeans::emmeans()` does not support it:

```
ERROR: Can't handle an object of class  "gamm4"
```

(Verified by execution. The backend works around this internally by handing emmeans the patched `$gam` stub — see backend R3.4 — but the dashboard passes the whole `gamm4` object.)

The error is caught by the enclosing `tryCatch` at `:566` and routed to `log_caught_error()`, so:
- `comp_type == "tukey"` → `all_pairwise_comparisons` is `NULL` → the tab renders `"Comparisons not available for tukey"`. There is no path by which it can ever populate for a GAM model.
- `comp_type == "reference"` → falls back to `rv$tg_results$pairwise_comparisons`, which for GAM is the backend's per-day smooth-difference table (one row per contrast × quantile day). That table is rendered under the header `"Comparisons: vs Reference Group"`, which does not tell the user the rows are per-day rather than one-per-contrast.

**Fix:** pass `result$model$gam` when `model_type_used == "gam"` (after applying the same class/`$call` patch the backend uses), or short-circuit the recomputation for GAM and render the backend's table with a header that says the contrasts are evaluated at five study-day quantiles.

---

## R3.D5 Survival per-animal reduction takes the last *row*, not the last *day* — **Major**
**File:** `R/mod_survival.R:645-657`

```r
animal_rows <- data[mask, , drop = FALSE]
survival_rows[[i]] <- animal_rows[nrow(animal_rows), , drop = FALSE]
```

The reduction to one row per animal — which is what makes the backend Cox fit valid (backend R3.14) — assumes the uploaded file is sorted by day within animal. There is no `order(animal_rows[[time_col]])` before the subset. A file exported in a different order (grouped by measurement session, sorted descending, or with a re-measure appended out of order) yields a survival time and event indicator taken from an arbitrary observation, and the analysis proceeds with no warning.

This is the same defect class as backend Round 1 §1.7 (`aggregate(..., FUN = function(x) x[1])` without an ordering guarantee), which was fixed in three backend files by filtering on `Day == min(Day)` first. The dashboard's reduction never got the same treatment.

**Fix:** `animal_rows <- animal_rows[order(as.numeric(animal_rows[[time_col]])), ]` before taking the last row. Consider also asserting that the resulting event indicator is the animal's maximum, since an event flagged at an intermediate day would otherwise be dropped.

---

## R3.D6 Displayed scale metadata is what the user *asked for*, not what the backend *did* — **Minor**
**File:** `R/mod_survival.R` n/a; `R/mod_tumor_growth.R:748-750`

```r
# (the backend result$methods list is not reliably populated for all paths)
result$transform_used  <- input$tg_transform   %||% "none"
result$model_type_used <- input$tg_model_type  %||% "lme4"
```

The comment is correct — no frequentist backend path returns `transform_used` (backend R3.29) — and reconstructing it from the input is a reasonable workaround. It is nonetheless recording the *request*, not the *outcome*, and the two diverge: the AUC path ignores `transform` entirely and computes on raw volumes (backend R3.16), so an AUC run with the default `transform = "log"` is tagged `transform_used = "log"`.

Today this is benign — both consumers of the field (`:1239` Fold_Change, `:1284` scale-note) are gated on `model_type_used == "lme4"`, so the AUC path escapes both. It is a latent trap: the next feature that reads `transform_used` without that gate will mislabel AUC output. Once backend R3.29 lands, read the field from the result instead.

---

## R3.D7 Analysis-parameter provenance is not captured in the HTML report — **Enhancement**
**File:** `R/helpers_report.R`

Round 2 §F.1 added HTML report export. For the report to serve the audit-artifact purpose stated there, it needs to record the settings that produced the numbers — reference group, transform, model type, random-effects specification, cage handling, endpoint day, prior strength and MCMC settings for Bayesian runs, and the multiplicity adjustment actually applied (which, per D1–D3, is currently not the one any argument names). A short "Analysis provenance" table at the head of the report would make each export self-describing and independently reproducible.

This becomes more valuable, not less, once the backend multiplicity work lands, because it is the natural place to state which correction was used.

---

## R3.D8 Round 2 §A.3 remains genuinely open — worth re-confirming against live NPM
**Cross-repo:** `mouseExperimentVPS/scripts/06-configure-npm.sh`

Round 2 §A.3 is marked "⚠ Partial — script fixed; live NPM still needs PUT". The provisioner fix is committed (`043419e`), but the live proxy host update was to be applied out-of-band via the API. Nothing in either repo records whether that was done. Since the symptom (Shiny's `#shiny-disconnected-overlay` appearing mid-fit) only reproduces on Bayesian runs longer than 60 s, it will not surface in routine testing.

**Suggestion:** add a one-line check to `scripts/99-diagnose-shiny-dashboard.sh` that reads the live NPM host config and asserts `proxy_read_timeout` ≥ 3600 and `allow_websocket_upgrade` is true, so the state is verifiable on demand rather than remembered.

---

## R3.D9 Status table — Round 3

| ID | Issue | Severity | File | Status |
|---|---|---|---|---|
| R3.D1 | Dashboard recomputes contrasts; p-values disagree with the package's | **Major** | `mod_tumor_growth.R:633` | ✅ Fixed v26.07.001 |
| R3.D2 | `adjust = "dunnett"` silently becomes `dunnettx`; help text claims exact Dunnett | Minor | `mod_tumor_growth.R:650` | ✅ Fixed v26.07.001 |
| R3.D3 | `p_adjust_method` not exposed in the UI at all | Minor | — | ✅ Fixed v26.07.001 |
| R3.D4 | "All Comparisons" tab permanently empty for GAM (emmeans can't take a `gamm4` object) | **Major** | `mod_tumor_growth.R:633` | ✅ Fixed v26.07.001 |
| R3.D5 | Survival reduction takes last row by position, not by max day | **Major** | `mod_survival.R:655` | ✅ Fixed v26.07.001 |
| R3.D6 | `transform_used` records the request, not what the backend did | Minor | `mod_tumor_growth.R:749` | Open |
| R3.D7 | HTML report lacks an analysis-provenance block | Enhancement | `helpers_report.R` | 🟡 Partial v26.07.001 (provenance in summary; report block pending) |
| R3.D8 | Round 2 §A.3 live-NPM state unverified and unrecorded | Process | VPS repo | Open |

---

## R3.D10 Cross-reference — backend findings that surface in the UI

These need no dashboard code change, but they change what users are looking at and should be tracked alongside the backend fix:

| Backend finding | What the dashboard currently shows |
|---|---|
| R3.1 (`p_adjust_method` inert) | Unaffected for TG (D1 recomputes), but no UI control exists to choose the correction |
| R3.3 (AUC over unequal windows) | The AUC results tab presents group AUCs as comparable; `Last_Day` is shown in the per-animal table but nothing flags the imbalance |
| R3.4 (BW GAM emmeans NULL) | Toxicity → GAM → marginal-means table renders empty, indistinguishable from "no effect" |
| R3.5 (survivor-biased TGI) | TWM, Synergy, and Efficacy–Toxicity tabs all display endpoint-day TGI with no denominator *n* |
| R3.7 (synergy point estimates) | Synergy tab shows a categorical "Strong Synergy / Additive / Antagonistic" call with no interval |
| R3.16 (AUC ignores `transform`) | Methods text in the results panel and HTML report states the transform was applied |
| R3.18 (TWM discontinuity) | TWM table is sorted by a metric that changes scale across the noise floor |

---

## R3.D11 Maintainer decisions (2026-07-29) — dashboard consequences

The backend §R3-F questions were answered; see `mouseExperiment/CODE_REVIEW.md` §R3-G. Three answers land directly on this repo.

**G.1 — "Users should be able to define the comparison of interest."** This settles R3.D1, R3.D2, and R3.D3 as one piece of work rather than a choice:

- The backend gains `comparison_family = c("vs_reference", "all_pairs", "custom")` plus a *working* `p_adjust_method`. The dashboard's existing "Post-hoc Comparisons" dropdown becomes the `comparison_family` control; add a second `selectInput` for the correction (closing R3.D3).
- **Delete the recomputation block at `mod_tumor_growth.R:633-653`** and render `result$pairwise_comparisons` directly. That removes the package/dashboard p-value divergence (R3.D1), removes the `emmeans` call that cannot handle a `gamm4` object (R3.D4, which then fixes itself), and removes the dependency on `return_model = TRUE` that currently forces full model objects into reactive values against Round 2 §3.3.
- When the backend implements the correction, use `adjust = "mvt"` for exact Dunnett rather than `"dunnett"`, which emmeans does not recognise and silently downgrades to the `dunnettx` approximation (R3.D2). Update the help text at `:1098` to name whichever is used.
- Display `result$comparison_family` / `result$p_adjust_method_used` in the table header, so the panel states what was applied instead of asserting a default.

**G.2 — cage structure.** The maintainer confirms the standard design is one treatment per cage with ≥2 cages per arm, which is estimable as `(1 | Cage)` and is currently being discarded (backend R3.17, upgraded to Critical). Once the backend fits it, the dashboard should surface the **cage ICC** in the model-summary panel — it is the number that tells a user how much the clustering mattered, and there is nowhere to see it today. The cage-handling selector should also be relabelled to match the new structure-aware `"auto"` default.

**G.3 — dropout.** Endpoint metrics move from raw survivor means to model-based estimated marginal means at day T. Dashboard work:
- Add **n at risk at day T** as a column beside every TGI / TWM / synergy figure. This is the number that currently makes survivor bias invisible.
- Add a dropout curve (time-to-removal per arm) to the Tumor Growth results panel, not only the Survival tab — the point is for it to be visible next to the efficacy result that depends on it.
- The Synergy tab's categorical verdict should not render at all until the backend supplies an interval for it (backend G.6).

**G.4 — extrapolation.** `extrapolation_points` is being removed from `tumor_growth_statistics()`. If any sidebar control is currently bound to it, remove that control. The `plot_tumor_growth()` extrapolation is explicitly retained, so plot-side projection controls stay — but projected segments should be rendered visually distinct (dashed, or a shaded band) so a reader cannot mistake a projection for measured data.

| ID | Item | Blocked on | Status |
|---|---|---|---|
| R3.D1 | Delete contrast recomputation; render backend output | backend `comparison_family` | Agreed, blocked |
| R3.D2 | Use `mvt`, or relabel as approximate Dunnett | — | Agreed |
| R3.D3 | Expose comparison family + correction in the UI | backend `p_adjust_method` | Agreed, blocked |
| R3.D4 | Resolved as a side-effect of R3.D1 | R3.D1 | Agreed |
| R3.D5 | Sort by day before taking the per-animal last row | — | Agreed, unblocked |
| — | Surface cage ICC in model summary | backend G.2 | New |
| — | n-at-risk column + dropout curve on efficacy panels | backend G.3 | New |


---

## R3.D12 Round 3 closed (v26.07.001, 2026-07-29)

R3.D1–R3.D5 closed; R3.D7 partial (the analysis-provenance block is in the model
summary, not yet in the HTML report export); R3.D8 remains a VPS-side item.

The theme of this round was **the dashboard silently disagreeing with the
package**. Four of the five findings traced to one cause: the UI re-derived a
statistic instead of rendering the backend's, so a backend defect (unadjusted
p-values) and a UI workaround (hardcoded Tukey/Dunnett) coexisted and produced
two different answers for the same contrast. Deleting the recomputation once the
backend was fixed closed R3.D1, R3.D2 and R3.D4 together.

**Worth carrying forward:** when the dashboard re-derives anything the package
already returns, that is a divergence waiting to happen. The remaining instances
are worth auditing on the same principle — grep for `emmeans::`, `lme4::` and
`survival::` calls in `R/mod_*.R` and ask whether the backend already provides
the answer.

### New dashboard work created by the backend changes

| Backend addition | Dashboard status |
|---|---|
| `comparison_family` / `p_adjust_method_used` | ✅ selector added, provenance shown |
| `cage_analysis$icc`, `cage_structure`, `cage_placement` | ✅ CAGE STRUCTURE block |
| `attrition` (N enrolled / at endpoint / % lost) | ✅ v26.07.004 — shown on the Synergy Metrics and Therapeutic Window tabs, with a warning when the survivor estimand is combined with real attrition |
| `endpoint_method` on TWM / synergy / weight-corrected TGI | ✅ v26.07.004 — selector in both the Synergy and Toxicity sidebars, defaulting to model-based, with the survivor-bias magnitude quoted in the help text |
| `synergy_ci`, `twm_ci` (bootstrap intervals) | ✅ v26.07.004 — interval table under the synergy verdict with an explicit "does the Bliss-excess interval exclude zero?" callout, and a ranking-resolvability note above the TWM table |
| `auc_permutations` / `perm_p_value` | Open — still no UI control; off by default so nothing changes silently |
| `cuminc` / `censoring_summary` (competing risks) | 🟡 v26.07.004 — the KM tab now warns that 1 − KM overstates incidence when animals were removed early, and breaks the censoring down by type. The Aalen-Johansen curve itself is still not plotted |
| `variance_test` (Brown-Forsythe on AUC) | ✅ v26.07.004 — on the ANOVA tab, framed as "was the unequal-variance choice needed" rather than a gate |

These are additions rather than defects — nothing displays a wrong number — but
the bootstrap intervals and the attrition counts are the two that most change how
a reader should interpret what is already on screen, and are the natural next
piece of dashboard work.


---

# Review Round 4 — Dependency declarations (2026-07-29)

Companion to `mouseExperiment/CODE_REVIEW.md` §R4. Two findings here.

## R4.D1 Async Bayesian dependencies were optional for a feature the app advertises — **Major**

`brms`, `bayesplot`, `future` and `promises` sat in `Suggests`, and five code paths
guarded on them to show a "Bayesian mode requires 'brms', 'future', and 'promises'.
Install with: install.packages(...)" danger banner.

That is a reasonable pattern for a genuinely optional feature. It is the wrong
pattern here: every analysis module offers a Bayesian model type in its UI, the
job registry (F.7) exists solely to keep `future` promises alive across
reconnects, and the VPS image installs all four. The dependency was optional in
`DESCRIPTION` and mandatory in practice.

All four moved to `Imports`; the five dead banner blocks removed. The
`mouseExperiment` floor was raised to `>= 0.10.0`.

## R4.D2 The test suite was running against a stale installed backend — **Major (process)**

`devtools::load_all()` resolves `Imports: mouseExperiment` against the **installed**
library, not the sibling source tree. The installed backend was **0.4.10**; the
source under review was 0.9.0 — six releases and the whole of Round 3 apart.

So "210 passing" was 210 passing against a backend from before any of the Round 3
work. Every dashboard test that exercises backend behaviour — comparison families,
cage handling, endpoint estimands, the result-object fields the render paths read
— was testing code that had not changed.

This is the dashboard's version of the brms skip in the backend §R3-L: a green
suite that is not testing what it appears to test, for a structural reason rather
than a coding mistake.

**What re-running actually showed.** After installing the backend from source
(0.10.0) the suite reports **210 passing / 0 failing / 0 errors — the same 210 as
against 0.4.10.** No break was hidden. But that identical result is itself the
finding: a suite whose outcome does not move across six releases of its primary
dependency — including Round 3's changes to comparison families, cage handling,
endpoint estimands, and the result-object fields the render paths read — is not
exercising that dependency. This corroborates §I.1 (`expect_true(TRUE)`
placeholders) and §I.2 (the three largest server modules untested) from a
different direction: the coverage is thin enough that the version skew could not
have mattered.

So R4.D2 did not hide a defect this time. It would hide the next one, and the
"210 either way" result means the suite currently cannot tell the difference
between a working backend and a six-release-stale one.

**Fix required (not yet applied):** install the backend from source before running
dashboard tests. Options, roughly in order of preference:

1. A `Remotes:` entry or a documented test preamble
   (`devtools::install("../mouseExperiment", upgrade = "never")`).
2. A CI step that installs the sibling source, so the version skew cannot recur
   silently.
3. At minimum, a `testthat` helper that compares
   `packageVersion("mouseExperiment")` against the `DESCRIPTION` floor and fails
   loudly on skew — the same "fail rather than skip" principle that Round 3 applied
   to brms.

| ID | Issue | Severity | Status |
|---|---|---|---|
| R4.D1 | brms/bayesplot/future/promises optional for an advertised feature | Major | ✅ Fixed v26.07.002 |
| R4.D2 | Test suite ran against installed backend 0.4.10, not the source | Major (process) | ✅ Fixed v26.07.003 |


---

## R4.D3 R4.D2 resolved (v26.07.003, 2026-07-30)

`test-backend-version.R` compares the installed `mouseExperiment` against the
sibling source's `DESCRIPTION` and **fails** when the installed copy is older.

It fails rather than skips, on purpose. Skipping is the exact failure mode Round 3
spent its effort removing (§R3-L): a suite that quietly declines to test something
reports green, and green against the wrong backend is worse than no run at all. The
failure message names the fix (`devtools::install("../mouseExperiment")`) so it is
actionable rather than merely alarming.

Verified in both directions: passes when the versions align, and fails when the
sibling `DESCRIPTION` is temporarily bumped ahead of the installed copy. When there
is no sibling checkout — a CI job installing from GitHub, say — it succeeds quietly,
since `load_all()` has already enforced the `DESCRIPTION` floor by that point and
there is nothing to compare against.

`CONTRIBUTING.md` now leads its Testing section with the install step.

**What this does not fix:** the coverage problem underneath. The same 210 tests pass
against backend 0.4.10 and 0.11.0, so the suite still cannot distinguish a working
backend from a stale one — it can now only detect that the *versions* differ. §I.1
(`expect_true(TRUE)` placeholders) and §I.2 (the three largest server modules
untested) remain the substantive gap, and closing them is what would make this
check redundant.


---

## R4.D4 Backend outputs surfaced (v26.07.004, 2026-07-30)

Round 3 and Round 5 added six result fields the UI ignored. Five are now shown; one
is partially done.

**Synergy Metrics tab** — bootstrap interval table under the verdict, with an
explicit callout stating whether the Bliss-excess 95 % interval excludes zero. That
callout is the point of the whole exercise: the verdict above it is a categorical
label derived from three group means, and a "Strong Synergy" whose interval spans
zero was previously indistinguishable on screen from one that does not. Plus an
attrition table and an endpoint-estimand selector.

**Therapeutic Window tab** — the table is sorted by TWM, so it reads as a ranking.
A note now compares the top two arms' intervals and says outright when they overlap,
because with typical group sizes they usually do and the order is not resolvable.
Plus attrition, endpoint day, and the estimand used.

**Weight-loss KM tab** — a competing-risks warning when animals left before the
study ended, stating that 1 − KM *overstates* weight-loss incidence in that case,
plus a censoring breakdown by type (event / administrative / early removal). The
Aalen-Johansen curve is available in the result object but not yet plotted, so this
row stays partial.

**ANOVA tab** — the Brown-Forsythe variance check, framed as "was the
unequal-variance choice needed" rather than as a gate, since the AUC path uses
Welch's tests either way.

**ANOVA tab, randomisation test** — the H.3 result, with three interpretive aids
that matter more than the p-value itself: a warning when the asymptotic p-value is
more than 100× smaller than the randomisation one (which is the norm at these group
sizes), a warning when the design cannot attain p < 0.05 at all, and the unit and
number of exchangeable units stated plainly.

### The part worth keeping

`test-backend-contract.R` — 39 assertions that call the backend the way the dashboard
does and check that every field the new renderers read exists with the expected
name. This is the first test in this repo that verifies the dashboard/backend
interface rather than module initialisation.

It closes a specific gap: the renderers above read roughly twenty field names off
result objects, and nothing checked them. A rename upstream would have silently
blanked a panel — the same class of failure as §R3.D4 (the GAM "All Comparisons" tab
that could never populate) and §R4.D2 (a suite passing against a six-release-stale
backend). Suite went from 212 to 251 passing.

**Still the substantive gap:** §I.1 and §I.2. The contract tests verify the
*interface*; they do not exercise the Shiny reactives, so the render logic itself
(the interval comparison, the overlap check, the floor warning) is still only
verified by reading it. Driving those through `shiny::testServer` is the next step
and is what would finally make this suite able to detect a broken backend.
