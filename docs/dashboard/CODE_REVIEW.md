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
| R3.D6 | `transform_used` records the request, not what the backend did | Minor | `mod_tumor_growth.R:749` | ✅ Fixed v26.07.006 |
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

---

# Review Round 8 (v26.07.005 — 2026-07-30)

**Scope:** The Power module, reached from the backend side. The backend's first
working coverage baseline named `apriori_power_simulation()` at 0.00 %; testing it
turned up two findings that live here rather than there.

## D.4 `%||%` fallbacks encoded the abandoned percent scale — Major

`mod_power.R` passed fallback values from the parameterisation the backend
abandoned in §2.2, when its data-generating model moved from the linear to the log
scale:

```r
control_growth_rate = input$pw_lmm_ctrl_rate     %||% 15,   # log rate; s/b 0.15
treatment_effect    = input$pw_lmm_tx_effect     %||% 10,   #            s/b 0.10
random_slope_sd     = input$pw_lmm_rand_slope_sd %||% 3,    #            s/b 0.05
random_intercept_sd = input$pw_lmm_rand_int_sd   %||% 10,   #            s/b 0.20
residual_sd         = input$pw_lmm_resid_sd      %||% 10,   #            s/b 0.10
```

All five input IDs exist and their UI defaults are correct, so a fallback fires
only when an input reads `NULL` — transiently before the client syncs, or if an ID
is ever renamed. Had one fired, `control_growth_rate = 15` makes `exp(15 · t)`
overflow to `Inf` by day 14 and the fit dies with a message pointing nowhere near
the cause.

Latent rather than live, but free to fix and now held down by
`test-power-params.R`, which additionally asserts each fallback equals the UI
default sitting beside it so the two cannot drift apart again.

## D.5 The Method panel pointed users at the one dial that does nothing — Major

The panel described `random_intercept_sd` as *"Between-animal variability in
baseline volume (log scale). **Higher values require larger N.**"*

The claim is false, and measurably so. Backend §K.14 establishes that power is
invariant to `random_intercept_sd` across a 30-fold range (0.05 → 1.50, power
0.6667 throughout), because the `Treatment:Day` slope contrast is orthogonal to
per-mouse intercepts, which `(Day | ID)` absorbs exactly. `baseline_sd` is inert
for the same reason and is aliased with it.

Meanwhile `random_slope_sd` — which moves power from 1.00 to 0.15 over
0.01 → 0.15 on the same design — was described only as "between-animal variability
in growth rate", with no indication that it dominates.

The practical harm is a user calibrating carefully against historical data, seeing
the answer not move, and concluding the tool is broken — or worse, not checking.
Both entries rewritten: the inert parameters now say so and explain why, and
`random_slope_sd` is labelled the dominant driver with the worked numbers.

## D.6 Dropout exposed in the UI (backend K.15)

The backend gained a `dropout_limit` argument so power can be planned against an
IACUC euthanasia limit rather than an idealised complete-data study. Added as
"Euthanasia volume limit (blank = none)", mapping blank/`NA`/non-positive to `Inf`
so the backend's validation never sees a malformed value.

The Summary tab now reports realised attrition, and — importantly — states the
assumption *even when no dropout is modelled*, so a complete-data figure is never
read as unconditional. Ignoring dropout costs roughly 6 points of power at 30 %
attrition and 12 at 40 %.

| ID | Issue | Severity | Status |
|---|---|---|---|
| D.4 | Fallbacks on the abandoned percent scale | Major | ✅ Fixed v26.07.005 |
| D.5 | Method panel named the wrong power driver | Major | ✅ Fixed v26.07.005 |
| D.6 | No way to plan against dropout | Major | ✅ Fixed v26.07.005 |

## R3.D6 closed: scale metadata now records the outcome

R3.D6 was left open with the note *"Once backend R3.29 lands, read the field from
the result instead."* It has landed — `me_result_meta()` (B7.1/B7.2) makes every
analysis function declare `transform_used` and `model_type_used`, with `"none"`
meaning none was applied rather than the field being absent.

The divergence R3.D6 predicted is real and now measured: an AUC run submitted with
the default `transform = "log"` returns `meta$transform_used == "none"`, because
that path computes on raw volumes and ignores the argument (backend R3.16).
Recording `input$tg_transform` labelled that output `"log"`.

Resolution order is now `result$meta$…` → top-level field → input, so an older
backend still degrades to the previous behaviour rather than erroring. Two
contract tests pin it: one asserts the backend's AUC/lme4 divergence directly (so
the fix cannot rot if the backend changes), the other asserts the module's
resolution order.

This was benign today — both consumers are gated on `model_type_used == "lme4"` —
and is fixed because that gate is exactly the kind of thing the next feature
forgets.


Suite went from 251 to 285 passing. §I.1 / §I.2 remain the substantive gap: these
are source-level assertions, not `testServer` exercises of the reactives.

---

# Review Round 9 (v26.07.007 — 2026-07-30)

**Scope:** I.2, attacked properly this time. Closing it turned up why it had
resisted four rounds of attempts.

## D.7 Modules wiped their own result slot at init — Major

Every module server opened with an unconditional clear:

```r
tumor_growth_server <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
  rv$tg_results <- NULL          # and the same in five other modules
```

`rv` is created empty in `app.R` (`rv <- reactiveValues()`) and each module mounts
exactly once, so this never cleared anything in production. What it did do is make
the modules **impossible to drive through `testServer()` with a pre-populated
result** — the slot was wiped the instant the module initialised, and every
attempt to assert render behaviour got the "run an analysis first" placeholder
back.

That is the reason I.2 has stayed open. `test-module_servers.R` asserts
`inherits(session, "session_proxy")` and that `session$ns()` returns a non-empty
string, and it has been reported as closing I.2 twice. Those assertions verify
that `testServer` built a mock session; they never touch the module. The
placeholder tests were not lazy — the harness made anything better impossible, and
nobody traced why.

Removed in all six modules (`tumor_growth`, `power`, `toxicity`, `survival`,
`synergy`, `dose_response`). Behaviour is unchanged: a missing key on a
`reactiveValues` reads `NULL` either way. The deliberate reset on *new data load*
in `mod_data_upload.R` is a different thing and is kept.

## D.8 The data-load reset named a slot that does not exist — Major

Found immediately once D.7 made the state inspectable:

```r
reset_analyses <- function(rv) {
  ...
  rv$pw_results <- NULL      # nothing has ever written this
```

`pw_results` appears in exactly one place in the repository — that line. The Power
module uses `ap_results` (18 references). And `tox_results` was absent from the
reset altogether.

So loading a new dataset reset every analysis tab **except Power and Toxicity**,
which kept displaying results computed from the previous dataset. For Toxicity
that is unambiguously wrong. For Power it is subtler but still wrong: the a priori
calculation is mostly data-independent, except that "Seed SD from loaded data"
ties the pooled SD to whatever was loaded at the time, so a stale power curve can
be calibrated against a dataset the user has since replaced.

Same failure mode as §K.2, §R3.D4, §R4.D5 and §R8.6: **a name that matches nothing
fails silently and looks like it works.** That is now five instances in this
review from one root cause — writing a name and never asserting anything reads it.

## I.2 closed: 14 render-logic assertions

`test-render-logic.R` sets `rv$<X>_results` to a realistic backend result and reads
the rendered output back, so the branches themselves are exercised rather than
read:

- the Power summary states the no-dropout assumption, reports realised attrition
  when dropout ran, and does not claim complete data when it did not;
- a pre-K.15 result with no `attrition` field still renders (resumed sessions);
- the randomisation-test panel distinguishes "not run" from a result, and never
  prints a p-value it does not have;
- the TWM ranking note separates a genuine ranking from arms whose intervals
  overlap — with identical point estimates in both fixtures, so only the interval
  logic can produce the difference.

| ID | Issue | Severity | Status |
|---|---|---|---|
| D.7 | Modules wiped their result slot at init, blocking testServer | Major | ✅ Fixed v26.07.007 |
| D.8 | `reset_analyses()` named a dead slot; `tox_results` never cleared | Major | ✅ Fixed v26.07.007 |
| I.2 | Largest server modules had no behavioural tests | Major | ✅ Fixed v26.07.007 |
| I.1 | `expect_true(TRUE)` placeholders | Major | ✅ Fixed — superseded by the above |

Suite: 300 → 314 passing, no skips.

---

# Review Round 10 (v26.07.008 — 2026-07-30)

**Scope:** D.1, the three files over 1,000 LOC — reframed after looking at what
was actually in them.

## D.1 revisited: length was the symptom, not the defect

D.1 read as "these files are too long". The repo's own history says that framing
does not pay. `mod_tumor_growth.R` had already been split once under C.1/C.3,
from 2,659 to 1,877 LOC, and the extraction was done well — a named-list contract
rather than lexical scope. It prevented nothing: **D.7 sat at line 7 of the file
that had just been refactored.** D.8 was in a 653-line file. R3.D6 was a one-line
resolution order.

Nor are the three files one problem. `ui_modules.R` (1,563) and `helpers.R`
(1,226) contain **zero** outputs, observers and reactives — static UI and pure
functions. The server modules are the inverse: 42 and 49 render blocks against 6
and 3 observers. Different shapes, different answers.

So the split was done on the one file where it bought something, and declined
where it did not.

## D.9 `helpers.R`: 29 of 39 functions had no test referencing them — Major

That is the finding length was hiding. Split by testability seam into
`helpers_data.R` (423), `helpers_format.R` (216), `helpers_plot.R` (501) and
`helpers_ui.R` (126), following the `helpers_*.R` convention already set by
`helpers_bayes.R` / `helpers_report.R`. The partition was verified lossless
before `helpers.R` was deleted, and all 39 functions were confirmed to still
resolve. Suite was unchanged at 314 across the move, which is what makes it a
refactor rather than a rewrite.

The split is the means. **164 new assertions are the point**, and two of the
functions they cover decide statistical results rather than presentation:

**`back_transform()`** — the inverse of the modelling transform. The backend fits
`log(Volume)`, not `log1p`, so `expm1()` here would be wrong by one unit at every
point: invisible on a 2,000 mm³ tumour, badly wrong on a 20 mm³ one. Now pinned
against `exp()` explicitly, with round-trip tests, and the silent fall-through for
an unrecognised transform name documented rather than left to be discovered.

**`synthesize_censor_column()`** — fabricates the event indicator the survival
analysis consumes, by the rule *last observation before study end ⇒ event*. It is
correct for the mechanism it models and it has a limitation that cannot be fixed
in this function: an animal that survived but whose final measurement was skipped
is coded identically to one that was euthanised. The information needed to
separate them is not in the data frame. What the function does do is report
`n_flagged`, and that contract is now pinned so the count stays available for the
caller to display. This is the same informative-dropout question the review opened
with, surfacing in the dashboard's data preparation.

Also newly covered: `pick_reference_group()` (which determines what every contrast
is measured against — its substring matching promotes an arm called "Uncontrolled"
over a genuine "Vehicle", asserted rather than silently changed, since requiring
exact equality would break the common "VEHICLE-1" style), `apply_subset_filter()`
(exclusions fail safe — a mis-mapped column drops nothing rather than everything),
`detect_date_format()` (a misread layout silently yields wrong elapsed days),
`format_pval()` (a p-value displayed as "0" is a claim no test supports),
`translate_error()` (unrecognised messages pass through verbatim rather than being
replaced by a generic shrug), and `js_string()` (escaping for values interpolated
into `runjs()`).

## D.10 Declined: mechanical splits of `ui_modules.R` and the server modules

`ui_modules.R` is 1,563 lines of static UI with no reactive coupling. Nothing can
hide there, and a split buys merge conflicts and a regression window in exchange
for a smaller number.

The server modules are dominated by render blocks. Moving a renderer to another
file does not make it more correct or more testable — `testServer()` does not care
which file it was defined in, and what blocked I.2 for four rounds was
`rv$x <- NULL`, not 1,877 lines. The move that would pay there is extracting the
*decision logic* out of `renderUI` — the TWM separability rule, the Bliss-excess
verdict, the interval-overlap checks — into pure functions, which shrinks the
modules as a side effect. Recorded as the next step rather than done here.

| ID | Issue | Severity | Status |
|---|---|---|---|
| D.9 | 29/39 helpers untested; `back_transform` / `synthesize_censor_column` unverified | Major | ✅ Fixed v26.07.008 — 164 assertions |
| D.1 | Files over 1,000 LOC | Architecture | ✅ Addressed — `helpers.R` split; see D.10 |
| D.10 | Mechanical splits of `ui_modules.R` / server modules | Architecture | Declined, with reasoning |
| D.11 | Extract render decision logic into pure functions | Enhancement | Open — the split that would pay |

Suite: 314 → 478 passing, no skips.

---

# Review Round 11 (v26.07.009 — 2026-07-30)

**Scope:** Reported by the user — demo data that ships but cannot be loaded, and
Info/Help tabs describing functionality that no longer exists while omitting what
does. Both are the same failure: **prose and packaged assets have no failure
mode.** They are not executed, so they cannot be wrong in a way a test suite
notices. They just quietly stop being true, which is worse than an error because
users act on them.

## D.12 Demo data that shipped but could never be loaded — Major

`master_synthetic_data` (619 × 13) and `synthetic_data` (120 × 7) are shipped by
the backend and appear in no CSV here and in no entry of the demo loader, so
neither could ever be selected. `master` is the most complete dataset in the
package — length/width/**height**, body weight, dose, cage, necrosis flags and
survival censoring across six arms — which makes it the only demo that reaches
every module without the user supplying their own data. It was unreachable.

Both are now exported to `inst/sample_data/` and wired into the loader and the
selector. `test-demo-data-reachable.R` (32 assertions) closes the loop in both
directions: every shipped file must be reachable or explicitly classified, and
every loader entry must resolve to a file that parses. It also asserts that the
master demo still carries every optional column, since a demo that silently loses
a column stops covering a module.

**Still open:** `combo_treatment_schedule.csv` and
`dose_levels_treatment_schedule.csv` ship in this package's `inst/sample_data/`
and nothing can open them. They are dosing tables (Treatment / Day / Dose), not
measurement data — loading one as the main dataset would break every analysis —
so they are deliberately excluded from the selector and classified as such in the
test rather than left ambiguous. The backend has `plot_treatments()`, which
consumes exactly this shape and is **not exposed anywhere in the dashboard**.
Wiring it up is a feature decision, not a documentation fix, and is left for the
user to call.

## D.13 Info and Help tabs describing an app that no longer exists — Major

The specific errors, all user-visible:

| Claim | Reality |
|---|---|
| "two statistical modeling approaches" | **five** — LME4, GAMM, AUC, Bayesian LMM, Bayesian GAMM |
| "Type III ANOVA F-test" | Wald **chi-square** (`Chisq` / `Df` / `Pr(>Chisq)`) |
| AUC compared by "Kruskal-Wallis / pairwise Wilcoxon" | **Welch one-way ANOVA** plus an equal-variance test |
| "emmeans with Tukey adjustment" | family and correction are **separate**; default is exact Dunnett |
| "LOCF can extrapolate AUC to a common endpoint" | extrapolation was **removed** from the statistical path; plots only |
| cage "included as a covariate … if not collinear" | structure-aware: crossed / nested-replicated / nested-confounded |
| Power offers "Bayesian simulation" | the selector offers **two** modes; there is no Bayesian power mode |

Missing entirely: the randomisation test, the five non-Dunnett corrections, GAMM
and its k-check/concurvity diagnostics, Bayesian GAMM, scale-location and
random-effects Q–Q plots, influence measures, Pareto-k, trace and rank plots, the
endpoint estimand, bootstrap intervals, TWM ranking bounds, per-arm attrition,
competing risks and Aalen-Johansen, dose-response trend tests, four auto-detected
column roles (Date, Sex, Necrotic, Survival_Time), the Height column, and five of
the six volume formulas.

All corrected or added. Where a limitation could not be fixed it is now stated
rather than omitted — AUC penalising early-removed animals, the observed
estimand flattering the control arm, the chi-square being anti-conservative in
small samples, and the randomisation test's dependence on the user correctly
naming the randomisation unit.

The Survival Info tab was audited and found **already accurate** — Firth
correction, the proportional-hazards test, competing risks and the natural-death
option were all current. Recorded so it is not re-audited.

## D.14 What makes this not recur

`test-docs-current.R` (47 assertions) pins the checkable claims to the source of
truth: the documented model count must equal the number of `selectInput` choices;
every selectable model must be named in the Info tab; every auto-detected column
role must appear in the Help table; every module offering an endpoint estimand or
bootstrap control must document it; and the removed terms ("Type III ANOVA
F-test", "Kruskal-Wallis", "Bayesian simulation", the old LOCF sentence) must not
reappear.

This cannot verify the prose is *good* — only that its falsifiable claims are
still true. That is the part that rots silently, and it is now the part that
fails a test run.

| ID | Issue | Severity | Status |
|---|---|---|---|
| D.12 | Two demo datasets shipped but unreachable | Major | ✅ Fixed v26.07.009 |
| D.13 | Info/Help describe removed functionality, omit current | Major | ✅ Fixed v26.07.009 |
| D.14 | No guard against documentation drift | Process | ✅ Fixed v26.07.009 — 47 assertions |
| D.15 | Treatment schedules unreachable; `plot_treatments()` unexposed | Enhancement | Open — needs a product decision |

Suite: 478 → 557 passing, no skips.

---

# Review Round 12 (v26.07.010 — 2026-07-30)

## D.15 resolved: dosing schedules moved onto the plots

The backend removed `plot_treatments()` and its two schedule datasets; this side
removes the schedule CSVs and adds the replacement UI. See the backend's Round 12
entry for why removal was right and for the measurement that shaped the design —
dosing timing is not recoverable from measurement data, so the user must supply
it, which makes a text field correct rather than lazy.

**Controls.** `dose_schedule_ui()` renders a toggle (off by default), a days
field, a style selector and a message slot, and is mounted on all three
time-axis modules. `apply_dose_schedule()` is one line at each plot's return
point. The toxicity trajectory has **two** branches — spaghetti and mean ± error
— and both are annotated; wiring only one would have made the marks vanish when
the user switched plot type, which is the sort of half-application a test that
counts call sites catches and a reading does not.

**Feedback rather than silence.** `resolve_dose_schedule()` captures parser
warnings and each module renders them beside the field. A mistyped arm name
dropping its row silently is the §K.2 pattern that has now accounted for five
findings in this review (K.2, R3.D4, R4.D5, R8.6, D.8); it is not repeated here.

**Failure is non-destructive.** An annotation error is logged and the plot is
returned unannotated. A decorative feature must never be able to take out the
analysis it decorates.

**The demo exemption list is now empty.** D.12 classified the two schedule CSVs
as "shipped but intentionally not offered" so the reachability test could tell
that apart from an oversight. With the files removed, `SCHEDULE_FILES` is
`character(0)` and every file in `inst/sample_data/` is loadable — which is the
stronger invariant and the one now enforced.

| ID | Issue | Severity | Status |
|---|---|---|---|
| D.15 | Schedule data unreachable; no dosing context on plots | Enhancement | ✅ Resolved v26.07.010 |

Suite: 557 → 621 passing, no skips.

---

# Review Round 13 (v26.07.011 — 2026-07-30)

## D.16 The dosing-schedule annotation is withdrawn

Added in 26.07.010, removed here — it did not work in practice. Removed from all
three modules: the sidebar controls, the message outputs, the four plot call
sites (growth, both body-weight trajectory branches, and the KM curve) and the
Info-tab documentation. `tumorGrowthUI` is back to its exact pre-feature size,
which is the cheap check that the removal was complete.

The schedule CSVs stay removed and the demo-reachability exemption list stays
empty: every file in `inst/sample_data/` is loadable from the selector.

**What this says about the dashboard tests.** The feature had 102 assertions
here, all passing. They verified that the controls existed in the UI source, that
the toggle gated the annotation, that both toxicity branches were wired, that
parse warnings reached a banner, and that failure degraded to the unannotated
plot. Every one of those was about wiring — whether a call site existed, whether
a layer count changed — and none about whether the result was usable.

This is §I.2 in a new costume. That finding was that the module tests asserted
`inherits(session, "session_proxy")` and so verified the harness rather than the
module. These assertions verified the plumbing rather than the picture. Counting
call sites did catch a real half-application (only one of the two toxicity
branches was wired at first), so the approach has value — it just cannot
substitute for looking at the output in the running app.

| ID | Issue | Severity | Status |
|---|---|---|---|
| D.16 | Dosing-schedule annotation non-functional in practice | — | ✅ Removed v26.07.011 |

---

# Review Round 14 (v26.07.012 — 2026-07-30)

## R14.3 The Tumor Growth Info tab claimed TGI it does not produce — Major

*"All modes compute Tumor Growth Inhibition (TGI) for each treatment group
relative to the reference."* `tumor_growth_statistics()` returns no TGI field and
the tab renders none.

The claim predates this review — it read "Both modes…" before Round 11 widened it
to "All" — and Round 11 propagated it without checking. That is precisely the
failure Round 11 was convened to fix, which is worth stating plainly rather than
quietly correcting: an audit that rewrites a claim without verifying it has
laundered the error, not caught it.

Fixed, and the replacement explains where TGI lives and why a trajectory model
does not report one. Pinned by `test-docs-current.R`.

## R14.5 The Synergy tab shows the more biased of two verdicts — Major

See the backend's Round 14 entry for the numbers. Summary: the linear-Loewe index
calls a Bliss-additive combination "antagonistic" across roughly FE 0.25–0.65,
and the tab renders only that index. The backend's `overall_assessment` — which
requires Bliss *and* Loewe to agree before claiming synergy — is returned and
never displayed.

Left to the owner: headlining the conservative label, showing both, or keeping
the current behaviour with a caveat are all defensible, and the choice is
scientific rather than technical.

| ID | Issue | Severity | Status |
|---|---|---|---|
| R14.3 | Info tab claims TGI the tab does not produce | Major | ✅ Fixed v26.07.012 |
| R14.5 | Conservative synergy label computed but never shown | Major | Open — owner's call |

---

# Review Round 15 (v26.07.013 — 2026-07-30)

## R14.6 The forest plot used unadjusted p-values on the AUC path — **Critical**

`mod_tumor_growth_plots.R` sniffed
`intersect(c("p.value", "p_value", "p.adj", "adj.p.value", "P_Value"), ...)` and
never listed `p_adjusted`. The AUC path's `p_value` is the raw t-test result, so
the Treatment Effects forest plot marked significance with no multiplicity
correction, while the GAM path — where `p_value` holds the adjusted value — used
corrected ones. Same figure, same legend, two multiplicity regimes, switching on a
model selection the reader cannot see from the plot.

This is the §K.2 sniff pattern again, and the sixth finding traceable to it. What
makes this instance the worst so far is that the sniff *succeeded*: it found a
column, so nothing was blank or missing. It simply found the wrong one.

Fixed by preferring the backend's new `P_Value_Adjusted`.

## R14.5 The consensus synergy verdict is now displayed — **Fixed**

Backend `overall_assessment` requires both models to agree before claiming
synergy, and was never rendered. Now shown above the CI verdict, with an explicit
disagreement warning when Bliss and Loewe point opposite ways.

## R14.7 Bayesian intervals were labelled median/HPD; they are mean/equal-tailed

See the backend's Round 15 entry for the verification against brms. Both interval
types appear under `Lower_95_CrI`/`Upper_95_CrI` — equal-tailed from
`summary(brmsfit)`, highest-density from emmeans — and the tab now distinguishes
them.

| ID | Issue | Severity | Status |
|---|---|---|---|
| R14.6 | Forest plot used unadjusted p-values on the AUC path | **Critical** | ✅ Fixed v26.07.013 |
| R14.5 | Consensus synergy verdict never shown | Major | ✅ Fixed v26.07.013 |
| R14.7 | Bayesian intervals mislabelled | Major | ✅ Fixed v26.07.013 |

---

# Review Round 18 (v26.07.014 — 2026-07-31)

## R17.2 The Combination Index is gone from the Synergy tab

Backend v0.21.0 removed it — see that repo's Round 18 entry for the measurement
(42 % of the effect grid disagreed with Bliss, always in the antagonistic
direction, because the implementation was response additivity rather than Loewe
additivity).

Removed here: the Loewe Info panel, the CI plot and its download handler, the CI
entry in the plot selector, the shaded-region checkbox that only that plot used,
the `Loewe Difference` / `Combination Index` metric descriptions, the
`Loewe_Expected_TGI` / `Combination_Index` over-time column notes, and the
`current_syn_ci_plot` reactive.

The verdict box was reworked. It previously described itself as a *consensus*
requiring both models to agree, with a warning when they disagreed. With one
model there is nothing to reconcile, so it now names Bliss explicitly, states the
formula, and says plainly that a threshold applied to three group means should be
read against the interval below it rather than on its own.

The Help entry became an explanation of *why* there is no CI rather than a
description of one, including the sham-combination argument. A user who expects a
Combination Index and finds none deserves the reason.

## R17.3 Three calls to backend functions that do not exist — **Major**

| call | why it failed |
|---|---|
| `mouseExperiment::plot_survival` | removed from the backend in `ace5b94` |
| `mouseExperiment::prepare_dose_data` | internal, never exported |
| `mouseExperiment::generate_summary_statistics` | internal, never exported |

All three sat inside `tryCatch`, so every survival and dose-response run threw,
logged, and fell back. Output appeared, so nothing looked broken.

**Dose-response:** `result$analysis_data` was never set, three renderers hit a
silent `req(ad)`, and the scatter, box and trend plots all came out blank.
**Survival:** the comment said *"Use the exported plot_survival() function
(matches demo notebook approach)"* — that plot has never been drawn; every user
has seen the fallback.

Fixed by reading `analysis_data` from the backend result (it now returns the
frame it analysed) and building the KM plot directly here.

## R17.4 A permanent guard against cross-repo staleness

`test-backend-calls-resolve.R` scans every `mouseExperiment::` call in `R/` and
fails if any does not resolve to an exported function. It strips comments first —
comments naming a removed function are wanted, since they explain what replaced
it; only live calls matter. That distinction is not hypothetical: an earlier
guard written without it produced false positives on its own documentation.

It also asserts the matcher can fail, because a guard that cannot fail is not a
guard, and that dose-response results carry `analysis_data` so the blank-plot
failure cannot return.

Same class as §R4.D2, where the suite passed against a six-release-stale backend.
Four lines, and it would have caught all three of these the day they broke.

| ID | Issue | Severity | Status |
|---|---|---|---|
| R17.2 | Combination Index removed from the UI | Major | ✅ Done v26.07.014 |
| R17.3 | Three dead cross-repo calls; blank dose-response plots | Major | ✅ Fixed v26.07.014 |
| R17.4 | No guard against cross-repo staleness | Process | ✅ Added — 11 assertions |

Suite: 558 → 569 passing, no skips.

---

# Review Round 20 (v26.07.015 — 2026-08-05) — production incident

## R19.6 Stale NAMESPACE imports broke the build; site returned 503

**Symptom:** after pushing, mouseexperiment.com returned 503 while the VPS was up.

**Cause:** v0.21.0 removed `plot_combination_index()` and
`plot_synergy_combined()` from the backend. Two `importFrom()` lines for them
remained in this package's NAMESPACE. `R CMD INSTALL` resolves imports at
lazy-load and fails:

```
Error: object 'plot_combination_index' is not exported by 'namespace:mouseExperiment'
ERROR: lazy loading failed for package 'mouseExperimentDashboard'
```

The Dockerfile's `remotes::install_github('sciOmics/mouseExperimentDashboard@staging')`
step therefore died, no image was produced, and NPM had no healthy upstream.

**Why nothing caught it — two compounding gaps, both mine.**

*First*, the §R17.4 guard I added specifically to prevent cross-repo staleness
scanned `mouseExperiment::` **call syntax**. `importFrom()` in NAMESPACE is a
different mechanism entirely, and the guard never looked at it. I wrote a check
for one half of the problem and treated it as covering the whole.

*Second, and more important:* **`devtools::load_all()` does not resolve NAMESPACE
imports the way installation does.** Every suite run in this entire review used
`load_all()`. A suite can be — and was — completely green against a package that
cannot be installed. 569 passing tests, and the artifact they were testing would
not build.

**Fix:** removed the two dead `importFrom()` lines. Verified by the path that
actually matters — `R CMD INSTALL` succeeds, `library()` loads, and
`runDashboard()` serves **HTTP 200** on a real port. That last check had never
been run in this review either; the container's entry point was never exercised.

**Guard:** `test-backend-calls-resolve.R` now also scans every
`importFrom(mouseExperiment, ...)` in NAMESPACE against the backend's exports,
with a companion test proving the scan can fail.

**The standing lesson.** Green tests are not a deployable artifact. Any release
that changes exports on either side of this pair needs `R CMD INSTALL` and a boot
check, not just `testthat`. Both are fast; neither was being run.

| ID | Issue | Severity | Status |
|---|---|---|---|
| R19.6 | Stale NAMESPACE imports broke install → 503 | **Critical** | ✅ Fixed v26.07.015 |
| — | `load_all()` never validates imports; suite green on unbuildable code | Process | ✅ Guarded + install/boot check documented |


---

# Review Round 21 (v26.07.015 — 2026-09-29)

**Scope.** Every file in `R/`, `inst/www/js`, the docs and the test suite, reviewed by two reviewers (analysis modules; infrastructure and reliability) plus a lead. Cross-repo findings came from the package's Round 20 and the VPS repo's Round 3. Evidence was executed with `shiny::testServer()`, with shinytest2 in headless Chrome (installed for this round), and by calling helpers directly. The lead re-ran or read every Critical.

**Baseline at HEAD (`b8a58c3`).**

| Check | Result |
|---|---|
| `devtools::test()` | 574 expectations green, but see R20.D11: `test-shiny_app.R` has never executed |
| `R CMD INSTALL` + `runDashboard()` on a real port | **HTTP 200**, 470 KB initial page, all six JS assets referenced (the Round 20 install-and-boot standard, re-run) |
| `R CMD check` | 6 WARNINGs, 5 NOTEs (R20.D26). The NOTEs hide three live defects (R20.D3, R20.D4, R20.D19) |

**How this relates to the package.** Most statistical Criticals surface here, because the dashboard calls the affected functions with their defaults. See package §R20-A (themes T1–T4) and R20.1–R20.13. This round records only what the dashboard itself does.

## R20-D.A Critical

### R20.D1 DD/MM/YYYY dates are silently read as MM/DD/YYYY, so every time-based result is wrong — **Critical**

`R/helpers_data.R:93` (`detect_date_format`):

```r
if (grepl("^\\d{1,2}/\\d{1,2}/\\d{4}$", first_val)) { fmt <- "%m/%d/%Y" ...
```

The format is chosen from the **first row only**, and any `d/d/yyyy` value is assumed to be US month-first. A European export (DD/MM/YYYY) matches, `fmt` is non-NULL, and no warning fires. The item 1.4 "date_fmt_warning" fix fires only when the detector returns NULL, which never happens for this input.

Lead re-run. The study starts 1 March 2025, with measurements on 01/03, 04/03, 08/03, 11/03, 15/03 and 18/03:

```
true days        0   3   7   10   14   17
dashboard Day  -57  33 155  247   NA   NA
```

Days 1–12 are read as months and become silently wrong values; days 13–31 become NA and are dropped with only a generic "Some dates could not be parsed" toast. Every growth rate, TGI and KM time downstream is wrong.

A related problem (lead, read-only): year-less formats (`15-Mar`, `March 15` plus the reference year) assume a single calendar year, so a study that crosses New Year gets negative days.

**Fix:** don't accept an ambiguous `d/d/yyyy` layout on faith. Scan the whole column: any first component above 12 forces day-first, any second component above 12 forces month-first. When the order is still ambiguous, ask the user. Warn whenever `calculate_dates()` produces NA or negative days relative to the reference date.

### R20.D2 The survival per-animal reduction keys animals by ID (plus cage) only, so IDs reused across arms are merged — **Critical**

`R/mod_survival.R:649-660` (column mode) and `:730-735` (volume mode); the Bayesian path uses the same frame.

```r
key_cols <- if (has_cage) c(id_col_surv, cage_col_surv) else id_col_surv
mask <- data[[id_col_surv]] == unique_animals[i, id_col_surv]
```

Ear tags reused per group are common. With no cage mapped, or cages numbered per arm, all of an ID's rows across arms become one record. That record takes the latest day, the maximum event across all arms, and the treatment of whichever row sorts last. The backend's one-row-per-animal check cannot detect this. Executed with the loop copied verbatim:

| Case | Dashboard | Correct reduction |
|---|---|---|
| 32 animals | 8 records, all "Combo" → error `Reference group 'Vehicle' is not present` | — |
| A silent case | 10 records (1 Vehicle, 9 DrugA); HR 0.197, p = 0.19; DrugA events 9/9 | 20 records; HR 0.123, p = 0.002; DrugA events 7/10 |

A related problem: in column mode the reduction takes the **last day** even when rows continue after a flagged event. An event on day 20 with rows to day 28 gives Day 28, Event 1 (executed; `:666-677`).

**Fix:** key by treatment + ID (+ cage). Better: move the reduction into the backend using `make_mouse_key()` (package T1), and use the first event day.

### R20.D34 On the bundled Weight demo with default settings, the toxicity tumour-mass correction is 1000× too large, and the warning is swallowed — **Critical**

`R/mod_toxicity.R:1316-1347,1377-1413` makes seven `mouseExperiment::` calls; none passes `volume_units`, and the UI has no unit selector. `R/mod_toxicity.R:1357-1358` re-derives net weight itself with a hard-coded `/1000`. The sidebar label (`:51`) says "Volume (mm³)". `R/helpers_ui.R:105-123`: `safe_analysis()` has only an `error =` handler, so warnings are dropped.

The package's unit heuristic classes the demo's genuinely-mm³ volumes (median 6.21) as cm³ (package R20.83, lead re-run). Three backend functions therefore subtract tumour masses 1000× too large, while four others and the dashboard's own trajectory plot assume mm³. One tab ends up showing mutually inconsistent, impossible numbers under "Analysis completed successfully!". Executed:

| Output | Shown | With `volume_units = "mm3"` |
|---|---|---|
| Net weight, LR/DMSO day 14 (22.2 g mouse, 1,092.6 mm³ tumour) | −1,070.4 g in the model; 21.1 g in the trajectory plot | 21.1 g |
| Body-weight EMMs (g) | DMSO −173.4, Drug A −166.1, Drug A+B 17.1, Drug B −18.7 | positive |
| ≥ 20 % weight-loss events | 20 of 20 animals | 0 of 20 |
| Top TWM arm | Drug A+B (2.13) | Drug B (6.44) |

The backend's own warning ("Estimated tumour mass exceeds 50 % of body weight (max 5309 %) … off by 1000x") never reaches the user.

**Fix:**
- Add a mm³/cm³ selector to Data Upload and pass `volume_units` to all seven calls (the package adds the argument to four functions under R20.27/R20.83).
- Plot the backend's `Net_Weight` rather than re-deriving it.
- Make `safe_analysis()` collect warnings with `withCallingHandlers()` and show them in the status box. That last change alone would have exposed this.

### R20.D35 Synergy and TWM "95 % intervals" are shown as uncertainty for estimates they do not describe — **Critical** (root cause package R20.2)

`R/mod_synergy.R:894-962` renders the CI table and an "interval excludes zero" verdict; `R/mod_toxicity.R:1584-1620` renders the TWM table and the ranking note. Under the default "Model-based" estimand, the backend's intervals come from last-observation arithmetic means (package R20.2).

- **Master demo, day 20:** TGI_A is 90.2 % while its "95 % CI" is [72.5, 85.6]. The Bliss excess is 0.018 ("Synergy"), yet its interval is [0.054, 0.171], and the tab says *"The Bliss-excess 95 % interval [0.054, 0.171] excludes zero, consistent with … synergistic direction."*
- **TWM table:** every treated arm has a point estimate outside its own bounds (Drug_A Low: TGI 87.97 %, CI [−1.34, 38.33]; TWM 6.43, CI [0, 2.93]). The ranking note orders arms by model TWM but separates them with last-observation bounds.
- **Control case:** with `last_obs` selected, the point sits inside the interval (79.6 within [72.5, 85.6]).

Executed.

**Interim dashboard fix:** until the package bootstraps the reported estimand, suppress the interval columns and the verdict note when `endpoint_method == "model"`, or relabel them as "interval for the last-observation estimand".

## R20-D.B Major

### R20.D3 All six Bayesian "legacy inline" fallbacks crash; the "fit running" modal then stays open forever — **Major**

`R/tg_resume.R:142,145`; `surv_resume.R:130,133`; `bw_resume.R:109,111`; `dr_resume.R:107,110`; `syn_resume.R:117,120`; `twm_resume.R:89,92`.

`%...>%` and `%...!%` are not imported (`exists("%...>%", asNamespace("mouseExperimentDashboard"))` is FALSE, and R CMD check flags it). `a %...>% b %...!% c` parses as `` `%...!%`(`%...>%`(a, b), c) ``, so R fails to find `%...!%` before the promise is ever built, and **no fit starts**. `show_analysis_modal()` has already run, and `removeModal()` lives only in the finalisers, so the session is wedged until reload.

The path is taken whenever `get_sid_or_none()` returns `SID_NONE`, i.e. no `mE_session` cookie on the websocket request. A normal first visit does get the cookie (verified in headless Chrome). Browsers that block or clear cookies do not. Executed:
- Each `.legacy_inline_*` called directly gives `could not find function "%...!%"`, and the stubbed promise builder was called 0 times.
- In headless Chrome with `Emulation.setDocumentCookieDisabled(TRUE)`, the "Bayesian tumor growth fit running…" modal was still open 12 s later, with `Error in %...!%` in the R log.

**Fix:** delete the legacy path and give `SID_NONE` sessions a server-generated `make_sid()`, so they use the registry like everyone else. Failing that, use `promises::then()` as `launch_job()` already does, or `@importFrom promises %...>% %...!%`.

### R20.D4 `NA_Date_` is undefined, so a non-US reference date crashes date preprocessing — **Major**

`R/helpers_data.R:89`: `error = function(e) NA_Date_`. That symbol is not visible from the namespace. The handler fires when the *reference date* matches none of the five `tryFormats`, e.g. `24/03/2025`, `2025/03/24`, `24.03.2025`, `24 March 2025` or `20250324`. It then raises `object 'NA_Date_' not found` instead of falling back to `Sys.Date()`. What the user sees:
- **Tumour Growth** (`mod_tumor_growth.R:90`, no surrounding `tryCatch`): the run observer errors, and "Running analysis…" stays up.
- **Survival** (`mod_survival.R:614`, a bare `ensure_preprocessed(rv)`; every other module wraps it): the raw message.
- **Data tab:** the error is swallowed and `rv$processed_data` stays NULL.

**Fix:** `as.Date(NA)`. The next line already handles NA.

### R20.D5 "Dismiss all" from a second tab strands the first tab's modal and discards its completed result — **Major**

`R/bayes_job_runner.R:43-46`, `R/resume_banner.R:107-112`, `R/job_registry.R:110-112`.

1. The banner renders in every tab that shares the cookie. **Dismiss all** retires every job, including running ones.
2. The first tab's poller then finds `get_job()` NULL, clears its job id and returns **without calling a finaliser**, so the modal stays open.
3. When the fit resolves, `onFulfilled` finds the job retired and drops the result.

The modal backdrop is click-through (`overlay.js:33-35`, `custom.css:397`), so this can also happen within a single tab.

Executed:
- **Browser (two tabs):** tab 1 stayed `modal_open: 1`, "Running analysis…", `jobs: <none>`, through +17.5 s.
- **Unit test:** retiring a job under a live poller fires no finaliser.

The banner's **Cancel** is handled correctly.

**Fix:** in the poller's `is.null(job)` branch, dismiss the modal and clear the running flag. Also let Dismiss clear terminal jobs only.

### R20.D6 Long analyses run synchronously on the one R process and freeze every user, and the 2-minute power simulation is reachable without an upload — **Major**

`future::plan(multisession, workers = 1L)` (`R/app.R:289`) serves only the six Bayesian promises. Everything else blocks the single Shiny process. Timings on the master demo with module defaults (executed):

| Analysis | Blocks every session for |
|---|---|
| A-priori LMM power simulation, defaults (500 sims × 5 N × 2 α) | **121 s**, scaling to the UI maximum of 2,000 sims |
| TG lme4 + permutation test, 999 / 9,999 permutations | 22 s / **145 s** |
| TG lme4 with diagnostics (package R20.76) | ~8 s here, **20–140 s** on realistic sizes |
| `analyze_body_weight` | 15 s |
| GAMM, AUC, survival, synergy (2,000 boot), TWM, dose-response | < 1.5 s |

The Power tab is enabled at startup with no data (`R/app.R:228`), so any anonymous visitor can freeze the server for about two minutes. Meanwhile every other session shows the disconnect overlay once the proxy timeout passes.

**Fix:** route the power simulation, permutation test and body-weight fit through the same `future` job path. Fix R20.76 in the package (per-animal influence). Consider a second worker with a global cap.

### R20.D7 `cancel_job()` does not stop the worker, one worker serialises every user's fit, and there is no per-browser cap — **Major**

`R/job_registry.R:174-181` only sets `status = "cancelled"`; its docstring says the future "continues to completion". With `workers = 1L`, a second user's fit queues behind the first for up to 12 minutes. A user who cancels and re-runs waits for both. Repeated clicks enqueue multiple full fits.

Read-only: this follows from `future_promise` queueing semantics and the implementation.

**Fix:** tell users in the UI that only one fit runs at a time. Add a per-sid concurrent-job cap. Use an interruptible backend (e.g. `callr`-based) or more workers under a global cap.

### R20.D8 Bayesian survival from volume-derived events always fails — **Major**

`R/mod_survival.R:746-750` creates `.__DerivedTime` and `.__DerivedEvent`. brms rejects names with double underscores (`Variable names may not contain double underscores or underscores at the end`), so the "derive from tumour volume threshold" + Bayesian combination errors on every input (executed).

**Fix:** use syntactic names without double underscores. In the package, map user columns to internal names before building any brms formula (package R20.43).

### R20.D9 `.rda` upload does not work, and it loads into the process-wide global environment — **Major**

`R/helpers_data.R:239-242`: `load(file_path, envir = .GlobalEnv); data <- get(ls()[1])`. Two problems:
1. **The feature is dead.** `ls()` inside the function lists the function's own frame, so `get(ls()[1])` returns a local variable, and every valid `.rda` is rejected as "Loaded file is not a data frame" (executed).
2. **Loaded objects persist and are shared.** The file's objects are written into `.GlobalEnv`, which every Shiny session in the container shares. They persist across sessions and can shadow functions that other code calls unqualified (`validate_uploaded_file()` calls `read.csv`; `logging.R:178` calls `tail`). Uploaded files must not be able to change the process-global environment.

**Fix:** load into a private environment, accept only data frames, and never touch `.GlobalEnv`. If `.rda` is not a feature anyone needs, remove it from `accept =` (`ui_modules.R:130`); that is the simpler fix. Recorded at design level; no demonstration is reproduced here.

### R20.D10 Session restore reads an uploaded `.rds` and copies every key into app state without validation — **Major (hardening)**

`R/mod_log_session.R:104-117`:

```r
snapshot <- readRDS(file_info$datapath)
...
for (key in names(snapshot)) { if (key != "saved_at") rv[[key]] <- snapshot[[key]] }
```

An uploaded file is untrusted input. Here it is deserialised with `readRDS()`, and every top-level key is written into `rv` with no key allow-list, no type or shape checks, and no check that the restored column mappings (`settings$col_*`, `tx_col`, `vol_col`, …) name columns that exist in the restored data. The CSV upload path is safe on that last point, because `read.csv(check.names = TRUE)` sanitises headers. Restore bypasses that, and the package pastes column names into model formulas (package R20.43).

The JSON **preset** path in the same file (`:153-196`) already does this correctly: it parses JSON and applies each setting through `updateSelectInput()`, which constrains values to the offered choices.

**Fix:** replace the `.rds` snapshot with a data-only format (the uploaded data as CSV plus settings as JSON, modelled on the preset path). Restore only an allow-listed set of keys with type checks, validate every `col_*` against `names(data)`, and cap the file size. Recorded at design level.

### R20.D11 `test-shiny_app.R` has never executed, and the only browser-level coverage is dead — **Major (process)**

`tests/testthat/test-shiny_app.R:16`: `skip_if_not_installed("chromote", message = "…")`. testthat 3.3.2's `skip_if_not_installed()` has no `message` argument.
- **Without shinytest2** (every previous run in this environment), the file skipped at line 14, at file level, so it never even appeared in the skip count.
- **With shinytest2 and chromote installed**, it **errors** at line 16 before any `test_that()`. The full suite reports 574 passing and 1 error, and zero tests run from this file.

**Fix:** drop the `message =` argument (use `skip_if()` for the Chrome check), install shinytest2/chromote in CI (package R20.53), and make a file-level error fail the run.

### R20.D12 Treatment names, reference group and column names are sent to the analytics server — **Major (confidentiality)**

`R/mod_tumor_growth.R:72-74` and `R/mod_survival.R:604-607` (`trackAnalysisStart`/`trackReferenceGroupSelected` carry the reference or treatment label, e.g. `"Drug_A Mid + Drug_B"`); `R/mod_data_upload.R:655-662` (`trackColumnMapping` carries the user's column names); `inst/www/js/umami-events.js:60-95`. Error payloads (`trackError`, up to 200 characters) can include a column name.

The endpoint is first-party (self-hosted Umami), but unpublished compound and study identifiers end up in an analytics database with no notice. File names are **not** sent.

Staging traffic is also recorded under the production website id, because the Umami tag in `app.R:127-131` is unconditional.

**Fix:** send only non-identifying facts (analysis type, method, counts), and gate analytics on consent or a notice. Use a separate website id for staging.

### R20.D13 The "Treat natural death as event (recommended)" checkbox is never read — **Major**

`R/mod_survival.R:76-87` is the only reference to the input; the Info text (`:205-209`) says it matters. The event derivation (`:717-754`) never reads it, and the results are identical with the box checked or unchecked.

Executed with a constructed natural death (animal 17 removed at day 10 at 248 mm³):
- **Volume mode:** Event = 0 even with the box checked ("recommended"), so toxic deaths in treated arms are censored, which favours the treatment.
- **Column mode:** Event = 1 even with it unchecked.
- **Consequence:** the Drug_A Mid HR is 0.0074 or 0.0389 depending only on the event source.

**Fix:** implement the option or remove it.

### R20.D14 Survival and default over-time synergy fail unless the ID column is literally named `ID` — **Major**

`mod_survival.R:805-815` and `mod_synergy.R:692-705` never pass `id_column`, while the column detector maps `mouse_id`, `animal_id`, `ear tag` and similar. Every bundled demo happens to use `ID`, which is why no test catches this. Executed on the master data with the ID column renamed "Animal":
- **Frequentist survival** fails with `undefined columns selected`. The Bayesian survival path does pass `id_column`, and the direct backend call with `id_column = "Animal"` succeeds.
- **Over-time synergy**, which runs on the default "All dates" time point, fails with the misleading "Could not calculate synergy for any time points". The same data works with a single day selected.

"All dates" also silently ignores the Endpoint-estimand and Bootstrap inputs: the over-time function accepts neither argument, and no interval table or attrition table is shown, contrary to the sidebar help. Internally it still runs 2,000 bootstrap resamples at each of 11 days and throws them away (package R20.15, R20.81). Every survival run also prints the full HR report to server stdout, because `verbose = FALSE` is not passed.

**Fix:** pass `id_column` (and cage) from both calls. The package gains the matching arguments in R20.15.

### R20.D15 Every table in the downloadable HTML report has all its column names in one header cell — **Major**

`R/helpers_report.R:44` applies `esc()` to `colnames(df)`, but `esc()` collapses vectors (`paste(xc, collapse = ", ")`). The header row is therefore a single `<th>` reading "Treatment, Estimate, SE", while each row has one `<td>` per column. Columns 2 onward are unlabelled in every table of every report, which is the artefact users send to colleagues.

Lead, executed: `th count: 1, td per row: 3`. Cells, captions, settings and the title *are* HTML-escaped, which was verified.

**Fix:** `vapply(colnames(df), esc, "")`.

### R20.D16 The combined efficacy–toxicity table and its download always fail — **Major**

`R/mod_toxicity.R:1741,1745`; `R/mod_toxicity_downloads.R:171,175`. The code selects `benefit_table[, c("Treatment","Benefit")]` and `comparison[, c("Treatment","Corrected_TGI","N_Excluded")]`, but the backend returns `Benefit_Score` and `TGI_Corrected`. The result is `undefined columns selected` whenever TWM and TBA succeed. `N_Excluded` lives in `excluded_summary`, not in `comparison`. Executed: every Efficacy–Toxicity run shows the error in the Combined Report tab, and `download_et_report` fails. `mod_toxicity.R:1479` `bt$Benefit` works only through `$` partial matching.

**Fix:** use the backend's names, and add these columns to `test-backend-contract.R`.

### R20.D17 The Power tab's analytic path ignores multiplicity and attrition — **Major**

`R/mod_power.R:494-530` calls `apriori_power_analysis()` without `p_adjust_method`, `n_comparisons` or `dropout_rate`, so the R3.15 and K.15 fixes are unreachable from the UI (lead, read). Combined with package R20.7 (√k), sample sizes the Power tab recommends for k ≥ 3 arms are too small twice over.

**Fix:** expose all three, defaulting to the k − 1 vs-control family.

### R20.D18 `result$model <- NULL; gc()` frees nothing, so every Bayesian result keeps its ~31 MB fit alive — **Major (memory)**

`tg_resume.R:100`, `bw_resume.R:69`, `surv_resume.R:84`, `syn_resume.R:72`, `dr_resume.R:68`. The result's ggplot objects hold the model in their `plot_env`, so the list still serialises to 30.8 MB without `$model` (0.01 MB once the plots are removed; package R20.80). `tg_results_prev` keeps an older fit alive as well. On a container limited to 4 GiB (1 GiB on staging, VPS V3.12), a handful of concurrent Bayesian users exhausts memory.

**Fix:** store plots that don't capture the model (built from data frames, or pre-rendered), then drop the model.

### R20.D19 The survival forest-plot download calls a function outside its scope — **Major**

`R/mod_survival_downloads.R:84` calls `build_surv_forest_plot()`, which is defined as a local inside `survival_server`'s module closure (`R/mod_survival.R:1346`). `add_surv_download_handlers()` is a separate top-level function whose enclosing environment is the package namespace, where the symbol does not exist (`exists(..., asNamespace(...))` is FALSE, and R CMD check flags it). Executed under `testServer`: the on-screen forest plot renders, but "Download Plot" on it returns HTTP 500 (`could not find function "build_surv_forest_plot"`), and the user gets a failed browser download with no in-app message. The KM, Cox CSV and risk-table downloads work. A `codetools` scan of every namespace function found no other closure-captured helper, apart from the promise operators in R20.D3.

**Fix:** move `build_surv_forest_plot()` to top level (`helpers_plot.R`), or pass it in as an argument, as the C.1/C.3 refactor convention prescribes for shared reactives.

### R20.D52 "Download All" zips a directory shared by every session, so one user's download can contain another user's files — **Major (confidentiality)**

`R/mod_tumor_growth_downloads.R:208-247`:

```r
outdir <- file.path(tempdir(), "mouseExperiment_results")
...
files <- list.files(outdir, full.names = TRUE)
zip(file, files, flags = "-j")
unlink(outdir, recursive = TRUE)
```

`tempdir()` is per process, and every user shares the one R process. `unlink()` runs only after a successful export (there is no `on.exit`), so any file left by another session's failed export is zipped into the next user's download. Executed: `OTHER_SESSION_growth_rates.csv` appeared in the ZIP. The Toxicity ZIP already uses a per-call directory correctly.

**Fix:** `outdir <- tempfile("results_")` per call, with `on.exit(unlink(outdir, recursive = TRUE), add = TRUE)`.

### Findings in the analysis modules

#### R20.D36 Synergy always fails after a Tumor Growth run when a dose column is mapped — **Major**

The Synergy selectors are built from the raw treatment column (`R/mod_synergy.R:512-537`). The analysis, however, uses `pp$tx_col` (`:627-630`), which the TG run has switched to "TreatmentGroup" (e.g. "Drug_A Mid 15"; `R/mod_tumor_growth.R:150-155`). The Synergy Info tab tells users to run TG first, and the Master demo auto-maps Dose, so the demo always ends in "The following specified groups do not exist in the treatment column" (executed).

**Fix:** build the choices from `ensure_preprocessed(rv)$data[[tx_col]]`, or pass `rv$settings$col_treatment` as `treatment_column`.

#### R20.D37 The Synergy "Tumor Growth Inhibition" plot recomputes TGI from surviving animals and contradicts the table beside it — **Major**

`R/mod_synergy.R:1175-1196`: `tgi_vals <- 100 * (1 - g_rows[[vol_col]] / ctrl_mean)`, where `ctrl_mean` is taken over the animals present that day. That is the "survivors" estimand the tab itself warns against; at day 20 one of eight Vehicle mice remains. This is the re-derivation pattern recorded in the maintainer notes (the dashboard computing a number the package already returns). Executed:

| Day | Arm | Backend table / trend plot | TGI plot |
|---|---|---|---|
| 20 | Drug_B | 76.0 | 55.7 |
| 20 | Drug_A | 90.2 | 83.5 |
| 18 | Drug_B | 72.0 | 56.8 |

**Fix:** plot the backend's `synergy_summary` TGIs, or state the estimand on the plot.

#### R20.D38 The Synergy Info tab documents the discredited response-additivity formula — **Major**

`R/mod_synergy.R:107-110,124,133-136,204` says *"Expected TGI = min(FE_A + FE_B, 1.0) × 100"* and "sum of individual TGIs, capped at 100 %". That is the null R17.2 removed as biased toward antagonism. The backend computes Bliss, `FE_A + FE_B − FE_A·FE_B`, and the verdict box on the same tab states it correctly. The Info tab's "ceiling" warning is wrong as a result (50 % + 50 % gives an expected 75 %, not 100 %). It also says "two reference models are evaluated", lists one, and refers to a "CI threshold" (read-only).

**Fix:** rewrite the panel, and add a Bliss-formula check to `test-docs-current.R`.

#### R20.D39 The AUC model fails on first click with the default correction — **Major**

The default is Dunnett (`R/ui_modules.R:401-409`). The substitution (`R/mod_tumor_growth.R:285-299`) only swaps Dunnett and Tukey, and the backend rejects both on the AUC path (`utils_contrasts.R:88-92`). Executed:
- **With AUC selected and default settings:** "Error: p_adjust_method = 'dunnett' requires comparisons from a single fitted model…".
- **With "All Comparisons":** Dunnett is swapped for Tukey, which fails the same way.

**Fix:** for `model_type == "auc"`, map Dunnett/Tukey to Holm (or restrict the choices) and say so in the help text.

#### R20.D40 The GAMM effect-size forest plot crashes; its download fails and the ZIP silently omits it — **Major**

`R/mod_tumor_growth_plots.R:549-589`. GAM contrasts are per day: 25 rows with 5 repeated labels. `factor(Contrast, levels = …)` fails with "factor level [3] is duplicated", and the plot area shows the error. `download_tg_plot_effect_forest` errors, and `download_tg_all` drops `effect_forest.png` without saying so (executed).

**Fix:** facet by day (or show the last day only), with unique labels.

#### R20.D41 The dashboard replaces the backend's working GAMM residual plot with a broken one — **Major**

`R/mod_tumor_growth.R:575-596`: `rf_df <- data.frame(Fitted = stats::fitted(result$model), Residuals = stats::residuals(result$model))`. On gamm4 objects both return NULL, so the frame is 0 × 0, and it overwrites the backend's valid `diag_resid_fitted_plot`. The Diagnostics panel shows "object 'Fitted' not found". The backend plot prints fine on the same fit (executed).

**Fix:** keep the backend's `diag_*` plots and delete the re-derivation (`:540-596`).

#### R20.D42 Residual and random-effect Q-Q plots always look non-normal — **Major**

The displayed residual Q-Q is built in the dashboard (`R/mod_tumor_growth.R:544-556`: `qqnorm(resid)` + `geom_abline(0, 1)`), and the backend's RE and body-weight Q-Q plots have the same flaw (package R20.84). The help text (`ui_modules.R:1235-1245,1297-1305`) says "points should fall along the red dashed line", which cannot happen when the values are unstandardised. On the master demo (residual SD 0.163), 4.7 % of points lie within ±0.05 of the line for a correctly specified model (executed).

**Fix:** standardise the values, or use a `qqline` reference (the package fix covers the backend plots), and delete the dashboard copy.

#### R20.D43 Forest plots show unadjusted 95 % CIs next to multiplicity-adjusted p-values — **Major**

`R/mod_tumor_growth_plots.R:533-544` builds `est ± qt(.975, df)·SE`; `R/mod_tumor_growth.R:614-624` uses the raw Welch CI with `p_adjusted` on the AUC path; `R/mod_toxicity_plots.R:114-136` uses a raw Welch CI with BH colouring. Master demo, interim data, Dunnett: "Drug_A Mid [−0.555, −0.031] p = 0.111". The interval clears zero while the annotation says not significant (executed).

**Fix:** use adjusted intervals from the package (R20.85), or label the interval "unadjusted".

#### R20.D44 The Tumor Growth Comparisons tab mislabels the correction, discards its own title, and reads the family at render time — **Major**

- **Mislabelled correction.** The help (`R/mod_tumor_growth.R:1197`) says "p.value — adjusted p-value (Dunnett…)" even when the correction is None. The p.value is then 0.0341 unadjusted, where Dunnett gives 0.111.
- **Discarded title.** `comparisons_title()` computes a correct title at `:1406-1418` and it is never used. `output$tg_comparisons_title` (`:1374-1389`) has no UI slot, although the GAMM Info tab says "The Comparisons header states this explicitly".
- **Family read at render time.** After an "All Comparisons" run, switching the dropdown to "vs Reference" shows all 15 pairwise rows. After a vs-reference run, switching to "All Comparisons" shows "Comparisons not available for tukey", with no hint to re-run.
- **AUC help.** It always says "Bonferroni … all pairwise".

Executed.

**Fix:** render `comparisons_title()` as the table caption; read the family from `rv$tg_results$comparison_family`; generate the help text from `p_adjust_method_used`.

#### R20.D45 The TG Data Summary shows log-scale means labelled as mm³ — **Major**

The table is at `R/mod_tumor_growth.R:1723-1732`, and the tab text (`ui_modules.R:1035-1047`) says "Mean — mean tumor volume (mm³)". The backend's `data_summary` is on the model scale, even on the AUC path (package R20.57). Drug_A High day 0 is shown as 4.31 where the raw mean is 75.9 mm³; Vehicle day 20 is shown as 8.29 where the raw mean is 3,991.6. The same table reaches the HTML report and the CSV (executed).

**Fix:** back-transform (and call it a geometric mean), or label the scale from `meta$transform_used`. Ideally the package returns raw-scale summaries.

#### R20.D46 The TG HTML report records the widgets' current values, not the settings that produced the numbers — **Major**

`R/mod_tumor_growth.R:1629-1636` reads `input$tg_*` at download time. Executed: run LME4 / log / Vehicle / Holm, change the dropdowns without re-running, then download. The report says "model_type: auc, transform: sqrt, reference_group: Drug_B" above the LME4 / log / Vehicle tables. This is the same failure R3.D6 fixed for the on-screen metadata, now in the audit artefact.

**Fix:** take settings from `result$meta` and fields stored at run time.

#### R20.D47 Dose-Response pools different drugs by dose value — **Major** (dashboard side of package R20.16)

`R/mod_dose_response.R:716-725` passes the whole dataset, and only the dose column defines groups. On the Master demo:
- **Dose 10** is Drug_B, a different compound.
- **Dose 15** pools Drug_A Mid with the Drug_A Mid + Drug_B combination (n = 16, SD 793).

EC50 12.9, the Tukey table and the ANOVA all describe the mixture, and there is no control for choosing the agent (executed).

**Fix:** add an agent filter (or require one treatment plus control), and warn when one dose maps to more than one treatment.

#### R20.D48 Dose-Response "All dates" means each animal's last observation, and its TGI plot contradicts the Toxicity tab — **Major**

- **"All dates" is mislabelled.** The choice (`R/mod_dose_response.R:640-649`) maps to the backend's `time_point = NULL`, "last measurement for each mouse", so Vehicle animals are compared at days 16–20 and treated animals at days 22–26 (package R20.17).
- **The TGI plot re-derives TGI.** `:913-929` computes `(1 - vol / mean(control)) * 100` from the live `input$dr_control_group`. Its Drug_A Low mean is 18.3 %, while the Toxicity TWM table (model estimand, day 26) reports 88.0 % for the same arm, and the plot changes whenever the Reference dropdown moves, without a re-run.

Executed.

**Fix:** rename the choice ("Each animal's last observation"); state the estimand on the TGI plot or remove it (R14.3's own policy); store the control group at run time.

#### R20.D49 Dose-Response never shows the Jonckheere–Terpstra test or the EC50 CI, although the backend returns both — **Major**

`R/mod_dose_response.R:765-813,980-997`. The Info tab (`:192-206`) calls JT "the primary trend result", yet `res$trend_test$jonckheere_test` (re-enabled by package R3.31), `statistics$ec50_ci` ([11.28, 14.53]), `overall_trend_pvalue` and `aic_comparison_note` are never rendered or exported (executed; grep finds no `jonckheere` in server code). Same class as R14.5. Package R20.18 and R20.20 should be fixed before these are surfaced.

**Fix:** render them in the Summary tab and `download_dr_stats`.

#### R20.D50 The Toxicity Plots-tab "Download" never produces a file, for all eight plot types — **Major**

`R/mod_toxicity_downloads.R:81,95`: `p <- tryCatch(output[[plot_id]](), error = function(e) NULL)`. In Shiny 1.13, reading an output aborts ("Can't read output …") unless `shiny.allowoutputreads` is set, and even then the value is not a function. The user sees the toast "Plot is not available." while the plot is on screen. `"bw_forest"` is also missing from the switch (executed: all 7 tested selections fail).

**Fix:** hoist each plot into a `reactive()` shared by `renderPlot` and the handler, as the TG module does.

#### R20.D51 Power's "Seed SD from loaded data" pools volumes across all days, inflating the required N 5–40× — **Major**

`R/mod_power.R:415-419` computes one SD per group over every observation, so repeated measures and the growth trend are both counted as noise. Master demo: the seeded SD is 617.8, while within-day pooled SDs are 15.9 (day 0), 94.5 (day 10) and 265.0 (day 14). For Δ = 300 mm³ at 80 % power, that gives N = 68 per group, against 14 with the day-14 SD (executed). Combined with R20.D17 and package R20.7, the Power tab's recommendations can be wrong in either direction depending on the path taken.

**Fix:** seed from one user-chosen endpoint day, using each animal's last observation up to it.

### R20.D20 The synergy bar plot is all grey — **Minor**

`plot_drug_synergy()` maps `fill` to `Type` (Observed/Expected), but the dashboard passes colours named by group (`mod_synergy.R:1036-1046`), giving "No shared levels" and grey50 bars (executed).

**Fix:** name the palette by `Type`, or map `fill` to the group.

## R20-D.C Minor

| ID | Location | Issue | Fix |
|---|---|---|---|
| R20.D21 | `R/helpers_ui.R:61,65-67` | `icon("spinner fa-spin")` logs "does not correspond to a known icon" on every info alert and renders the regular-weight glyph | `icon("spinner", class = "fa-spin")` |
| R20.D22 | `R/app.R:215`, `mod_data_upload.R:127,167,655`, `mod_dose_response.R:747`, … | about ten `runjs(sprintf("…'%s'…", value))` calls interpolate values inside JS quotes instead of using the existing `js_string()` (`helpers_ui.R:124`). Effect is limited to the user's own browser: a group or column name containing an apostrophe breaks the call | use `js_string()` everywhere |
| R20.D23 | `R/helpers_data.R:235-281` | uploads: `shiny.maxRequestSize` is never set, so Shiny's 5 MB default applies while `DEPLOYMENT.md` says 100 MB; semicolon + decimal-comma CSVs are rejected with a raw error; Latin-1 and UTF-16 are rejected; tab-delimited `.csv` loads *silently* as one column; an unquoted comma in a value silently shifts columns. UTF-8 BOM works, and `.xlsx` is cleanly rejected (executed) | set the limit; sniff delimiter, decimal mark and encoding; actionable messages |
| R20.D24 | `global.R`, `app_standalone.R`, `shiny-server-app.R` | dead entry points: `global.R` sources `R/helpers.R` and `R/log_session_modules.R`, neither of which exists; none is in `.Rbuildignore` | delete, or build-ignore |
| R20.D25 | `R/logging.R:172-204` | every log line is teed into a process-global 500-line buffer that nothing reads (the Log tab was removed) | drop the appender |
| R20.D26 | `DESCRIPTION`, `NAMESPACE`, `.Rbuildignore` | R CMD check: `later`, `posterior`, `survival` used via `::` but undeclared; Imports `tidyr` and `wesanderson` unused; no `man/`, so `runDashboard` is undocumented; `\|>` in `mod_tumor_growth_plots.R` needs R ≥ 4.1 (`Depends: R (>= 4.0)`); non-ASCII in 22 files; ~100 "no visible global" NOTEs, most of them `stats`/`utils` functions resolved by the search path but three real (R20.D3, R20.D4, R20.D19). `mouseExperiment (>= 0.12.0)` is stale (the dashboard needs ≥ 0.21.0 for `analysis_data`); `Date` is 2026-06-18; `.Rbuildignore` misses `.waylog` (AI transcripts end up in built tarballs), `CODE_REVIEW.md`, `CHANGELOG.md`, `docs` | declare, trim, raise the floor to 0.22.0, build-ignore |
| R20.D27 | `docs/DEPLOYMENT.md`, `CONTRIBUTING.md` | `DEPLOYMENT.md` describes nonexistent Dockerfiles, a "pinned R so a CRAN bump can't break us" (the opposite of VPS V3.3), an HTTP-200 gate, Slack notification and SHA rollback that do not exist, `ssh root@…` (root login is disabled; port 2222; path `/opt/vps`), and `SHINY_FUTURE_WORKERS`/`SHINY_MAX_UPLOAD_MB` variables no code reads. `CONTRIBUTING.md` still says "backend ≥ 0.4.7" and "25-item review across Round 1 + Round 2", and its test table is incomplete | rewrite against VPS Round 3 |
| R20.D28 | `R/app.R:31-43,294` | the proxy-timeout reminder never fires in production (keyed on `SHINY_PORT`, which is never set; VPS V3.22); `launch.browser = TRUE` in a headless container | fire unconditionally; `launch.browser = FALSE` |
| R20.D29 | `mod_synergy.R:749,1126`; `mod_survival.R:848`; `mod_dose_response.R:661` | over-time plot `ylim(0, 100)` silently drops negative TGIs; the KM plot shows survminer's unadjusted omnibus log-rank p beside Bonferroni-adjusted Cox p-values; dose-response uses `rv$processed_data` rather than `ensure_preprocessed()`, so subset-filter changes after the first preprocess are not applied (read-only) | `coord_cartesian`; label the p; `ensure_preprocessed()` |
| R20.D30 | `mod_tumor_growth.R:1237-1258`; `ui_modules.R:659` | explainer text: the Wald χ² is called a "likelihood-ratio chi-square"; the LMM Type III text is shown for the GAM smooth table (package R20.60); the AUC ANOVA help describes `aov` with equal variances while the backend returns Welch; a crossed cage "enters as a random intercept" while the backend's `auto` uses a fixed effect | correct the text |
| R20.D31 | `mod_toxicity.R:1358,1378-1383` | trajectory plot hard-codes `/1000` (mm³, package R20.27); `cage_column` is not passed to `therapeutic_window_metric()` | resolve units; pass cage |
| R20.D32 | NUTS alert | tells users to "raise adapt_delta", which no argument allows (package enhancement) | expose `control` or change the advice |
| R20.D33 | `mod_data_upload.R:587,590` | `showFeedbackWarning`/`showFeedbackSuccess` with empty text spam the browser console with `Cannot read properties of undefined` | pass text or `hideFeedback()` |
| R20.D54 | `mod_survival.R:146,204-205,233-234` | survival help: says P_Value is "for the individual group coefficient" when it is Bonferroni-adjusted (0.000841 vs unadjusted 0.000168); "interval excludes 1.0 ⇒ significant" is false under adjustment; promises KM log-rank p-values (not shown on the Cox/Firth paths) and censor-column auto-detection (not implemented) | correct the text |
| R20.D55 | `mod_tumor_growth_plots.R:57-85` | the displayed growth-curve plot drops the necrotic × markers (47 on the master demo; the observer's correct plot is only a fallback), and "Extrapolation Points" has no effect without a mapped censor column | use the observer's plot; wire extrapolation |
| R20.D56 | GAMM path; Bayesian synergy | GAMM silently applies hidden LME4-only controls (transform `sqrt`, cage `never_include` were applied while hidden), although the Info tab says GAMM uses log; Bayesian synergy uses the hidden `syn_time_point` | reset or ignore hidden inputs |
| R20.D57 | `mod_tumor_growth.R:1580-1586`; `mod_tumor_growth_plots.R:178-180,565-572`; `ui_modules.R:995,1001` | drill-down merges cage-mates that share an ear tag (combo demo: 29 rows from two cages plotted as one mouse); the "geometric mean" bar clamps negative log-slopes to 1e-9 (−0.044 shown as 0.000); the all-pairs forest adds a spurious "Vehicle (Reference)" row; the growth-rate help says log1p where the backend uses `log()` | key by composite ID; plot the arithmetic mean of the slopes; fix labels |
| R20.D58 | `mod_toxicity.R:1568-1571`; `mod_toxicity_downloads.R:207-210` | BW AUC BH adjustment covers 5 vs-reference rows on screen but all 15 pairs in the ZIP (3.91e-03 vs 5.34e-03); the "Tumor Volume" covariate (checked by default) is silently dropped when "Subtract estimated tumor weight" is also on (the default) | one family; surface the backend message |
| R20.D59 | `helpers_bayes.R:58-94`; `tg_resume.R:94`; `mod_tumor_growth.R:1067-1092` | Bayesian panel shows a green "NUTS sampler diagnostics clean: NA divergence(s)" when the backend returns its failure default (the backend's `clean` flag is ignored); `output$tg_posterior_summary` has no UI slot; the finaliser overwrites `model_type_used` ("bayes_tg_gam" → "bayes"); the Posterior Check tab is inserted for "bayes" but not "bayes_gam"; Bayesian results never get a forest plot | read `clean`; add the slot; keep the field |
| R20.D60 | `mod_synergy.R:858-871`; `mod_power.R:665,705-716` | "Bliss Difference" is an unrounded fraction ("0.0179841522745299") on the single-day path but percentage points elsewhere; the Power summary reads a nonexistent `input$pw_ap_target_power` and always says "at target power = 80 %" (it said so with 0.95 selected); the method label is always "analytic (analytic)" | one unit; read the real input |
| R20.D61 | `ui_modules.R:692,716`; misc | Info tabs promise outputs that are never rendered (`treatment_effects_over_time`, "EMMs at five quantile days"); seven unused helpers; about 150 lines per TG run build plots that are never read; `detect_columns` maps the Necrotic demo's ID to "Tumor.Width..mm."; `translate_error`'s case-insensitive `"NaN\|NA\|missing values"` matches any "na" substring (latent); with `ragg` installed, drawing the RE Q-Q of a singular fit hangs the process (dev only; the production image has no ragg) | remove or render; tighten the regex |

## R20-D.D Efficiency

| ID | Where | Cost (executed) | Fix |
|---|---|---|---|
| R20.D53 | `mod_tumor_growth_plots.R:33-138` | every keystroke in the title, subtitle and axis-label boxes rebuilds the growth plot via `mouseExperiment::plot_tumor_growth()` and redraws it: typing a 33-character title caused 35 rebuilds and 19.5 s of rendering, all on the shared process | separate the data layer (keyed on data and extrapolation) from cosmetic layers; `debounce()` text inputs; `bindCache()` |
| — | `mod_tumor_growth.R:314` | every frequentist TG run requests `include_diagnostics = TRUE`: the case-deletion refits take 6.4 s of 6.6 s on the 48-mouse master demo (package R20.76), synchronously for every user (R20.D6) | compute on demand, or per animal |
| — | `mod_synergy.R:692-705` | over-time synergy runs 11 endpoint LMMs and 11 × 2,000 bootstrap resamples whose output is discarded (package R20.81) | fixed together with R20.D14 |

## R20-D.E Verified correct (recorded so it is not re-audited)

- **Install and boot at HEAD** (lead): `R CMD INSTALL` of both packages into an empty library, then `runDashboard()` on a real port, serves HTTP 200.
- **Cross-sid isolation holds.** `list_jobs`/`get_job` are scoped by sid; sid B cannot see or fetch sid A's job even given its job id.
- **No double finalisation.** Two pollers on one completed job run `finalize_success` exactly once.
- **Registry cleanup.** `retire_job` frees memory, and `evict_stale_jobs` handles running jobs (`completed_at = NA`) without error.
- **Cookie timing.** A fresh first visit gets the `mE_session` cookie before the websocket opens (headless Chrome, `HTTP_COOKIE` present on session 1).
- **Cancel and resume.** Banner **Cancel** finalises correctly. Same-day resume of a completed job works end to end in a browser.
- **Logging.** It records Shiny's own `session$token` and the user agent, never the `mE_session` value.
- **Uploads.** CSV upload sanitises column names (`read.csv(check.names = TRUE)`), so formula-building from column names cannot be reached through an upload.
- **Report escaping.** The HTML report escapes cells, captions, settings and the title; only the header is broken (R20.D15).
- **Backend calls resolve correctly.** All 34 `mouseExperiment::` calls use argument names that exist in `formals()`, and every UI choice value is valid for the backend's `match.arg`, except AUC with dunnett/tukey (R20.D39). A `codetools` scan finds no closure-captured helper other than `build_surv_forest_plot` and the promise operators.
- **Output slots.** Only `tg_comparisons_title` and `tg_posterior_summary` lack UI slots.
- **Earlier fixes still hold.** R3.D1 (comparison tables are the backend's `pairwise_comparisons`, `p.value == P_Value_Adjusted`); R14.6 (AUC forest labels use adjusted p); R3.D5 (survival reduction takes the last day and the maximum event flag).
- **Weight-loss KM.** The survminer p on the toxicity KM plot equals the backend's `log_rank` (5.80e-07).
- **Cox path.** The standard Cox path renders its C-index, PH banner, Schoenfeld plot and Cox table.
- **Dose-mapped labels.** They flow correctly through TG and Survival.
- **Working downloads.** Every TG download on LME4 and AUC, every Synergy and Dose-Response download, and the KM, Cox CSV and risk-table survival downloads.

## R20-D.F Open questions for the maintainer

1. Is `.rda` upload a supported feature? If not, removing it closes R20.D9 outright.
2. Should analytics carry any study-specific strings? Recommended: no (R20.D12).
3. What is the intended upload limit? 5 MB applies today; the docs say 100 MB.
4. The session snapshot: is carrying the uploaded data and results in one file a requirement, or would the JSON preset plus re-uploading the CSV be acceptable (R20.D10)?
5. Should volume units be declared by the user rather than inferred? Recommended: yes (R20.D34, package R20.83).
6. "Treat natural death as event": implement it (it needs a removal-cause column), or remove it (R20.D13)?
7. Should Dose-Response require choosing one agent's dose series (R20.D47)?
8. Should the frequentist fits move off the main process as the Bayesian ones do (R20.D6)?

## R20-D.G Status — Round 21

| ID | Issue | Severity | Status |
|---|---|---|---|
| R20.D1 | DD/MM dates read as MM/DD; wrong Day values | **Critical** | ✅ Fixed v26.09.003 (R20-D.K) |
| R20.D2 | Survival reduction merges reused IDs across arms | **Critical** | ✅ Fixed v26.09.003 (R20-D.K) |
| R20.D3 | Legacy inline Bayesian path crashes; modal stuck | Major | ✅ Fixed v26.09.002 (R20-D.J) |
| R20.D4 | `NA_Date_` crashes on non-US reference dates | Major | ✅ Fixed v26.09.002 (R20-D.J) |
| R20.D5 | Dismiss-all strands the modal and drops the result | Major | Open |
| R20.D6 | Synchronous analyses freeze all users (power sim 121 s, no upload needed) | Major | Open |
| R20.D7 | Cancel doesn't stop the worker; one worker for everyone | Major | Open |
| R20.D8 | Bayesian survival from volume-derived events always fails | Major | Open |
| R20.D9 | `.rda` upload dead, and loads into the shared global env | Major | ✅ Closed by removal v26.09.001 (R20-D.I) |
| R20.D10 | Session restore deserialises and copies unvalidated keys | Major (hardening) | ✅ Closed by removal v26.09.001 (R20-D.I) |
| R20.D11 | `test-shiny_app.R` has never executed | Major (process) | ✅ Fixed v26.09.002 — the tests run and pass; CI is package R20.53 (R20-D.J) |
| R20.D12 | Study strings sent to analytics | Major (confidentiality) | Open — the Tumour Growth and Survival reference names no longer sent (R20-D.K) |
| R20.D13 | "Natural death as event" checkbox never read | Major | ✅ Fixed v26.09.003: removed; the assumption is stated (R20-D.K). The working option comes with the removal-reason column (step 5) |
| R20.D14 | Survival/over-time synergy need an ID column named `ID` | Major | ✅ Fixed v26.09.003 (R20-D.K) |
| R20.D15 | HTML report headers merged into one cell | Major | ✅ Fixed v26.09.002 (9d5f674) |
| R20.D16 | Combined E-T table and download always fail | Major | ✅ Closed by removal v26.09.001 (R20-D.I) |
| R20.D17 | Power tab ignores multiplicity and attrition | Major | ✅ Fixed v26.09.002 (R20-D.J) |
| R20.D18 | Bayesian results keep their fit alive via `plot_env` | Major (memory) | Open |
| R20.D19 | Survival forest-plot download calls an out-of-scope function (HTTP 500) | Major | ✅ Fixed v26.09.002 (R20-D.J) |
| R20.D34 | Toxicity mass correction 1000× on the bundled demo; warning swallowed | **Critical** | ✅ Fixed v26.09.003 (R20-D.K) |
| R20.D35 | Synergy/TWM intervals shown as uncertainty for estimates they don't describe | **Critical** | Open (package R20.2) |
| R20.D36 | Synergy fails after TG when a dose column is mapped | Major | ✅ Fixed v26.09.003 (R20-D.K) |
| R20.D37 | Synergy TGI plot re-derived from survivors; contradicts the table | Major | Open |
| R20.D38 | Synergy Info tab documents response additivity, not Bliss | Major | Open |
| R20.D39 | AUC fails on first click with the default correction | Major | ✅ Fixed v26.09.002 (R20-D.J) |
| R20.D40 | GAMM forest plot crashes; download fails; ZIP omits it | Major | ✅ Closed by removal v26.09.001 (R20-D.I) |
| R20.D41 | Dashboard overwrites the backend's GAMM residual plot with a broken one | Major | ✅ Closed by removal v26.09.001 (R20-D.I) |
| R20.D42 | Q-Q plots against y = x on unstandardised values | Major | Open |
| R20.D43 | Forest plots: unadjusted CI with adjusted p | Major | Open |
| R20.D44 | Comparisons tab: mislabelled correction, discarded title, live family | Major | Open |
| R20.D45 | Data Summary log-scale means labelled mm³ | Major | Open |
| R20.D46 | HTML report records widget state, not run settings | Major | Open |
| R20.D47 | Dose-Response pools drugs by dose value | Major | Open |
| R20.D48 | DR "All dates" = own last observation; TGI plot contradicts Toxicity | Major | Open |
| R20.D49 | JT test and EC50 CI returned but never shown | Major | Open |
| R20.D50 | Toxicity plot downloads never produce a file | Major | Open |
| R20.D51 | Power "seed SD" pools all days (N inflated 5–40×) | Major | Open |
| R20.D52 | Download-All ZIP can include another session's files | Major (confidentiality) | ✅ Fixed v26.09.002 (9d5f674) |
| R20.D20 | Synergy bar plot all grey | Minor | Open |
| R20.D21–R20.D33, R20.D54–R20.D61 | Minor | Minor | Open. Fixed: R20.D21, R20.D33 (R20-D.J); R20.D31 (R20-D.K). In part: R20.D29, R20.D54, R20.D57 (R20-D.K) |
| R20.D53 | Plot-label keystrokes rebuild the growth plot | Efficiency | Open |

**Suggested order of work:**
1. **The one-liners.**
   - `NA_Date_` → `as.Date(NA)` (R20.D4)
   - the promise import (R20.D3)
   - `vapply(colnames(df), esc, "")` (R20.D15)
   - `id_column = id_col_surv` (R20.D14)
   - the real column names in the E-T report (R20.D16)
   - `tempfile()` for the ZIP (R20.D52)
   - move `build_surv_forest_plot` to top level (R20.D19)
   - drop the `message =` argument in `test-shiny_app.R` (R20.D11)
2. **Correctness of what users read.**
   - date-order detection (R20.D1)
   - the volume-unit selector plus a warning-surfacing `safe_analysis()` (R20.D34)
   - the survival reduction key (R20.D2)
   - suppress or relabel the model-estimand intervals (R20.D35)
3. **Security hardening.** Remove `.rda` upload (R20.D9); move session restore to a data-only format (R20.D10); strip study strings from analytics (R20.D12).
4. **Delete the dashboard's statistical re-derivations** (R20.D37, R20.D41, R20.D42, R20.D48) in favour of the backend's outputs. That is the recurring source of contradictions here.
5. **Move the long synchronous fits into the job path** (R20.D6, R20.D7).

## R20-D.H Maintainer decisions (2026-09-29)

- **Volume units are declared by the user** (R20.D34). Data Upload gains a required mm³/cm³ selector, and the value is passed as `volume_units` to every backend call that converts volume to mass. The dashboard's own `/1000` re-derivation is removed.
- **Survival: the user declares the unit of randomisation** (mouse or cage) when a cage column is mapped. Mouse gives Cox/Firth; cage gives a cage-level permutation log-rank with its resolution floor displayed (package R20-K). With 2–3 cages per arm a cage-randomised comparison cannot reach p < 0.05, and the tab must say so.
- **Synergy and TGI are reported only on days with enough animals on study in every arm** (package R20-K). The "All dates" view becomes "all evaluable days". **Rule confirmed 2026-09-30:** an arm is evaluable on a day when at least 50 % of its enrolled animals, and at least 3 animals, are still on study. **It must be relayed to users.** The Synergy, Toxicity (Therapeutic Window) and Dose-Response Info tabs state it, and every result that uses an evaluable day names the day used and lists the days and arms the rule excluded. This lands with the rule itself (step 4); stating it earlier would describe behaviour the code does not yet have.
- **`.rda`/`.RData` upload is removed** (CSV upload stays), which closes R20.D9.
- **Tumour growth:** the default random effects become correlated random slopes with fallback. **GAMM is removed** (with R20.D40 and R20.D41). **AUC is kept but must be redefined** (package R20.4 and R20.D39).
- **Bayesian scope:** only **Bayesian tumour growth (LMM)** remains, to be made robust (compiled-model caching, the R20.D3/D5/D7/D18 fixes, resources). **Bayesian synergy and Bayesian dose-response are removed.**
- **Toxicity:** Body-weight AUC, Total Benefit Area, Efficacy–Toxicity bivariate and weight-corrected TGI are removed, along with the Combined Report (R20.D16). The tab keeps the body-weight model, time to weight-loss threshold, and the therapeutic window.
- **Bayesian scope (revised 2026-09-29):** keep **tumour growth and survival**; remove body weight and therapeutic window as well. The F.7 resume layer keeps two job kinds (`bayes_tg`, `bayes_surv`); `bw_resume.R`, `syn_resume.R`, `dr_resume.R` and `twm_resume.R` go.
- **Therapeutic window:** a two-axis view (TGI vs worst weight loss, each with intervals, plus a tolerability flag) replaces the ratio (package R20-K).
- **Removal reasons:** an optional removal-reason column on Data Upload. The "natural death as event" checkbox appears only when it is mapped; otherwise the tab states its assumption.
- **Bayesian survival families:** Weibull and log-normal only; exponential and gamma are removed from the selector.
- **Session Save/Restore: removed** (decided 2026-09-29). The UI claims it saves "settings, results, plots" and that "filter state survives session save/restore". The snapshot (`mod_log_session.R:81-96`) actually holds only the uploaded data, column mappings and plot-group settings: no results, no plots, no subset filters. Re-uploading the CSV and loading the JSON settings preset gives the same result without deserialising an uploaded R object. Keep the preset and extend it with the study-design settings and subset filters. This also closes R20.D10.


---

## R20-D.I Implementation log — v26.09.001, step 1: removals (2026-09-29)

This step applies the scope decisions in R20-D.H and needs backend 0.23.0 (package R20-L).

**Removed:**
- **Bayesian:** every fit except tumour growth and survival. Bayesian synergy, dose-response, body weight and therapeutic window go, with their settings, tabs, diagnostics, downloads and resume handlers.
- **Tumour growth:** GAMM and Bayesian GAMM.
- **Survival:** the exponential and gamma families.
- **Toxicity:** body-weight AUC (with its nadir and comparison tables, the forest plot and the comparison-type selector), the safety-efficacy scatter, total benefit score, weight-corrected TGI and the Combined Report. The efficacy-toxicity mode is relabelled **Therapeutic Window**, since the TWM is all it now holds.
- **Upload and session:** `.rda`/`.RData` upload, session save/restore, and the dead `mod_efficacy_toxicity.R` stubs.

**Closed by removal:** R20.D9, R20.D10, R20.D16, R20.D40 and R20.D41 (Major), and R20.D56 (Minor: both halves concerned GAMM and Bayesian synergy).

**Partly closed; the rest stays open:**

| ID | Gone with the removals | Still open |
|---|---|---|
| R20.D3 | four of the six legacy inline fallbacks | tumour growth and survival (`tg_resume.R`, `surv_resume.R`) |
| R20.D30 | the Type III text shown for the GAM smooth table | the other three explainer errors |
| R20.D50 | five of the eight plot types | the trajectory, KM and TWM downloads still read an output and fail |
| R20.D58 | the BW AUC adjustment mismatch | the silently dropped "Tumor Volume" covariate |
| R20.D61 | the GAMM "EMMs at five quantile days" promise | the other items |

**Verification (executed):**
- `devtools::test()` against backend 0.23.0 installed in a scratch library: 184 tests, 564 expectations, 0 failures. The drop from 574 is the four removed resume kinds (8 expectations) and the two GAMM labels in `test-docs-current.R` (2). `test-shiny_app.R` still errors at file level whenever shinytest2 is installed (R20.D11); the fix is in step 2.
- `testServer()` smoke on the master demo, with results, plots and downloads read back: 37 of 37 checks pass. The runs were body weight (7 outputs, 3 plots, 4 diagnostics, 4 downloads including the ZIP), Therapeutic Window (summary, table, ranking note, attrition, plot, plot groups, 2 downloads), tumour growth LME4 and AUC, survival, synergy (single day and over time, with CSV download) and dose-response.
- Install plus `runDashboard()` gives HTTP 200. The served page has 604 element IDs and no duplicates. None of the removed input IDs or labels appears; the only `exponential` values left are the tumour-growth prior-distribution choices.
- Headless Chrome ran the same script against this build and 26.07.015:
  - Toxicity now offers `bw,et`.
  - Body weight inserts Model Results and Weight Loss KM. Therapeutic Window inserts only its own tab.
  - The TWM table and plot render.
  - The TG choices are `lme4,auc,bayes`, the survival families `weibull,lognormal`, and no Bayesian selectors remain on Synergy or Dose-Response.
  - Both builds log the same 20 console errors, all shinyFeedback on demo load (R20.D33), so none are new.
- `R CMD check --no-tests --no-manual --as-cran`: 6 WARNINGs, 5 NOTEs, as at baseline. The detailed output differs only by lines that disappeared with the removed files.

**Harness notes (no new findings):**
- With a dose column mapped, tumour growth relabels arms "Treatment Dose" and offers those labels as the reference choices (`mod_tumor_growth.R:20-56`). A `testServer()` run must therefore pass "Vehicle 0", not "Vehicle". The same relabelling, written back into `rv$tx_col`, is what breaks synergy in R20.D36.
- AUC with the default Dunnett correction fails exactly as R20.D39 describes.
- The Toxicity diagnostics stay empty after switching back from Therapeutic Window to Body Weight until the analysis is re-run. That is the same in 26.07.015: results belong to the mode that produced them.

**Next:** step 2, the one-line and near-one-line bugs.


---

## R20-D.J Implementation log — v26.09.002, step 2: one-line and near-one-line bugs (2026-09-30)

This step needs backend 0.24.0 (package R20-M).

**Fixed:**

| ID | Fix | Verified by |
|---|---|---|
| R20.D3 | Sessions without the identity cookie get a registry key, `make_local_sid()`: "local-" + 64 hex. `get_sid()` accepts only 64 hex characters from a Cookie header, so no browser can present the key. Both `.legacy_inline_*` paths are deleted, and the launchers fall back to a local key when called without one. | testServer: a cookie-less launch lands in the registry under a local key. Browser: see the note below this table. |
| R20.D4 | An empty or unparseable reference date stops with "Enter it as YYYY-MM-DD", instead of `NA_Date_`. **Deviation from the suggested fix:** the review proposed `as.Date(NA)` and keeping the `Sys.Date()` fallback, but that would count every Day from today without a word. Tumour Growth and Survival now report the error and reset their run flags, and the Data tab shows it under the reference-date box. | testServer (TG run with "24/03/2025"); browser (message appears, and clears once corrected). |
| R20.D11 | The skip call no longer passes `message =`. The test apps are now app directories that load the package themselves: AppDriver cannot carry a `load_all()` package into its process. Also fixed: the pre-module signatures, test hooks missing from the page, `get_html(timeout =)`, and a workflow hook that passed on a failed fit. | 8 browser tests run and pass. CI remains package R20.53. |
| R20.D17 | The analytic Power tab passes `p_adjust_method` (Bonferroni by default), `n_comparisons` (follows k − 1) and `dropout_rate`. The summary states the family, the attrition and the backend's method note. The submode labels (`find_n`/`find_power`) were never matched and always printed "analytic"; now fixed. | testServer: k = 4, d = 1 gives 23 per group (was 7), Enroll_N 29 at 20 % attrition. Browser: the same, plus the comparison count following the group count. |
| R20.D19 | `build_surv_forest_plot()` moved to `helpers_plot.R`. | testServer and browser: the download is a real PNG (it returned HTTP 500 before). R CMD check no longer lists the symbol as undefined. |
| R20.D21 | `icon("spinner", class = "fa-spin")`. | The test fails on 26.09.001 (the notice is a message, not a warning). |
| R20.D33 | **Different root cause from the one recorded.** The empty feedback text was not the problem: feedback was sent before the mapping selectors were bound, and shinyFeedback cannot look up a missing binding. The validator now skips inputs that read NULL. | Browser: 0 console errors on demo load (20 before, on both builds). |
| R20.D39 | `tg_p_adjust_choices()` drops Dunnett/Tukey for AUC. Holm is selected on the way into AUC, and Dunnett returns on the way out; the server maps a stale Dunnett or Tukey to Holm. | testServer (a stale Dunnett runs with `p_adjust_method_used = "holm"`); browser (choices, run, restore). |
| R20.8 (dashboard side) | The synergy summary says "Not evaluable" instead of "Synergistic: NA", reports evaluable timepoints, and has no peak line when there is no peak. The Info tab documents `Evaluable`. | testServer on real non-evaluable results, single day and over time. |

**R20.D3 end to end.** In headless Chrome with cookies disabled through the DevTools protocol (`document.cookie` empty), a Bayesian tumour-growth fit (1 chain, 200 + 200 draws) ran through the registry. The modal closed, and the tab reported "Bayesian analysis completed successfully!" with a Bayesian summary. That run took 28 s; an earlier run in the same worker had compiled the model. Before the fix this case left the modal open with `Error in %...!%` in the log.

Chrome treats `http://127.0.0.1` as a secure context and stores the Secure cookie there. The README's claim that local sessions have no cookie is therefore only sometimes true, but either way they now work.

**Found and fixed during this step:**
- **A start-up crash introduced and caught within the step.** The first R20.D3 edit read `rv$session_id` in the top-level server body. Every real session closed at start-up, while the page returned HTTP 200 and `testServer()` passed; it runs the body inside a reactive context. The headless browser check caught it. `runDashboard()` now builds its app with `dashboard_app()`, and a new browser test starts the full app. That test fails on the broken line (session not connected, no disclaimer, demo load times out) and passes on the fix.
- **The forest-plot axis label** rendered "log₁₀" as "log" plus two boxes on screen and in the download (no subscript glyphs in the PNG font). It now reads "log10".
- **A test-suite hang.** `.reset_ttl_sweeper()` reset a flag but left its 9,999 s callback on the global `later` loop. `MockShinySession$getOutput()` waits for that loop to empty (`shiny:::wait_for_it`), so every later test file that read an output hung at 0 % CPU. The existing files sorted before the sweeper test; `test-round21-fixes.R` did not. The reset now cancels the callback, and a test asserts the loop is empty after reset.
- The Bayesian survival modal said "Weibull" for log-normal fits.

**Tests.**
- `test-round21-fixes.R` holds 9 tests and 34 expectations. Run against the 26.09.001 source with backend 0.24.0, 8 of 9 fail or error. The one that passed there was the icon test; it now uses `expect_no_message()` and fails too.
- `test-shiny_app.R` holds 8 browser tests, including the new full-dashboard test; `test-resume_banner.R` gains the sweeper-reset test.

**Verification (executed):**
- Full suite: 21 files, 201 tests, 612 expectations. No failures, no skips; the one warning is pre-existing.
- R CMD check: 6 WARNINGs, 5 NOTEs, as at baseline. Against 26.09.001 the only differences are notes that disappeared: `%...>%`, `%...!%`, `NA_Date_` and `build_surv_forest_plot` are no longer undefined.
- Install plus boot gives HTTP 200, and the headless-Chrome script passes all checks with 0 console errors.

## R20-D.K Implementation log — v26.09.003, step 3: study design, one animal key, declared units, dates (2026-09-30)

This step needs backend 0.25.2. That is package R20-N plus two follow-ups found while wiring the dashboard, both logged there:
- v0.25.1: `make_mouse_key()` accepts a frame with no rows.
- v0.25.2: two spurious warnings raised on every tumour-growth fit are gone.

**The Study Design panel (maintainer decision R20-D.H).** Data Upload gains a panel below the column mapping, and every tab uses it:
- **Volume units, mm³ or cm³.** The demos, and uploads whose volume is calculated, start at mm³; an upload with a volume column starts undeclared.
- **Control group.** It replaces the five per-tab reference selectors.
- **Unit of randomisation.** Shown when a cage column is mapped; it replaces the Tumour Growth permutation-unit selector.

`study_reference()` maps the control into each tab's labels ("Vehicle" → "Vehicle 0" with a dose column). Each tab states the reference it used, and says so when the control is excluded by the subset filter. The panel renders once per dataset, so remapping a column does not reset the declared units. The preset (format 2) saves the design and the subset filter, and refuses to load before a dataset exists; before, its updates were lost.

**Fixed:**

| ID | Fix | Verified by |
|---|---|---|
| R20.D1 | `detect_date_format()` reads the whole column. A first number above 12 means day-first, a second above 12 means month-first; with neither, the declared order (new selector) is used, or the user is asked. Contradictory data, or data contradicting the declared order, are errors. Year-less dates take the year nearest the reference date. Month names are matched without the locale. `study_days()` hands the backend YYYY-MM-DD and reports unreadable values, mixed layouts, a New Year crossing and days before Day 0. Data Upload states the layout and the reason. All three date callers (preview, `ensure_preprocessed()`, Tumour Growth) use it. | The review's case: 0, 3, 7, 10, 14, 17 (was −57, 33, 155, 247, NA, NA). Ambiguous, declared, contradicting, New-Year and negative-day cases; the declared order through `ensure_preprocessed()`; the Data Upload note. Browser: the Combo demo reads MM/DD with its reason, and declaring DD/MM is reported as "would need month 24". |
| R20.D2 | `survival_animal_rows()` keys animals by arm + ID (+ cage). The event time is the first flagged day; when every row is flagged, the flag is an animal-level status and the last day is used. `volume_threshold_events()` uses the same key. The demos flag only the event row, so their results do not change. | testServer: 10 + 10 animals with tags restarting per arm and no cage now give 20 records and 14 events (10 before), and the animal flagged on day 20 but measured to day 28 has time 20. Unit test of both flag encodings. |
| R20.D13 | The unread checkbox is removed. The sidebar states the assumption of the chosen event source: in column mode events are the rows coded 1; in volume mode other removals are censored, which favours the treatment if they are treatment-related. The Info tab says the same. | testServer and UI-source test. |
| R20.D14 | Survival passes `id_column`, `cage_column`, `randomisation_unit` and `verbose = FALSE`. "All dates" synergy passes `id_column`, `cage_column`, `endpoint_method`, `n_boot` and `boot_seed`; the bootstrap costs 0.7 s for 11 days at 200 resamples on the Master demo. | testServer: survival and over-time synergy with the ID column named "Animal" (both failed before). |
| R20.D34 | `volume_units` goes to `analyze_body_weight()`, `weight_loss_threshold()` and `therapeutic_window_metric()`. With units undeclared, the correction is blocked with a message. The sidebar formula follows the declared unit. The trajectory plot draws the model's `weight_data$Net_Weight` per animal key and labels net or gross weight. `safe_analysis()` collects warnings with `withCallingHandlers()`, logs them and attaches them; `completion_status()` lists up to five. | testServer on the Weight demo: 0 of 20 threshold events (20 of 20 before), all net weights above 10 g, and the plot data equal to the model's net weights. With undeclared units the run is blocked. Browser: declaring cm³ is refused with the backend's "exceeds 50 % of body weight". |
| R20.D36 | Synergy's selectors and analysis both use `arm_labels()` ("Treatment Dose" when a dose column is mapped), independent of whether Tumour Growth has run. | testServer: TG runs first, then Synergy runs with its own rendered choices ("groups do not exist" before). Browser: the same on the Master demo. |
| R20.D31 | `cage_column` goes to `therapeutic_window_metric()`, and the trajectory re-derivation is gone (see R20.D34). | Mocked backend call receives `cage_column = "Cage"` and `volume_units = "mm3"`. |
| Package R20.6, dashboard half | Survival passes the declared unit of randomisation. Model Results shows the cage-permutation table (cages per arm, assignments, p, smallest attainable p) and the backend's caveat when animals were randomised but housed one arm per cage. The floor warning reaches the status box. | testServer and browser on the Master demo, where 2 cages per arm give a floor of 1/3. |

**Fixed in part:**
- **R20.D12:** the Tumour Growth and Survival analytics events no longer send the reference name.
- **R20.D29:** Toxicity and Dose-Response use `ensure_preprocessed()`, so the current subset filter applies. The KM p-value label and `coord_cartesian` remain.
- **R20.D46:** the Tumour Growth report records the reference of the displayed run. The other settings still come from the widgets.
- **R20.D48:** the Dose-Response TGI plot uses the control stored at run time. The "All dates" label and the TGI estimand remain.
- **R20.D54:** the P_Value help now says the value is adjusted, and P_Method is documented. The other help errors remain.
- **R20.D57:** the drill-down matches treatment + ID + cage (`drilldown_rows()`). The geometric-mean bar, the forest row and the log1p text remain.

**Found and fixed during this step:**
- **The subset filter excluded an ear tag from every arm.** It now lists and excludes animals by arm, ID and cage, for example "3 (Vehicle 0, cage C1)".
- **Stale preprocessing cache.** Changing the reference date, date order, volume formula, height, cage or dose column did not clear the cache. Tabs reading it then disagreed with Tumour Growth, which recomputes on every run.
- **The completion message rendered in fragments** ("Analysis completed ⏎ with ⏎ 1"). It is now built as one string.
- **A dead Reset handler.** It had no button, and it replaced the mapping, summary and preview renderers with permanent NULLs; a Reset button added later would have blanked the next dataset. It was removed rather than fixed.

**Status rows brought up to date.** R20.D15 (report headers) and R20.D52 (private ZIP directory) were fixed in commit 9d5f674, after the R20-D.J log was written; their rows now say so.

**Surfaced but not fixed.** Showing warnings exposed two backend defects, left for the package:
- **Weight-loss Cox fit with few or no events.** `weight_loss_threshold()` still fits Cox/Firth when no animal crosses the threshold. `coxphf` then reports "NAs introduced by coercion" and "Numerical problem in parameter estimation" on the Weight and Master demos.
- **"NaNs produced" from Dose-Response on the Master demo.** To be traced in step 5, which reworks Dose-Response.

**Tests.**
- `test-round21-study-design.R` holds 24 tests and 83 expectations. Run against the v26.09.002 source with backend 0.25.1, all 24 fail or error; the module tests set the old inputs, so they fail on the defect, not on a missing input.
- `test-shiny_app.R` gains a browser test of the panel and the date note (9 browser tests).
- Updated to the new contracts:
  - `test-helpers_data.R`: the date layouts are asserted as day values, and filter exclusions are animal keys; it gains a reused-tag exclusion test.
  - `test-round21-fixes.R`: the Tumour Growth tests set the study's control instead of `tg_reference`, and the R20.D4 tests use a date that does not exist, since "24/03/2025" is now read.

**Verification (executed):**
- **Full suite:** 22 files, 228 tests, 704 expectations. No failures or skips; the one warning is pre-existing.
- **Package suite at 0.25.2:** 286 tests, 0 failed, the 1 legitimate skip. It records 12 fewer warnings than at step 2 and the same assertions.
- **R CMD check:** 6 WARNINGs and 5 NOTEs, in the same categories as the step 2 baseline. No new undefined symbols, and `get_groups` (the removed survival selector) is gone.
- **Install and boot:** HTTP 200.
- **Headless-Chrome script:** 18 checks across three sessions, all passing, with 0 console errors:
  - Master: every tab runs from the study design; survival by animal and by cage; synergy after Tumour Growth; the preset JSON.
  - Weight demo: mm³ correct, cm³ refused.
  - Combo: date reading, the contradiction message, and Tumour Growth against Control.

