# Architecture

How the dashboard is organized and how a request flows through it.

This document is the entry point for new contributors. After reading it, you should be able to (a) find the file that owns any given UI surface, (b) trace what happens when a user clicks "Run Analysis", and (c) understand why the F.7 session-resume layer exists and what files participate in it.

For statistical method choices, see [`docs/METHODS.md`](METHODS.md) in the backend repo.
For the data-format spec, see [`docs/DATA_FORMAT.md`](DATA_FORMAT.md).
For the F.7 isolation guarantees and threat model, see [`docs/SECURITY.md`](SECURITY.md).

---

## High-level layout

The dashboard is a Shiny `navbarPage` with one tab per analysis module. Every module is a self-contained `shiny::moduleServer()` pair (UI + server). A small set of cross-cutting helpers (`helpers.R`, `helpers_bayes.R`, `helpers_report.R`) provides services every module uses. A separate F.7 layer (`job_registry.R` + per-module `*_resume.R` files) keeps long-running Bayesian fits resumable across tab disconnects.

```
┌─────────────────────────────────────────────────────────────────────┐
│  navbarPage (R/app.R)                                                │
│  ┌──────────────┬──────────┬──────────┬──────────┬──────────┐        │
│  │ Data Upload  │ TG       │ Survival │ Toxicity │ DR / Syn │ Power  │
│  └──────────────┴──────────┴──────────┴──────────┴──────────┘        │
│                                                                      │
│  Each tab = one module:                                              │
│    UI:  fooUI(id)            (sidebar + results tabsets)             │
│    Srv: foo_server(id, rv)   (observers + fit dispatch + renderers)  │
│                                                                      │
│  Resume banner sits above the tabs (renders nothing when sid has     │
│  no jobs).                                                           │
└─────────────────────────────────────────────────────────────────────┘
              ▲
              │
       rv (reactiveValues) flows between modules:
         uploaded_data, processed_data, detected_cols, settings,
         <X>_results, <X>_running, <X>_job_id, session_id, bayes_models
```

`rv` is created once in `server()` (in `R/app.R`) and passed to every `<module>_server(id, rv)` call. Modules write their results back to `rv$<X>_results`.

---

## Directory layout

```
R/
├── app.R                          # entry point: runDashboard()
├── helpers.R                      # core utilities
├── helpers_bayes.R                # Bayesian diagnostics panel + LOO helpers
├── helpers_report.R               # standalone HTML report builder
├── logging.R                      # structured JSON logging
├── ui_modules.R                   # shared UI components
├── zzz.R                          # package onLoad hooks
│
├── mod_data_upload.R              # Data tab
├── mod_log_session.R              # Log, filters and settings-preset panel
├── mod_power.R                    # Power Analysis tab
├── mod_tumor_growth.R             # TG lifecycle (observers, fit dispatch)
├── mod_tumor_growth_plots.R       #   TG plot reactives + renderPlots
├── mod_tumor_growth_downloads.R   #   TG download handlers
├── mod_survival.R                 # Survival lifecycle
├── mod_survival_downloads.R       #   Survival downloads
├── mod_toxicity.R                 # Toxicity lifecycle (BW + ET)
├── mod_toxicity_plots.R           #   Toxicity plot renderers
├── mod_toxicity_downloads.R       #   Toxicity downloads
├── mod_synergy.R                  # Synergy lifecycle
├── mod_synergy_downloads.R        #   Synergy downloads
├── mod_dose_response.R            # DR lifecycle
│
├── job_registry.R                 # F.7 process-level job registry
├── session_identity.R             # F.7 server-side cookie reader
├── resume_banner.R                # F.7 banner UI/server + dispatch
├── bayes_job_runner.R             # F.7 shared poller + resume registration
├── tg_resume.R                    # F.7 per-kind handlers for the two
├── surv_resume.R                  #     Bayesian fits (TG, survival)
└── ttl_sweeper.R                  # F.7 periodic registry housekeeping

inst/www/
├── css/                           # custom + dark-mode styles
└── js/
    ├── session-id.js              # F.7 cookie setter (first connect)
    ├── navbar.js
    ├── sidebar.js
    ├── darkmode.js
    ├── overlay.js                 # Bayesian "fit running" overlay
    └── umami-events.js            # analytics tracking

tests/testthat/
├── helper-app.R
├── test-demo_data.R
├── test-helpers.R
├── test-modal-and-cachebust.R
├── test-module_servers.R
├── test-new-helpers.R
├── test-session_isolation.R       # F.7 security audit (39 assertions)
├── test-resume_banner.R           # F.7 banner module (10 assertions)
├── test-tg_resume.R               # F.7 per-kind helpers (11 assertions)
├── test-resume_e2e.R              # F.7 disconnect/reconnect cycle (26 assertions)
└── test-shiny_app.R               # shinytest2 (Data tab smoke test)
```

### File-naming convention (CODE_REVIEW.md C.1/C.3)

Module files follow a three-suffix convention so a 2,000-LOC `moduleServer` body doesn't all live in one file:

| Suffix | What lives there |
|---|---|
| `mod_<X>.R` | Lifecycle: observers, fit dispatch, status messages, tab insertion, table renderers |
| `mod_<X>_plots.R` | Plot reactives + their `renderPlot` bindings + group-order/palette UI + diagnostic plots |
| `mod_<X>_downloads.R` | All `downloadHandler()` bindings (CSV / PNG / ZIP) |

Each helper is `add_<X>_<concern>(input, output, session, rv, ...)` called once at the bottom of the main module's `moduleServer` body. Plot reactives needed by both `_plots.R` and `_downloads.R` are passed explicitly as arguments rather than captured by closure, so the cross-file dependency is visible at the call site.

Only `mod_tumor_growth.R` and `mod_toxicity.R` have all three (they're the largest); `mod_survival.R` and `mod_synergy.R` have downloads-only; `mod_dose_response.R` has none (it was small enough to leave whole).

---

## What happens when the user clicks "Run Analysis"

Frequentist and Bayesian paths differ. The frequentist path is synchronous and contained inside the module. The Bayesian path goes through the F.7 registry.

### Frequentist (e.g. LME4 tumor growth)

```
User clicks Run
       │
       ▼
observeEvent(input$tg_run_analysis, { ... }) in mod_tumor_growth.R
       │
       ├─ rv$tg_running <- TRUE
       ├─ show_analysis_modal("Fitting tumor growth model…", ...)
       ├─ on.exit(shiny::removeModal())          # auto-dismiss
       ├─ result <- safe_analysis(module = "tumor_growth", {
       │      mouseExperiment::tumor_growth_statistics(df = data, ...)
       │   })
       ├─ rv$tg_results <- result
       ├─ updateTabsetPanel(session, "tg_results_tabs", selected = "Summary")
       ├─ rv$tg_running <- FALSE
       └─ output$tg_status_msg <- renderUI(status_alert("…", "success"))
```

`safe_analysis()` (in `helpers.R`) wraps any expression in `tryCatch` + structured logging so errors don't blow up the session. Total runtime: 5-30 seconds for typical preclinical datasets.

### Bayesian (tumor growth or survival)

```
User clicks Run
       │
       ▼
observeEvent(input$<X>_run_analysis, { ... }) in mod_<X>.R
       │
       ├─ Build bayes_params <- list(data, col_names, prior_strength, …)
       │
       └─ launch_<X>_bayes_job(rv$session_id, bayes_params,
                               input, output, session, rv)
                                       │
                                       ▼
                          R/<X>_resume.R
                          ┌─────────────────────────────────────────┐
                          │ 1. show_analysis_modal("Bayesian …")    │
                          │                                         │
                          │ 2. promise <- promises::future_promise( │
                          │      mouseExperiment::bayesian_<X>(...) │
                          │    )                                    │
                          │                                         │
                          │ 3. job_id <- launch_job(                │
                          │      sid, kind = "bayes_<X>", promise,  │
                          │      context = list(bayes_params=…))    │
                          │                                         │
                          │ 4. rv$<X>_job_id <- job_id              │
                          │                                         │
                          │ 5. attach_bayes_poller(input, output,   │
                          │      session, rv,                       │
                          │      job_id_field    = "<X>_job_id",    │
                          │      finalize_success = .finalize_<X>…) │
                          └─────────────────────────────────────────┘
                                       │
                                       ▼
                       Promise is now held by .job_registry
                       (process-level env). The future runs in a
                       multisession worker independently of the
                       Shiny session.

                       The polling observe() in the module:
                       - reads rv$<X>_job_id
                       - calls get_job(sid, job_id) every 2 s
                       - on terminal status, dispatches to
                         .finalize_<X>_success / _failure / _cancellation
                       - retires the registry entry
                       - clears rv$<X>_job_id to stop polling
```

The finalizers are the same code path whether the fit completed inline or was resumed from disconnect — they read the result from the registry job entry, populate `rv$<X>_results`, switch to the Summary tab, and remove the modal.

For Bayesian fits, runtime is 3-12 minutes (Stan compile + MCMC).

---

## The F.7 session-resume layer

When a Bayesian fit is running and the user's tab disconnects (NPM proxy timeout, browser sleep, network blip), three things happen:

1. The Shiny session ends. `rv` is garbage-collected.
2. The future continues running in its worker (workers don't know about Shiny).
3. The `.job_registry` env (package-level, process-scoped) still holds the promise.

On reconnect:

1. The user's browser sends the `mE_session` cookie they got on first connect (256-bit CSPRNG token, `SameSite=Strict; Secure; max-age=86400`).
2. The server reads the cookie from `session$request$HTTP_COOKIE` (server-side, not from any `input$`).
3. `rv$session_id` is set to the cookie value.
4. The resume banner (mounted in `R/app.R`) renders an alert when `list_jobs(rv$session_id)` returns any jobs.
5. The user clicks Resume on the relevant row.
6. The banner calls `dispatch_resume(job)`, which looks up the handler registered for `job$kind`.
7. That handler (registered at the bottom of each `<X>_server()` via `register_<X>_resume()`) sets `rv$<X>_job_id` and reattaches the poller — exactly as if the fit had been started in this session.
8. The poller's next 2-second tick picks up the completed result and populates `rv$<X>_results`.

```
┌──────────────────────────────────────────────────────────────────┐
│  Shiny container (single R process)                              │
│                                                                  │
│  ┌────────────────────┐   ┌─────────────────────────────────┐    │
│  │ Session A (alive)  │   │  .job_registry (package env)    │    │
│  │   rv$session_id=X  │   │  sessions: list(                │    │
│  │   rv$tg_job_id ────┼──▶│    "X" = list(                  │    │
│  └────────────────────┘   │      "job_…" = list(            │    │
│                           │        kind, status, result,    │    │
│  ┌────────────────────┐   │        promise, context,        │    │
│  │ Session B (dropped)│   │        started_at)),            │    │
│  │   was: rv$sid=Y    │   │    "Y" = list(...))             │    │
│  └────────────────────┘   │                                 │    │
│           ▲               │  future workers (multisession)  │    │
│           │ reconnects    │  ↑ run to completion regardless │    │
│           │ with cookie Y │    of which session started them│    │
│           ▼               └─────────────────────────────────┘    │
│  ┌────────────────────┐                                          │
│  │ Session B' (new)   │ banner: "You have 1 TG fit running"      │
│  │   reads cookie Y   │ → click Resume → register_tg_resume      │
│  │   → finds job(s)   │   handler reattaches the poller          │
│  └────────────────────┘                                          │
└──────────────────────────────────────────────────────────────────┘
```

The registry is multi-user-isolated via four independent layers; see [`SECURITY.md`](SECURITY.md) for the full audit.

### Where each F.7 piece lives

| File | Responsibility |
|---|---|
| `R/job_registry.R` | The `.job_registry` env. `launch_job()`, `get_job()`, `list_jobs()`, `cancel_job()`, `retire_job()`, `evict_stale_jobs()`. |
| `R/session_identity.R` | `get_sid(session)` — reads the cookie from `session$request$HTTP_COOKIE`. Never reads any Shiny input. |
| `inst/www/js/session-id.js` | Client-side cookie setter on first connect (256-bit CSPRNG). |
| `R/resume_banner.R` | Banner UI + server. Reads `list_jobs(sid)`, renders rows, dispatches Resume/Cancel clicks via `register_resume_handler` / `dispatch_resume`. |
| `R/bayes_job_runner.R` | Generic `attach_bayes_poller()` + `register_bayes_resume()`. Per-module files compose these. |
| `R/tg_resume.R` etc. | Per-kind: `launch_<X>_bayes_job()`, `register_<X>_resume()`, the three finalizers, the `.legacy_inline_<X>_bayes()` fallback for when no sid is available. |
| `R/ttl_sweeper.R` | `start_ttl_sweeper()` — `later::later`-driven periodic eviction. Called once per R process from `runDashboard()`. |

---

## How async actually works

`runDashboard()` sets up the future plan once:

```r
future::plan(future::multisession, workers = max(2L, future::availableCores() - 1L))
```

Each `promises::future_promise({...})` evaluates its body in one of the worker processes. Workers are persistent across jobs (they don't re-launch per fit). The body of the promise captures whatever it needs via lexical scope — in our case, just `bayes_params` (a list of plain R objects).

When a worker completes:
- The result is serialized back to the parent process.
- The `promises::then(...)` callbacks fire on the main event loop.
- For F.7-routed fits, that callback is the registry-write hook inside `launch_job()`, which sets `status = "completed"` and stashes the result.

The polling `observe()` in the parent module periodically checks the registry and consumes the result. It runs every 2 seconds via `invalidateLater(2000)`. This is intentional: polling decouples the parent module from `then()` callbacks tied to a specific session, so resume from a different session works the same way.

### Worker count tuning

Default `availableCores() - 1` leaves one core for the main Shiny event loop. On the VPS this is typically 3-7 workers. If a user kicks off N+1 fits (where N is the worker count), the (N+1)th queues in the future scheduler — its registry entry stays `status = "running"` from the user's perspective. Workers are not preempted; cancellation marks the registry entry but lets the future complete.

---

## Shared cross-cutting helpers

### `R/helpers.R`

| Function | Purpose |
|---|---|
| `safe_analysis(module, expr)` | tryCatch wrapper around any analysis expression; logs failures, returns `list(error = TRUE, message = ...)` on error |
| `status_alert(msg, type)` | Bootstrap-styled `<div class="alert ...">` for status banners |
| `show_analysis_modal(title, detail)` | Generic "fit running" modal with spinner |
| `translate_error(msg)` | Heuristic mapping from R error messages to user-friendly text |
| `js_string(x)` | Safely-escape a string for inclusion in inline JS |
| `format_df_for_display(df)` | Round numerics, format p-values, prepare a data.frame for DT |
| `format_pval(x)` | "p = 0.013" style formatting with scientific notation for small p |
| `build_color_vector(input, prefix, groups)` | Reads per-group color pickers; returns named character vector |
| `apply_plot_labels(p, title, subtitle, ...)` | Adds title/subtitle/axis labels to a ggplot |
| `apply_legend_pos(p, pos)` | Handles "right" / "top" / "inplot_tl" etc. legend positions |
| `observe_palette_preset(input, session, prefix, groups_rv)` | Wires colorblind-friendly palette preset selector |
| `ensure_preprocessed(rv)` | Memoizes data preprocessing across modules |
| `apply_subset_filter(df, rv)` | F.5 — drops rows by IDs / Treatments / Cages exclusion lists |

### `R/helpers_bayes.R`

| Function | Purpose |
|---|---|
| `bayes_diagnostics_panel(result)` | Renders the unified Bayesian diagnostics block: NUTS banner + Bayes R² + PPC coverage + LOO. Used by every Bayesian module's MCMC Diagnostics tab |
| `bayes_influential_obs(result, threshold)` | Filters LOO Pareto-k values for influential observations; returns a data.frame |

### `R/helpers_report.R`

| Function | Purpose |
|---|---|
| `build_html_report(result, title, settings, file)` | Composes a standalone HTML report from an analysis result list — base-R + manual HTML, no external dependencies. Called from the F.1 Summary-tab download handler. |

### `R/logging.R`

Structured JSON logging via the `logger` package. Every analysis click, every download, and every session start/end emits a JSON line to stdout. Buffered to an in-process ring buffer for the Help tab's "Recent activity" view.

---

## Adding a new analysis module

The conventional shape:

1. Create `R/mod_<new>.R` with `<new>UI(id)` and `<new>_server(id, rv)`.
2. If it's > ~1,500 LOC, split per the C.1/C.3 convention into `mod_<new>_plots.R` and `mod_<new>_downloads.R`.
3. Add a `tabPanel("New Module", <new>UI("new"))` to the navbarPage in `R/app.R`.
4. Add `<new>_server("new", rv)` in the server body.
5. If the new module has a long-running Bayesian fit:
   - Create `R/<new>_resume.R` mirroring `tg_resume.R`. Define a `kind` tag like `"bayes_<new>"`.
   - Wire the Bayesian branch to call `launch_<new>_bayes_job(rv$session_id, params, ...)`.
   - Call `register_<new>_resume(input, output, session, rv)` at the bottom of `<new>_server()`.
   - Add the `kind` to `.kind_to_label()` in `R/resume_banner.R` so the banner labels it correctly.
6. Add tests in `tests/testthat/`:
   - `test-module_servers.R`: an init test asserting the module mounts cleanly.
   - For resumable fits: a `test-<new>_resume.R` asserting the per-kind finalizer populates rv correctly and a registration assertion in `test-resume_e2e.R`.

---

## What's deliberately NOT here

A few patterns intentionally aren't used; they came up during the C.1/C.3 / F.7 work and were rejected.

| Anti-pattern | Why we don't do it |
|---|---|
| `R/<X>/` subdirectories per module | R CMD INSTALL flattens `R/`. Subdirectories don't survive the install step, so the convention is the `_suffix` flat naming instead. |
| Putting plot reactives in `_plots.R` and accessing them via closure from `_downloads.R` | Closure capture across files is invisible to readers. We pass plot reactives explicitly as function arguments — the cross-file contract is then visible at the call site. |
| Read `input$mE_session_id` to get the F.7 sid | Forgeable from any websocket message. The sid is always read server-side from `session$request$HTTP_COOKIE`. The audit test `test-session_isolation.R::"get_sid is the only sid source"` greps every R/ file to assert no exported code reads that input. |
| Inline `if/else` chains across module files for resume dispatch | The `register_resume_handler` / `dispatch_resume` registry decouples the banner from any specific module. Modules register their `kind` at server startup. |
| Persisting the F.7 registry to disk | Explicit design choice: in-memory only. Container restart loses in-flight fits. A SQLite-backed registry would survive but adds a service dependency. Tracked for a future iteration if the need arises. |

---

## Common code paths to know about

| When you need to … | Read this file |
|---|---|
| Add a new column to the data-format spec | `R/mod_data_upload.R` (column detection) + [`DATA_FORMAT.md`](DATA_FORMAT.md) |
| Change what shows in the Bayesian MCMC Diagnostics tab | `R/helpers_bayes.R::bayes_diagnostics_panel()` |
| Add a download button | The module's `_downloads.R` file (or `R/<X>.R` if no `_downloads.R` exists) |
| Surface a new Bayesian diagnostic | Backend `mouseExperiment::bayesian_*` adds the field; `helpers_bayes.R` renders it; module tab insertion observers wire it up |
| Change the resume-banner UI | `R/resume_banner.R::.build_resume_banner_ui()` |
| Fix a sid leak | Grep for `session_id` across `R/`, audit `test-session_isolation.R` |
| Add a new Bayesian fit type | New `R/<new>_resume.R` + new `kind` tag + module integration (see section above) |
| Change the future worker count | `R/app.R::runDashboard()` |
| Change the registry TTL | `R/ttl_sweeper.R::start_ttl_sweeper(max_done_age, max_run_age)` |

---

## See also

- [`docs/DATA_FORMAT.md`](DATA_FORMAT.md) — long-format CSV spec, column detection rules
- [`docs/SECURITY.md`](SECURITY.md) — F.7 isolation guarantees, threat model
- [`docs/DEPLOYMENT.md`](DEPLOYMENT.md) — Docker / NPM / Cloudflare deployment
- [`CONTRIBUTING.md`](../CONTRIBUTING.md) — code-style + testing conventions
- Backend `docs/METHODS.md` — statistical methods reference
- Backend `docs/BAYESIAN.md` — Bayesian diagnostics interpretation guide
- `CODE_REVIEW.md` (repo root) — the 25-item review tracking document; explains why a lot of the current structure exists
