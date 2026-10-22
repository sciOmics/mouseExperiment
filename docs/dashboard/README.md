# The mouseExperiment dashboard

The interactive Shiny dashboard in the `mouseExperiment` package, started with `run_dashboard()`. It is a point-and-click interface to the package's analysis functions, covering the full analysis pipeline from raw data upload through publication-ready figures.

This folder is the dashboard's developer documentation. Until v0.34.0 the dashboard was a separate private package, `mouseExperimentDashboard`. Its history is now part of this repository, and its changelog and review log are kept here as [CHANGELOG.md](CHANGELOG.md) and [CODE_REVIEW.md](CODE_REVIEW.md).

## Overview

The dashboard wraps the statistical and visualization capabilities of `mouseExperiment` in a Shiny web application. Users can:

- Upload tumor growth data (CSV) and map columns interactively
- Analyze tumor growth with LME4 mixed-effects models, AUC, or Bayesian LMM (brms)
- Run survival analysis with Kaplan-Meier curves, Cox/AFT models, and Bayesian AFT
- Analyze drug synergy via Bliss Independence, with bootstrap intervals
- Fit dose-response curves (Hill/Emax)
- Assess toxicity via a body-weight mixed-effects model and time to a weight-loss threshold, with a Therapeutic Window mode
- Compute a priori power via analytic formulas or LMM simulation
- Customize and download all plots and tables (PNG / CSV / HTML report / ZIP)
- **Resume Bayesian fits after tab disconnect** — the registry survives an NPM proxy timeout or browser sleep; the banner on reconnect offers Resume / Cancel

The Bayesian fits (tumor growth and survival) run asynchronously via `future` / `promises`, so the app remains responsive during Stan compilation and MCMC sampling.

## Current status

| Item | State |
|---|---|
| Version | the package's (see `DESCRIPTION`); the About tab shows it |
| Tests | with the package's, in `tests/testthat/`; some run in a headless browser (shinytest2) |
| Deployment | Docker + Nginx Proxy Manager (NPM) on a self-hosted VPS (see `sciOmicsVPS` repo) |

## Local installation & development

### Prerequisites

- **R** >= 4.1
- The dashboard's own packages (DT, shinyjs, future, promises and others) are Suggests of `mouseExperiment`, so installing the package for its analysis functions does not pull them in. `run_dashboard()` checks for them and, in an interactive session, offers to install any that are missing.

### Setup

```r
# Installed
mouseExperiment::run_dashboard()

# From a clone of the repository
devtools::load_all()
run_dashboard()
```

The app opens at `http://localhost:3838`. Browsers treat `http://localhost` and `http://127.0.0.1` as secure contexts, so the `Secure` session-identity cookie is usually set even without HTTPS. Where it is not (cookies blocked, or another host over plain HTTP), the session gets a local registry key: Bayesian fits run through the job registry as usual but cannot be resumed after a disconnect. Before v26.09.002 that case took a legacy inline path that failed before the fit started (CODE_REVIEW.md R20.D3).

### Test suite

```r
devtools::test()       # full suite
devtools::test(filter = "session_isolation")   # the F.7 security tests
devtools::check()      # R CMD check
```

## Usage

### Basic workflow

1. **Data Upload** — Select a CSV file in long format
2. **Column Mapping** — Map ID, time, volume, treatment, and optional columns; preview validates
3. **Analysis Setup** — Pick a module, choose method (frequentist / Bayesian), set sidebar parameters
4. **Run Analysis** — Click Run; Bayesian fits run in the background and surface in the resume banner if you disconnect
5. **Explore Results** — Per-module result tabs: summary, model output, plots, pairwise comparisons, MCMC diagnostics (Bayesian)
6. **Download** — Per-tab PNG / CSV / "Download All Results" ZIP. HTML report from the Tumor Growth Summary tab

### Analysis modules

| Module | Methods |
|--------|---------|
| **Tumor Growth** | LME4 mixed model, AUC, Bayesian LMM (brms) |
| **Survival** | Kaplan-Meier, Cox PH, log-rank, `cox.zph()` proportional-hazards check, Bayesian AFT (Weibull / log-normal), concordance / C-index |
| **Toxicity** | Body-weight LME4 mixed model; time to a weight-loss threshold (Kaplan-Meier, competing risks flagged); Therapeutic Window Metric |
| **Drug Synergy** | Bliss Independence with bootstrap intervals; synergy over time |
| **Dose-Response** | Hill/Emax with per-dose centre + error overlays |
| **Power Analysis** | A priori analytic (t-test / ANOVA), LMM simulation |

### Resumable Bayesian fits (F.7)

When a Bayesian fit is running and the user's tab disconnects:

1. The future continues running in its worker (separate from the user's Shiny session).
2. The user's `mE_session` cookie (256-bit CSPRNG, 24-hour expiry) re-identifies them on reconnect.
3. On reconnect, a banner above the navbar lists each outstanding job with kind label, elapsed time, status, and Resume / Cancel buttons.
4. Clicking Resume switches to the matching module tab and populates the result reactively (same UI surface as if they never disconnected).
5. Completed jobs are evicted from the registry 30 minutes after they finish; abandoned running jobs after 4 hours.

Multi-user isolation is enforced at four layers:
- Registry keyed by sid (cross-sid lookup returns NULL — the registry has no `list_all_jobs()`)
- Sid is read server-side from the HTTP Cookie header, never from `input$mE_session_id`
- 256-bit cookie entropy + sid is never logged, echoed, or rendered in HTML
- Worker processes don't `eval()` user-supplied strings, so worker memory isn't a cross-job leak vector

See `R/job_registry.R` and `tests/testthat/test-session_isolation.R` for the full security audit.

## Data format

Long format: one row per observation/measurement.

### Required columns

| Column | Notes |
|---|---|
| ID / Mouse_ID / Animal_ID | Unique identifier per animal |
| Day / Time / Timepoint | Study day or date — date detection auto-converts |
| Volume **or** Length + Width | Tumor measurement(s) |
| Treatment / Group | Treatment group assignment |

### Optional columns

| Column | Used by |
|---|---|
| Cage / Cage_ID | LME4 / Bayesian random intercept |
| Survival_Censor / Event | Survival analysis |
| Dose | Dose-response |
| Weight / Mass | Toxicity (body weight) |
| Necrotic | Necrosis handling (exclude / covariate / none) |

### Tumor volume from length + width

| Formula | Equation |
|---------|----------|
| Ellipsoid (default) | V = (L × W² × π) / 6 |
| Modified Ellipsoid | V = (L × W²) / 2 |
| Cylinder | V = (π × W² × L) / 4 |
| Sphere | V = (π × W³) / 6 |

Dimensions are auto-corrected so the longer measurement is always L.

## Architecture

The dashboard is structured as a Shiny `navbarPage` with one tab per analysis module. Each module is a self-contained `moduleServer()` with its UI + server pair. Long-running Bayesian fits go through a process-level job registry so they survive disconnect.

```
                    Shiny session (per user tab)
                    ┌───────────────────────────────┐
                    │  reactiveValues (rv)           │
                    │    session_id, *_results,      │
                    │    *_job_id, ...               │
                    └───────────────┬────────────────┘
                                    │
       ┌────────────────────────────┼────────────────────────────┐
       │ Per-module moduleServer    │  Resume banner             │
       │  Run button →              │   uiOutput("banner")       │
       │  launch_<X>_bayes_job()    │   if list_jobs(sid) > 0    │
       │  rv$<X>_job_id ← job_id    │     show one row per job   │
       │  attach_bayes_poller       │   on Resume click →        │
       │   (polls registry every    │     dispatch_resume(job)   │
       │    2 s; on terminal status │                            │
       │    → finalize_<X>_success) │                            │
       └────────────┬───────────────┴────────────────────────────┘
                    │
                    ▼
          Process-level .job_registry env
          ┌────────────────────────────────────┐
          │  sessions: list(                    │
          │    "<sid>" = list(                  │
          │      "<job_id>" = list(             │
          │        kind, status, result, error, │
          │        started_at, completed_at,    │
          │        context = list(...))))       │
          └─────────────────┬───────────────────┘
                            │
                            ▼
                  future workers (multisession)
                  ↑ run to completion regardless
                    of which session started them
```

### Key files

| File | Purpose |
|------|---------|
| `R/app.R` | App entry point; UI layout, server dispatch, async plan setup, banner mount, sweeper start |
| `R/mod_data_upload.R` | Data upload, column mapping, preprocessing cache |
| `R/mod_tumor_growth.R` | TG analysis UI + server (LME4, AUC, Bayesian) |
| `R/mod_tumor_growth_plots.R` | 5 TG plot reactives + their `renderPlot` bindings + diagnostics + group-order UI |
| `R/mod_tumor_growth_downloads.R` | 16 TG download handlers + necrosis renderDT |
| `R/mod_survival.R`, `mod_survival_downloads.R` | Survival analysis (KM, Cox, Bayesian AFT) |
| `R/mod_toxicity.R`, `mod_toxicity_plots.R`, `mod_toxicity_downloads.R` | Body weight (LMM, weight-loss threshold) and therapeutic window |
| `R/mod_synergy.R`, `mod_synergy_downloads.R` | Drug synergy |
| `R/mod_dose_response.R` | Dose-response (Hill/Emax) |
| `R/mod_power.R` | Power analysis |
| `R/mod_data_upload.R`, `R/mod_log_session.R` | Cross-cutting modules (upload; log, filters and settings presets) |
| `R/job_registry.R` | F.7 process-level job registry |
| `R/session_identity.R` | F.7 server-side cookie reader |
| `R/resume_banner.R` | F.7 banner UI/server + resume dispatch |
| `R/bayes_job_runner.R` | F.7 shared poller + resume registration helpers |
| `R/tg_resume.R`, `surv_resume.R` | F.7 resume helpers for the two Bayesian fits (tumor growth, survival) |
| `R/ttl_sweeper.R` | F.7 periodic registry housekeeping |
| `R/helpers.R` | Core utilities: validation, error translation, formatting, plot helpers |
| `R/helpers_bayes.R` | Bayesian diagnostics panel + LOO influential observations |
| `R/helpers_report.R` | Standalone HTML report builder |
| `R/logging.R` | Structured JSON logging to stdout |
| `R/ui_modules.R` | Shared UI components / layout helpers |
| `inst/www/css/` | Custom styles (light + dark modes) |
| `inst/www/js/` | `session-id.js` (F.7 cookie), `navbar.js`, `darkmode.js`, `overlay.js`, `umami-events.js` |
| `inst/sample_data/` | Example datasets |
| `tests/testthat/` | The package's tests, the dashboard's among them; the shinytest2 browser tests use the apps in `tests/testthat/apps/` |

## Dependencies

`shiny` is an Import of `mouseExperiment`. The dashboard's other packages are Suggests: `dashboard_packages()` in `R/app.R` lists them, and `run_dashboard()` checks them. The Bayesian fits use the package's own `brms` stack. See `DESCRIPTION` for versions.

## Deployment

Production deployment is handled by the companion [sciOmicsVPS](https://github.com/sciOmics/sciOmicsVPS) repository — Docker + Nginx Proxy Manager (NPM) on a self-hosted VPS. Key NPM settings the dashboard depends on:
- Websocket support enabled (required for Shiny's reactivity layer)
- Proxy read / send timeouts ≥ 15 minutes (so Bayesian fits don't trigger the NPM disconnect overlay; F.7 covers the gap)
- HTTPS (Let's Encrypt) so `Secure` cookies work for F.7

## License

MIT, as for the rest of `mouseExperiment`. See [LICENSE](../../LICENSE).

## Citation

If you use the dashboard in your research, please cite the package: `citation("mouseExperiment")`.

## Support

For questions or issues, contact: david@sciomics.bio
