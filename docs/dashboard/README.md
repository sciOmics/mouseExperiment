# mouseExperimentDashboard

An interactive Shiny web application for analyzing mouse tumor growth experiments. Provides a point-and-click interface to the `mouseExperiment` R package, covering the full analysis pipeline from raw data upload through publication-ready figures.

**Note**: This is a proprietary application for private use. See [LICENSE](LICENSE) for details.

## Overview

`mouseExperimentDashboard` wraps the statistical and visualization capabilities of `mouseExperiment` in a Shiny web application. Users can:

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
| Version | 26.09.001 |
| Backend (`mouseExperiment`) requirement | ≥ 0.23.0 (the release that removed the analyses this dashboard no longer offers) |
| `CODE_REVIEW.md` (Round 1 + Round 2) | 24/25 dashboard items closed — only **I.6** (VPS-side CI smoke test) remains |
| Test suite | 210/210 testthat, 1 shinytest2 (Data tab); 75 F.7-specific isolation/dispatch assertions |
| Deployment | Docker + Nginx Proxy Manager (NPM) on a self-hosted VPS (see `mouseExperimentVPS` repo) |

## Local installation & development

### Prerequisites

- **R** ≥ 4.1
- **mouseExperiment** ≥ 0.23.0 (see [mouseExperiment](https://github.com/sciOmics/mouseExperiment))
- For Bayesian analysis: `brms` ≥ 2.19, `future` ≥ 1.33, `promises` ≥ 1.2, `bayesplot` ≥ 1.10

### Setup

```r
setwd("/path/to/mouseExperimentDashboard")
devtools::load_all()
runDashboard()
```

The app opens at `http://localhost:3838`. Note that on plain `http://localhost` the session-identity cookie can't be `Secure`, so F.7 resume falls back to the legacy inline path — fits still run, just without disconnect-survival.

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
| `R/logging.R` | Structured JSON logging to stdout + in-process buffer |
| `R/ui_modules.R` | Shared UI components / layout helpers |
| `inst/www/css/` | Custom styles (light + dark modes) |
| `inst/www/js/` | `session-id.js` (F.7 cookie), `navbar.js`, `darkmode.js`, `overlay.js`, `umami-events.js` |
| `inst/sample_data/` | Example datasets |
| `tests/testthat/` | 8 test files; 210 testthat tests + 1 shinytest2 |

## Dependencies

| Package | Min version | Purpose |
|---------|-------------|---------|
| R | 4.0 | Language runtime |
| shiny | 1.7.0 | Web framework |
| mouseExperiment | 0.23.0 | Statistical backend |
| ggplot2 | 3.4.0 | Visualization |
| dplyr, tidyr | 1.1, 1.3 | Data manipulation |
| emmeans | 1.8.0 | Estimated marginal means |
| DT | 0.28 | Interactive tables |
| shinyjs, shinycssloaders, shinyFeedback | recent | UI helpers |
| logger | 0.3 | Structured logging |
| jsonlite | 1.8 | Logger JSON serialization |
| brms (Suggests) | 2.19 | Bayesian models |
| future, promises (Suggests) | 1.33, 1.2 | Async execution |
| bayesplot (Suggests) | 1.10 | MCMC diagnostic plots |
| shinytest2, chromote (Suggests) | 0.2, recent | Headless-browser integration tests |

See `DESCRIPTION` for the complete list.

## Deployment

Production deployment is handled by the companion [mouseExperimentVPS](../mouseExperimentVPS) repository — Docker + Nginx Proxy Manager (NPM) on a self-hosted VPS. Key NPM settings the dashboard depends on:
- Websocket support enabled (required for Shiny's reactivity layer)
- Proxy read / send timeouts ≥ 15 minutes (so Bayesian fits don't trigger the NPM disconnect overlay; F.7 covers the gap)
- HTTPS (Let's Encrypt) so `Secure` cookies work for F.7

## License

**Proprietary** — All Rights Reserved. See [LICENSE](LICENSE).

`mouseExperiment` (the underlying R package) is licensed under the MIT License.

## Citation

If you use the underlying `mouseExperiment` package in your research, please cite it from there. `mouseExperimentDashboard` is a proprietary private application and is not for citation or redistribution.

## Support

For questions or issues, contact: david@sciomics.bio
