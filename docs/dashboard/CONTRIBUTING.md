# Contributing

How to work on `mouseExperimentDashboard`. This document covers the conventions, testing requirements, and review process. For the technical layout, read [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) first.

---

## Setup

```r
# Clone, then in R:
setwd("/path/to/mouseExperimentDashboard")
devtools::load_all()

# The backend is a hard dependency; DESCRIPTION's Imports gives the minimum
# version (mouseExperiment >= 0.30.0 as of 26.10.002)
# If you're developing both, sync the local backend:
devtools::install("/path/to/mouseExperiment", quiet = TRUE, upgrade = FALSE)

# Run the dashboard
runDashboard()
```

For Bayesian feature development you also need `brms` ≥ 2.19, `future`, `promises`, and `bayesplot`.

---

## Branch / commit conventions

The repository follows trunk-based development with two long-lived branches:

| Branch | Purpose |
|---|---|
| `staging` | Active development. Push freely. Triggers a webhook deploy to `https://staging.mouseexperiment.com` |
| `production` | Production deploys only. Triggers a webhook deploy to `https://mouseexperiment.com` |

Feature branches: `feature/<short-description>`, merged into `staging` first, then promoted to `production` via a fast-forward merge once verified.

### Commit message style

Single-line summary in conventional-commits format, followed by an empty line and a body explaining the *why* (not the *what* — the diff shows the what). Body wraps at ~72 columns.

```
feat: F.7 — job registry + session-identity foundation

First of six commits implementing CODE_REVIEW.md F.7 (session resume).
Lays the security-critical foundation: a process-level job registry
that holds long-running futures outside the per-session reactive scope...
```

| Type | Use for |
|---|---|
| `feat:` | New user-visible behavior |
| `feat(<mod>):` | New behavior in a specific module (`feat(tg):`, `feat(dr):`) |
| `refactor:` | Internal restructuring, no behavior change |
| `fix:` | Bug fix |
| `test:` | Test-only changes |
| `docs:` | Documentation only |
| `chore:` | Version bumps, CHANGELOG, CODE_REVIEW updates |

Reference the CODE_REVIEW.md item ID in the body when applicable: "Closes F.7", "Partial fix for C.1".

Every commit message ends with the trailer:

```
Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
```

when AI assistance was used in writing the commit.

### Atomic commits

A commit should leave the test suite passing. Don't push WIP commits to shared branches. Use a feature branch and squash before merge if you need intermediate breakpoints.

For multi-commit feature work (e.g., a six-commit F.7 implementation), each commit should be **independently functional** — the tests pass at every point in the history.

---

## File organization

The C.1/C.3 refactor established a per-concern naming convention. Adhere to it when adding or modifying modules.

| Suffix | What lives there |
|---|---|
| `mod_<X>.R` | Lifecycle: observers, fit dispatch, status messages, tab insertion, table renderers |
| `mod_<X>_plots.R` | Plot reactives + their `renderPlot` bindings + group-order/palette UI + diagnostic plots |
| `mod_<X>_downloads.R` | All `downloadHandler()` bindings (CSV / PNG / ZIP) |

A module file should not exceed ~1,500 LOC. If it does, split per the suffix convention.

### When to split

Trigger | Action
---|---
File > 1,500 LOC | Extract `_plots.R` (typically the largest chunk)
File > 2,500 LOC | Extract `_downloads.R` as well
Multiple unrelated download handlers | Extract `_downloads.R` even at smaller file sizes

### Adding a new analysis module

See [`docs/ARCHITECTURE.md::"Adding a new analysis module"`](docs/ARCHITECTURE.md#adding-a-new-analysis-module).

### F.7-specific conventions

If your new module has a long-running Bayesian fit, route it through the F.7 registry. Don't write a new ad-hoc `promises::future_promise({...}) %...>% (function(result) {...})` chain.

Required steps for a resumable Bayesian fit:

1. Create `R/<X>_resume.R` mirroring `R/tg_resume.R`:
   - `launch_<X>_bayes_job(sid, params, input, output, session, rv)`
   - `register_<X>_resume(input, output, session, rv)`
   - `.<X>_bayes_promise(p)`, `.finalize_<X>_success`, `.finalize_<X>_failure`, `.finalize_<X>_cancellation`
   - no separate cookie-less path: `R/app.R` gives such sessions a `make_local_sid()` key, so the registry path covers them (R20.D3)
2. Wire your module's Bayesian branch to call `launch_<X>_bayes_job(rv$session_id, params, ...)` instead of the inline future.
3. Call `register_<X>_resume(input, output, session, rv)` at the bottom of `<X>_server()`.
4. Add the `kind` tag to `.kind_to_label()` in `R/resume_banner.R`.
5. Add a test in `tests/testthat/test-resume_e2e.R` asserting `register_<X>_resume` exists and is callable.

---

## Code style

The codebase is intentionally permissive on lint warnings — many existing files have lines over 80 characters and non-standard hanging indents because matching ggplot constructions readably wins over strict 80-column compliance. Don't reformat an existing file purely for lint cleanup.

### Rules we enforce by convention

- **snake_case** for function and variable names: `launch_tg_bayes_job`, `rv$session_id`. Never camelCase except where it matches Shiny conventions (`tabPanel`, `actionButton`).
- **Roxygen for every exported function.** Non-exported helpers may use a minimal `@keywords internal` block.
- **No new emoji in code or comments** unless a UI element explicitly calls for one (e.g., status icons via `shiny::icon()`).
- **No `eval(parse(text = user_input))`** anywhere. This is non-negotiable — F.7 worker isolation depends on it.
- **Pure helpers go in `helpers.R`** (or `helpers_bayes.R` / `helpers_report.R` if the topic is large enough to warrant its own file).
- **Tests for new behavior.** PRs that add user-visible behavior without tests will be asked for tests before merge.

### When in doubt, match the existing style of the file you're editing.

---

## Testing

**Install the sibling backend first.** `devtools::load_all()` resolves
`Imports: mouseExperiment` against the **installed** library, not
`../mouseExperiment`, so testing against a stale backend is the default failure
mode rather than an unlikely accident:

```r
devtools::install("../mouseExperiment", upgrade = FALSE)
devtools::test()
```

`test-backend-version.R` **fails** — it does not skip — when the installed backend
is older than the sibling source, and tells you the command above. That is
deliberate: during the Round 3 review this suite reported 210 passing against a
backend six releases behind the code under review, and a green run against the
wrong backend is worse than no run. See `CODE_REVIEW.md` §R4.D2.

Targeted runs:

```r
devtools::test(filter = "session_isolation")    # F.7 security audit
devtools::test(filter = "resume_e2e")           # F.7 disconnect/reconnect
devtools::test(filter = "module_servers")       # module init smoke tests
devtools::test(filter = "helpers")              # helper functions
```

### Test files

| File | What it covers |
|---|---|
| `test-helpers.R` | Pure helper functions |
| `test-new-helpers.R` | More recent helpers (subset filter, palette presets, etc.) |
| `test-modal-and-cachebust.R` | Modal display + cache-busting query string |
| `test-module_servers.R` | Module init smoke tests (`shiny::testServer`) |
| `test-session_isolation.R` | F.7 four-layer isolation audit. **Do not weaken these tests** |
| `test-resume_banner.R` | F.7 banner UI/server |
| `test-tg_resume.R` | F.7 TG-specific helpers |
| `test-resume_e2e.R` | F.7 disconnect/reconnect cycle |
| `test-demo_data.R` | Demo dataset integrity |
| `test-shiny_app.R` | shinytest2 headless-browser integration (skips if `chromote` missing) |
| `test-backend-contract.R` | The named fields the dashboard reads from backend results (R3.D9, I.2, R4.D2) |
| `test-backend-calls-resolve.R` | Every `mouseExperiment::` call names a function that exists (R17.3) |
| `test-backend-version.R` | Installed backend vs sibling source skew; **fails**, does not skip (R4.D2) |
| `test-demo-data-reachable.R` | The backend's demo datasets load from the app (D.12) |
| `test-docs-current.R` | Info / Help tabs match the app and backend (D.13) |
| `test-helpers_data.R`, `test-helpers_format.R`, `test-helpers_plot_ui.R` | Data, formatting, plot and UI helpers (D.1, D.9) |
| `test-power-params.R` | Power module parameter plumbing (D.4, K.15) |
| `test-render-logic.R` | Render paths driven through `testServer()` (I.2) |
| `test-round21-fixes.R` | Round 21 fixes (R20.D*), labelled by finding ID |
| `test-round21-study-design.R`, `test-round21-step4.R`, `test-round21-step5.R`, `test-round21-step6.R`, `test-round21-step7.R` | Steps 3–7 of the Round 20/21 implementation; each test fails on the previous version |
| `helper-app.R` | Test fixtures (`make_test_rv()`, etc.) |

### F.7 security tests (specifically)

The four-layer F.7 isolation tests in `test-session_isolation.R` are non-negotiable. They include grep audits that fail if anyone:

- Reads `input$mE_session_id` from any R/ file
- Logs `session_id` via `log_*` / `message` / `stop` / `warning` / `showNotification`
- Includes the sid in any banner output

If you legitimately need to bypass any of these (you almost never do), discuss in the PR before making the change.

### testServer hangs

`shiny::testServer` doesn't play well with `invalidateLater` observers. If you add a periodic polling observer to a module and the test hangs:

- Check whether `invalidateLater` is being scheduled inside a reactive that depends on a `reactiveVal` that the observer updates (recursive loop)
- Either guard with a flag the test can set, or remove the periodic observer (the F.7 banner originally had one; we removed it because the gain was small and the test debt was real)

### Adding new tests

For Shiny modules: use `shiny::testServer()` for in-process testing. Use `shinytest2::AppDriver` for end-to-end browser tests — these are heavier so use them sparingly.

For helpers: plain `testthat::test_that()` blocks.

For F.7 isolation: any change to `R/job_registry.R`, `R/session_identity.R`, or `R/resume_banner.R` requires updating the corresponding `test-*` file.

---

## Documentation

Roxygen docs are required for every exported function. The reference list in `man/` is regenerated by:

```r
devtools::document()
```

User-facing docs live in `docs/`:

| File | What it covers |
|---|---|
| [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) | Module pattern, F.7 flow, async |
| [`docs/DATA_FORMAT.md`](docs/DATA_FORMAT.md) | Long-format CSV spec, column detection |
| [`docs/SECURITY.md`](docs/SECURITY.md) | F.7 threat model, isolation guarantees |
| [`docs/DEPLOYMENT.md`](docs/DEPLOYMENT.md) | Docker / NPM / Cloudflare setup |

Keep these in sync with the code. When you change cookie attributes, update SECURITY.md. When you add a new analysis module, update ARCHITECTURE.md and DATA_FORMAT.md.

### CHANGELOG

`CHANGELOG.md` follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/). Every version bump gets an entry. For substantial features, the entry explains *what changed* and *why*, not just the diff.

### CODE_REVIEW.md

`CODE_REVIEW.md` (repo root) tracks every review round, with item IDs per round (for example D.13, R20.D27). When you close an item:

1. Edit the item's entry to add `✅ Fixed v<version>` with a one-line summary
2. Update the status table at the bottom
3. Reference the item ID in your commit message

---

## Pull request process

1. Branch off `staging`
2. Make changes; tests pass locally
3. Open a PR with:
   - Title: matches the commit-message convention (`feat:`, `fix:`, etc.)
   - Body: same shape as the commit body (the *why*); for multi-commit PRs, summarize the sequence
   - Linked CODE_REVIEW.md item ID if applicable
4. CI runs (test suite + R CMD check)
5. Self-review for:
   - Tests cover new behavior
   - No code or comments leak the F.7 sid
   - No new lines longer than ~120 chars unless matching pattern in the file
   - Docs updated for any user-visible change
6. Merge to `staging` (rebase or squash; we don't strictly enforce one over the other)
7. Verify on staging.mouseexperiment.com
8. Open a follow-up PR `staging` → `production` for the deploy

---

## When you're stuck

| Question | Where to look |
|---|---|
| How does the F.7 banner work? | `R/resume_banner.R` + [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) + [`docs/SECURITY.md`](docs/SECURITY.md) |
| What columns does the dashboard accept? | [`docs/DATA_FORMAT.md`](docs/DATA_FORMAT.md) |
| What does the Bayesian X mean? | Backend `docs/BAYESIAN.md` and `docs/METHODS.md` |
| Why is X organized this way? | `CODE_REVIEW.md` (the relevant item explains the history) |
| Deployment broke | `sciOmicsVPS/docs/troubleshooting.md` |
| A test is failing | Read the test + the code path it exercises; many tests have detailed comments |

---

## See also

- [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md)
- [`docs/DATA_FORMAT.md`](docs/DATA_FORMAT.md)
- [`docs/SECURITY.md`](docs/SECURITY.md)
- [`docs/DEPLOYMENT.md`](docs/DEPLOYMENT.md)
- `CODE_REVIEW.md` (repo root) — design-decision history
- Backend `CONTRIBUTING.md` — for backend-side conventions
