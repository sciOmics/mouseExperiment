# Contributing

How to work on `mouseExperiment`. This document covers the conventions, testing requirements, and review process. For the statistical methods reference, read [`docs/METHODS.md`](docs/METHODS.md); for the Bayesian-specific guide, [`docs/BAYESIAN.md`](docs/BAYESIAN.md).

---

## Setup

```r
# Clone, then in R:
setwd("/path/to/mouseExperiment")
devtools::load_all()

# Run the test suite
devtools::test()

# Regenerate man/ from roxygen
devtools::document()

# Full R CMD check
devtools::check()
```

For Bayesian feature development, you also need `brms` ≥ 2.19. For the `cmdstanr` backend, also `cmdstanr` + a working CmdStan toolchain (`cmdstanr::install_cmdstan()`).

---

## Branch / commit conventions

Two long-lived branches:

| Branch | Purpose |
|---|---|
| `staging` | Active development. Dashboard's staging environment installs from this branch |
| `production` | Production deploys. Dashboard's production environment installs from this branch |

Feature branches: `feature/<short-description>`, merged into `staging` first, promoted to `production` once verified.

### Commit message style

Conventional-commits format: single-line summary, empty line, body explaining the *why*. Body wraps at ~72 columns.

```
feat: close D.2, D.3, K.10 — config helpers + plots/data split + covr (v0.4.7)

Closes the three long-deferred architecture / process items from
CODE_REVIEW.md Round 2.

D.2 — config helpers (additive, no breaking change):
  * R/bayesian_config.R exports `tg_priors()` and `tg_mcmc()` —
    bundle the five prior arguments and four MCMC arguments accepted
    by every bayesian_* entry point.
```

Types:

| Type | Use for |
|---|---|
| `feat:` | New exported function or substantial behavior change |
| `fix:` | Bug fix |
| `refactor:` | Internal restructuring, no behavior change |
| `test:` | Test-only changes |
| `docs:` | Documentation only |
| `chore:` | Version bumps, CHANGELOG, CODE_REVIEW updates |

Reference the CODE_REVIEW.md item ID in the body when applicable.

### Co-author trailer

Every commit message ends with:

```
Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
```

when AI assistance was used.

---

## File organization

The package's R/ directory follows a per-concept layout. Each statistical function lives in its own file, and complex multi-path analyses (like `tumor_growth_statistics()`) split helper functions into separate files prefixed with the parent's abbreviation.

| Pattern | Examples |
|---|---|
| Top-level user function | `tumor_growth_statistics.R`, `bayesian_tumor_growth.R`, `survival_statistics.R` |
| Helpers for a specific top-level function | `tgs_path_auc.R`, `tgs_gam.R` (helpers for `tumor_growth_statistics`) |
| Cross-cutting utilities | `utils_bayes.R`, `utils_synergy.R`, `utils_auc.R`, `utils_necrosis.R` |
| Plot-only files | `plot_tumor_growth.R`, `plot_auc.R`, `plot_drug_synergy.R`, etc. |
| Config helpers | `bayesian_config.R` (CODE_REVIEW.md D.2) |
| Result class | `me_result.R` |

### When to split

| Trigger | Action |
|---|---|
| Top-level function > 600 LOC | Extract paths into `<abbrev>_path_*.R` helpers (e.g., `tgs_path_auc.R`) |
| Top-level function > 1,000 LOC | Definitely split |
| Shared utility used by ≥ 3 callers | Extract to `utils_*.R` |
| Plot generation > 100 LOC inside an analysis function | Extract to `plot_*.R` (CODE_REVIEW.md D.3) |

### The data/plot separation (CODE_REVIEW.md D.3)

Every analysis function takes `plots = TRUE` / `FALSE`. When `FALSE`:
- The returned list has all data frames populated
- Plot fields (e.g., `pp_check_plot`) are `NULL`
- The function body skips all ggplot generation, importing no ggplot2 in that call path

Maintain this invariant. New analysis functions must respect `plots = FALSE`. New diagnostic plots must be guarded by `if (isTRUE(plots))` blocks.

---

## Code style

### Conventions

- **snake_case** for function and variable names: `tumor_growth_statistics`, `n_chains`. Never camelCase
- **Roxygen for every exported function**, including `@param` for every argument, `@return`, `@examples` for non-trivial functions, and `@references` for published methods
- **`@keywords internal`** for non-exported helpers; they still get roxygen docstrings
- **No `eval(parse(text = user_input))`** — security and clarity
- **Match the existing style of the file you're editing.** The codebase has some files with strict 80-column compliance and others where ggplot constructions overflow. Don't reformat for style alone

### Bayesian conventions

When adding a new Bayesian function:

1. **Accept `priors = NULL, mcmc = NULL`** parameters and resolve via `.resolve_priors()` / `.resolve_mcmc()` (in `R/bayesian_config.R`). This keeps all `bayesian_*` entry points uniformly callable via the config helpers
2. **Surface the full diagnostics suite:** every Bayesian fit returns `mcmc_diagnostics`, `nuts_diagnostics`, `bayes_R2`, `ppc_coverage`, `loo_diagnostics` (when computable). The dashboard's `bayes_diagnostics_panel()` consumes these
3. **Respect `plots = FALSE`** — guard plot generation with `if (isTRUE(plots))`
4. **Return a flat list with documented fields.** Don't nest result objects more than one level deep unless there's a strong reason (sub-results for AUC + LMM paths is the established pattern)
5. **`return_model = TRUE`** by default so the dashboard can re-extract diagnostics later. Document that `return_model = FALSE` reduces memory at the cost of post-hoc analysis

### Frequentist conventions

- Match the existing return-list shape (`treatment_effects`, `pairwise_comparisons`, etc.) so the dashboard's per-module table renderers work without per-module customization
- Use `emmeans` for marginal means (consistent treatment of contrasts across the package)
- Surface `data_summary`, `growth_rates`, `diag_qq_plot`, `diag_resid_fitted_plot` when applicable

---

## Testing

Run the suite:

```r
devtools::test()                                    # full suite, ~300 tests
devtools::test(filter = "tumor_growth")             # TG-related only
devtools::test(filter = "bayesian_tumor_growth")    # specific Bayesian path
```

Bayesian tests skip when `brms` isn't installed. This is intentional — `R CMD check` should pass even on minimal R installs.

### Test files

Roughly one test file per analysis function, plus:

| File | Covers |
|---|---|
| `test-additional_functions.R` | Cross-cutting helpers |
| `test-toxicity_functions.R` | BW, AUC, TWM, ET bivariate, weight-corrected TGI |
| `test-plot_functions.R` | Plot helpers |
| `test-utils_and_me_result.R` | Utilities + `me_result` S3 class |
| `test-param_sensitivity.R` | "Parameter actually changes output" tests (CODE_REVIEW.md K.4) |
| `helper-fixtures.R` | Shared fixtures (`make_bw_simple()`, etc.) |

### When adding a new function

Add a `test-<function>.R` file with:

1. **Smoke test:** the function runs on the canonical fixture without erroring
2. **Output schema test:** the returned list has the documented fields with the documented types
3. **Edge cases:** small N, single group, NA handling, unusual arg combinations
4. **Parameter sensitivity:** changing a key argument actually changes the output (CODE_REVIEW.md K.4 pattern)

For Bayesian functions, also:

5. **`brms` skip:** wrap with `skip_if_not_installed("brms")` so the test is non-fatal without brms
6. **Diagnostics surface:** assert all expected diagnostic fields are populated
7. **`plots = FALSE`** test: assert plot fields are `NULL` and the function still completes

### Coverage measurement (CODE_REVIEW.md K.10)

```bash
Rscript coverage.R                    # per-file summary
Rscript coverage.R --html             # HTML report
```

Known caveat: pre-existing stale tests in `test-post_power_analysis.R` (function removed in v0.3.4) and `test-toxicity_functions.R` (`efficacy_metric = "log_cell_kill"` no longer valid) block `covr` from completing. These are tracked under CODE_REVIEW.md K.11 — fix them first to get a clean coverage baseline.

---

## Documentation

### Roxygen

Required for every exported function. Required content:

- `@title`, `@description` (the first paragraph is the title; the rest is description)
- `@param` for every argument, including defaults and accepted values
- `@return` describing the list shape (or single object)
- `@examples` for non-trivial functions (use `\dontrun{}` for examples that need brms / Stan)
- `@references` for published statistical methods (DOI when available)
- `@seealso` linking to companions (e.g., the frequentist version)
- `@section`s for substantive topics (Assumptions and Limitations, MCMC diagnostics, etc.)

### NAMESPACE

Regenerate via `devtools::document()`. Don't edit `NAMESPACE` by hand. The legacy `build.R` script in the repo root is *not* the source of truth — it's a vestige and is bypassed by `devtools::document()`.

### Vignettes

Three vignettes are maintained in `vignettes/`:

| File | Topic |
|---|---|
| `mouseExperiment.Rmd` | Package overview |
| `mouseExperiment_combo_demo.qmd` | Worked combination-treatment example |
| `mouseExperiment_dose_demo.qmd` | Worked dose-response example |

Build the Rmd with `devtools::build_vignettes()`. Build the qmd files with `quarto render`.

### `docs/` directory

| File | Topic |
|---|---|
| [`docs/METHODS.md`](docs/METHODS.md) | Statistical methods reference; per-module |
| [`docs/BAYESIAN.md`](docs/BAYESIAN.md) | Priors, diagnostics, MCMC interpretation |

Keep these in sync with the code. When you add a new diagnostic field, update BAYESIAN.md. When you add a new analysis function, update METHODS.md.

### CHANGELOG.md

Follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/). Every version bump gets an entry. Substantial features get an entry that explains *what changed* and *why*.

### CODE_REVIEW.md

`CODE_REVIEW.md` (repo root) tracks the multi-round code review. When you close an item, mark it `✅ Fixed v<version>` with a one-line summary and update the status table.

---

## Pull request process

1. Branch off `staging`
2. Make changes; full test suite passes locally
3. Run `devtools::document()` and commit any updated `man/*.Rd` files
4. Open a PR with:
   - Title: conventional-commits format
   - Body: explains the *why*; references CODE_REVIEW.md item ID if applicable
   - For multi-commit PRs, the body summarizes the sequence
5. CI runs (test suite + `R CMD check`)
6. Self-review:
   - Tests cover the new behavior
   - Roxygen `@param`, `@return`, `@examples` updated
   - For Bayesian: diagnostics surface intact, `plots = FALSE` honored
   - Vignettes still build if you changed an example used in them
7. Merge to `staging`
8. Verify the dashboard installs and runs against the staging branch (it builds on every push)
9. Open a follow-up PR `staging` → `production` for the release

---

## Releases

A new release is a version bump in `DESCRIPTION` + a CHANGELOG entry + a tag.

```r
# Bump version (e.g., 0.4.7 → 0.4.8)
# Edit DESCRIPTION manually
# Add CHANGELOG.md entry

# Document + check
devtools::document()
devtools::check()

# Tag and push
system("git tag -a v0.4.8 -m 'Release v0.4.8'")
system("git push --tags")
```

The dashboard's `DESCRIPTION` pins backend version via `mouseExperiment (>= X.Y.Z)`. After a release, update that constraint in the dashboard repo if any new exports are required.

---

## When you're stuck

| Question | Where to look |
|---|---|
| What does method X do? | [`docs/METHODS.md`](docs/METHODS.md) |
| What does Bayesian diagnostic Y mean? | [`docs/BAYESIAN.md`](docs/BAYESIAN.md) |
| Why is this organized this way? | `CODE_REVIEW.md` |
| Function arg reference | Roxygen / `?function_name` |
| Worked example | `vignettes/` |
| Coverage gap | `Rscript coverage.R` then read the gaps |

---

## See also

- [`docs/METHODS.md`](docs/METHODS.md) — statistical methods reference
- [`docs/BAYESIAN.md`](docs/BAYESIAN.md) — Bayesian diagnostics + priors
- `CODE_REVIEW.md` — design-decision history
- Dashboard `CONTRIBUTING.md` — dashboard-side conventions
- Dashboard `docs/ARCHITECTURE.md` — how the dashboard exposes these methods
