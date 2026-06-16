# Model diagnostics

What the package's analyses check, what the dashboard renders, what's appropriate for each method, and where the gaps are. Pairs with [`METHODS.md`](METHODS.md) (when to choose a method) and [`BAYESIAN.md`](BAYESIAN.md) (how to interpret the Bayesian diagnostics suite in detail).

This document is part audit, part reference. Sections 1-5 inventory the diagnostics per analysis; sections 6 and 7 cover what's missing and the known bugs / gaps in the current surfacing.

---

## TL;DR

- **Bayesian path:** comprehensive and well-rendered. Rhat / ESS / NUTS / Bayes R² / PPC coverage / PSIS-LOO + Pareto-k + influential observations + PPC density overlay + prior-posterior overlay + MCMC trace. Unified `bayes_diagnostics_panel()` helper renders them consistently across modules.
- **Cox survival:** PH assumption (`cox.zph` per-term + global) and concordance / C-index. Both surfaced in dedicated tabs.
- **Frequentist Body Weight (`analyze_body_weight`):** Q-Q plot + residuals-vs-fitted ggplot. Rendered properly.
- **Frequentist Tumor Growth LME4:** **the dashboard's TG Diagnostics tab is broken** — the backend doesn't produce the fields the dashboard tries to render. See section 7.
- **GAM / GAMM paths:** no rendered diagnostics. `mgcv::k.check()` and `mgcv::concurvity()` are the missing essentials.
- **Frequentist Dose-Response:** no goodness-of-fit, no residual plots, no Hill-fit diagnostics. `drc` exposes plenty (`modelFit`, `mselect`); none surfaced.
- **Frequentist Synergy:** no model-assumption diagnostics for the t-test path.
- **All LME4 paths:** no random-effects QQ, no influence diagnostics, no homogeneity-of-variance check.
- **Repeated Measures ANOVA:** no diagnostics returned at all.

If you only read one section: [Section 7 — Bugs and gaps to fix](#7-bugs-and-gaps-to-fix).

---

## 1. Bayesian diagnostics (all `bayesian_*` functions)

Every Bayesian fit returns a uniform diagnostics surface. The backend's `R/utils_bayes.R` helpers (`make_nuts_diagnostics()`, `bayes_loo()`, `bayes_ppc_coverage()`, `bayes_r2_summary()`) compute the fields; the dashboard's `R/helpers_bayes.R::bayes_diagnostics_panel()` renders them.

### Convergence — Rhat, ESS

| Field | Type | Where computed | Where rendered |
|---|---|---|---|
| `mcmc_diagnostics$Rhat` | numeric per parameter | `R/utils_bayes.R::make_mcmc_diagnostics()` | `output$<x>_mcmc_diag_table` (DT) |
| `mcmc_diagnostics$Bulk_ESS`, `Tail_ESS` | numeric per parameter | same | same table |
| `mcmc_diagnostics$Converged` | logical per parameter | `Rhat <= 1.01` test | same table |

**Assessment:** the threshold `Rhat <= 1.01` matches modern guidance (Vehtari et al. 2021; older guidance of 1.1 was too lax). ESS thresholds aren't enforced; they're left for the user to read.

**Recommended addition:** an "ESS adequacy" column flagging parameters with `Bulk_ESS < 400` or `Tail_ESS < 400`, mirroring the existing `Converged` column. Currently the user has to eyeball the numbers.

### Sampler — NUTS divergences, max_treedepth, E-BFMI

| Field | Type | Where computed |
|---|---|---|
| `nuts_diagnostics$n_divergent` | integer per chain | `R/utils_bayes.R::make_nuts_diagnostics()` |
| `nuts_diagnostics$n_max_treedepth` | integer per chain | same |
| `nuts_diagnostics$ebfmi_min` | numeric per chain | same |

Rendered by `bayes_diagnostics_panel()` as the colored alert at the top of the MCMC Diagnostics tab: warning when any divergences, max_treedepth hits, or E-BFMI < 0.3; success otherwise.

**Assessment:** correct. All three are the right thresholds (zero divergences ideal; one is a red flag at small N; > 0 always worth investigating). E-BFMI 0.3 is the Stan reference threshold.

### Model fit — Bayes R²

| Field | Type | Where computed |
|---|---|---|
| `bayes_R2$Estimate` | numeric | `R/utils_bayes.R::bayes_r2_summary()` |
| `bayes_R2$Lower_95_CrI`, `Upper_95_CrI` | numeric | same |

Rendered by `bayes_diagnostics_panel()` as a dl-horizontal row.

**Assessment:** correct (uses `brms::bayes_R2()` which uses the Gelman et al. 2019 definition). Always reported with CrI rather than just point estimate — appropriate.

### Predictive coverage — PPC coverage

| Field | Type |
|---|---|
| `ppc_coverage$cov_50`, `cov_80`, `cov_95` | numeric per nominal interval |
| `ppc_coverage$n_obs` | integer |

Rendered by `bayes_diagnostics_panel()` as a dl-horizontal row.

**Assessment:** correct conceptually. Reports empirical coverage rates at 50/80/95% nominal. Off by ≥ 5pp indicates model misspecification.

**Recommended addition:** a small visual — even a bar chart of empirical vs nominal — instead of just a text row. Easier for users to spot under-coverage at a glance.

### Out-of-sample — PSIS-LOO + Pareto-k

| Field | Type |
|---|---|
| `loo_diagnostics$elpd_loo` | numeric |
| `loo_diagnostics$se_elpd` | numeric |
| `loo_diagnostics$n_high_k` | integer |
| `loo_diagnostics$pareto_k` | list-column carrying the per-observation k vector |

Rendered:
- Summary row in `bayes_diagnostics_panel()`
- Influential-observations table via `bayes_influential_obs()` filtered to `k > 0.7`, displayed as a DT below the panel

**Assessment:** comprehensive. The threshold `0.7` matches the Vehtari et al. PSIS guidance.

**Recommended addition:** display the full Pareto-k distribution as a small histogram with the 0.7 cut-off line. Currently users see only the high-k observations; the distribution gives context (e.g., "5 observations > 0.7 but the next-highest is 0.6" vs "5 observations > 0.7 and the next 30 are > 0.65").

### Plots

| Field | What it shows | Where rendered |
|---|---|---|
| `pp_check_plot` | `brms::pp_check(type = "dens_overlay")` — posterior predictive density vs observed | "Posterior Check" tab per Bayesian module |
| `prior_posterior_plot` | Prior (grey) vs posterior (blue) density for each treatment-effect coefficient | "Prior vs Posterior" tab |
| `mcmc_trace_plot` | `bayesplot::mcmc_trace()` for fixed-effect coefficients | "MCMC Diagnostics" tab |

**Assessment:** the three most-useful plots are present. Implementation is straight `bayesplot` / `brms` calls, which is the correct approach.

### Plots produced but NOT rendered

The Bayesian functions also compute and return:

| Field | What it would show |
|---|---|
| `residuals_plot` | Posterior mean residuals vs day, faceted by group, loess smoother — non-linearity check |
| `posterior_dist_plot` | Posterior density areas for treatment parameters (`bayesplot::mcmc_areas`) |
| `credible_intervals_plot` | Forest plot of group EMMs with 95% CrI |

**None of these are rendered in the dashboard.** They're computed (so they cost compile/render time when `plots = TRUE`) and silently dropped. See section 7.

---

## 2. Survival — frequentist (Cox)

### Proportional-hazards check — `cox.zph` (Schoenfeld residuals)

| Field | Type | Where computed |
|---|---|---|
| `ph_test` | `cox.zph` object — per-term + global tests on scaled Schoenfeld residuals | `R/survival_statistics.R:508` |

Rendered:
- Global banner: success when global p ≥ 0.05; warning when p < 0.05 (`output$surv_model_quality_panel`)
- Per-term table in the dedicated PH Assumption tab (`output$surv_ph_test_table`)

**Assessment:** correct method. `cox.zph` is the standard PH check. The banner makes the global result hard to miss.

**Recommended addition:** the **Schoenfeld residuals scatter plot** (`plot(cox.zph_object)`) — visual diagnostic that shows how each covariate's effect changes over time. Currently the user only sees the test table; the plot reveals shape (linear trend in log HR, sinusoidal effect, etc.) that the table can't.

### Concordance / C-index

| Field | Type |
|---|---|
| `c_index$concordance` | numeric — Harrell's C |
| `c_index$var` | numeric — variance (for CI computation) |

Rendered as a dl row in `surv_model_quality_panel` with computed 95% CI (`val ± 1.96*sqrt(var)`).

**Assessment:** correct. Harrell's C is the survival analog of AUC-ROC. CI computation is asymptotic (Wald-style) which is fine at typical sample sizes.

### What's not assessed (and arguably should be)

- **Martingale residuals plot** — diagnoses non-linearity in continuous covariates (e.g., dose). Standard `residuals(model, type = "martingale")`. Not surfaced.
- **Deviance residuals** — identify outliers / high-influence observations. `residuals(model, type = "deviance")`. Not surfaced.
- **`dfbeta`** / `dfbetas` — influence on individual coefficient estimates. `residuals(model, type = "dfbeta")`. Not surfaced.
- **Influence on linear predictor** — `residuals(model, type = "score")`. Not surfaced.
- **Firth correction status** — when `firth_correction = TRUE` triggers because of separation, the user isn't told which covariate had separation. Worth surfacing.

For typical preclinical datasets (N ≤ 30 per group, one categorical covariate = treatment), the missing diagnostics are nice-to-haves rather than essentials. If you ever add continuous covariates (dose, animal age) the martingale plot becomes important.

---

## 3. Survival — Bayesian (parametric AFT via brms)

Inherits the full Bayesian diagnostics suite (section 1). No PH check because the parametric AFT family doesn't assume proportional hazards.

### What's not assessed

- **Survival function calibration** — does the posterior survival function match the Kaplan-Meier estimate? The package returns `survival_curve_plot` which overlays them, so visually checkable. Numeric calibration (e.g., Hosmer-Lemeshow-style) would add rigor but isn't standard for AFT.
- **Shape parameter posterior** — for Weibull/gamma, the shape parameter's posterior tells you whether the hazard is increasing (shape > 1) or decreasing (shape < 1). Currently buried in `posterior_summary`; could be surfaced as a banner like "Posterior median shape = 1.34 → increasing hazard."

---

## 4. Tumor Growth — frequentist

### LME4 path

| Field | Type | Where set |
|---|---|---|
| `diagnostics$residuals$fitted` | numeric vector | `R/tumor_growth_statistics.R:962-966` |
| `diagnostics$residuals$residuals` | numeric vector | same |
| `diagnostics$residuals$qq_plot` | `qqnorm()` output (data, **not a ggplot**) | same |
| `diagnostics$random_effects$intercepts` | data frame of BLUPs | line 970 |
| `diagnostics$random_effects$slopes` | numeric vector (when `random_effects_specification = "slope"`) | line 972 |
| `diagnostics$variance_components` | `lme4::VarCorr()` output | line 976 |

**The data is computed and stored, but NO rendered ggplot is produced.**

**Dashboard impact (BUG):** the dashboard's `output$tg_diag_qq` and `output$tg_diag_resid_fitted` look for `rv$tg_results$diag_qq_plot` and `rv$tg_results$diag_resid_fitted_plot` — fields the backend never sets. Every LME4 TG run shows "Q-Q plot not available" and "Residuals vs Fitted not available" in the TG Diagnostics tab. See section 7.

### AUC path (`tgs_path_auc.R`)

| Field | Type |
|---|---|
| `qq_plot` | `qqnorm()` output (data, not a ggplot) |

Same issue — data computed, no rendered ggplot.

### GAM path (`tgs_gam.R`)

| Field | Type |
|---|---|
| `qq_plot` | `qqnorm()` output |
| Various `gam_diagnostics` fields | from `tgs_gam_diagnostics()` |

**What's missing for GAM specifically:**

- **`mgcv::k.check()`** — the basis-dimension check. Tells you whether the smoother had enough basis functions to capture the wiggliness in the data. Without this, GAM users can't tell if their `k = default` was too small. This is the #1 GAM-specific diagnostic and it's absent.
- **`mgcv::concurvity()`** — the GAM analog of multicollinearity. Detects when smooths of different covariates explain overlapping variance. Important when GAM is fit with multiple smooths.
- **Smooth-specific Q-Q** via `mgcv::qq.gam()` — better-calibrated than the generic residuals Q-Q.

---

## 5. Other frequentist analyses

### Body Weight (`analyze_body_weight`)

The only frequentist path with properly rendered diagnostic ggplots:

| Field | Type | Where set |
|---|---|---|
| `diag_qq_plot` | ggplot — Q-Q with reference line | `R/analyze_body_weight.R:299-307` |
| `diag_resid_fitted_plot` | ggplot — residuals vs fitted with loess | line 310-318 |

Rendered correctly in the Toxicity Diagnostics tab.

**Assessment:** correct implementation. The reference line on the QQ plot helps users see deviations from normality. The loess on residuals-vs-fitted catches non-linearity that a raw scatter wouldn't.

Added in response to CODE_REVIEW.md J.12 (the BW model didn't expose any frequentist diagnostics before — fixed).

### Synergy (`analyze_drug_synergy`)

Returns t-test p-values for combo-vs-monotherapy. **No model-assumption diagnostics.**

- **Welch's t-test** doesn't require equal variances, so a Levene check isn't strictly needed.
- **Normality of the AUC distribution** — would be checkable with a QQ; not produced.
- **Outlier identification** — boxplot of per-mouse AUC by group would surface candidates; not produced.

The synergy module is the lightest on diagnostics. For exploratory work this is acceptable; for reporting, the absence of a "does the t-test assumption hold" check is a gap.

### Dose-Response — frequentist (`dose_response_statistics`)

Returns the fitted `linear_model`, `anova_model`, `statistics$ec50`, etc. **No goodness-of-fit diagnostics.**

`drc` (the package used for Hill / Emax fitting) exposes:
- `drc::modelFit()` — lack-of-fit test against a smoother
- `drc::mselect()` — model comparison (4PL vs 3PL vs 2PL)
- `plot(drm_model, type = "residuals")` — residuals plot
- `drc::EDcomp()` — comparison of EC50s across groups

**None of these are surfaced.** For a regulatory dose-response analysis, the lack-of-fit test is a near-essential check. Currently the user only sees parameter estimates and ANOVA — they can't tell whether the Hill equation is the right functional form for their data.

### Repeated Measures ANOVA (`repeated_measures_anova`)

Returns only `model`, `anova_table`, `transform`. **No diagnostics returned.**

Mixed-effects models don't need sphericity (see [`METHODS.md`](METHODS.md)) but the standard LMM diagnostics still apply:
- Residual Q-Q
- Random-effects Q-Q (do the BLUPs look normal?)
- Residuals vs fitted

None are produced. Worth adding even a minimal `diagnostics` field for symmetry with `analyze_body_weight`.

---

## 6. Frequentist LME4 / LMM diagnostics — what's broadly missing

Across `tumor_growth_statistics()`, `analyze_body_weight()`, `repeated_measures_anova()`, the LMM fits share the same underlying model. Several useful LMM-specific diagnostics aren't computed anywhere:

| Diagnostic | What it shows | R implementation |
|---|---|---|
| **Random-effects Q-Q plot** | Whether the BLUPs (per-animal intercepts/slopes) look normally distributed. Required by lme4 model assumption | `lattice::qqmath(ranef(model))` or `ggplot2::stat_qq` on `lme4::ranef(model)` |
| **Scale-location plot** (residuals vs fitted, sqrt-abs residuals) | Heteroscedasticity check | `ggplot(aes(fitted, sqrt(abs(residuals))))` |
| **Cook's distance / leverage** | High-influence observations | `influence.ME::influence(model)`, then `cooks.distance(infl)` |
| **DFBETAS** | Influence on each fixed-effect coefficient | `influence.ME::dfbetas(infl)` |
| **Levene's test (per group)** | Variance equality across treatment levels (relevant when `random_effects_specification = "intercept_only"`) | `car::leveneTest()` on residuals |
| **VIF (variance inflation factor)** | Multicollinearity among fixed effects | `car::vif(model)` |

For preclinical datasets with one categorical covariate, VIF and influence diagnostics are less load-bearing. The random-effects Q-Q is the most-important missing one.

---

## 7. Bugs and gaps to fix

A prioritized list. (1)-(3) are user-visible bugs; (4)-(8) are gaps; (9)-(13) are nice-to-haves.

### (1) TG Diagnostics tab is broken

**Symptom:** the Tumor Growth → Diagnostics tab shows "Q-Q plot not available" and "Residuals vs Fitted not available" on every LME4 fit.

**Cause:** field-name mismatch. The dashboard reads:

```r
rv$tg_results$diag_qq_plot
rv$tg_results$diag_resid_fitted_plot
```

But `tumor_growth_statistics()` returns:

```r
diagnostics$residuals$qq_plot   # qqnorm output, not a ggplot
# (no diag_resid_fitted_plot field at all)
```

**Fix:** in `R/tumor_growth_statistics.R` around line 962, mirror the pattern from `R/analyze_body_weight.R:291-318`:

```r
diag_qq_plot <- NULL
diag_resid_fitted_plot <- NULL
if (include_diagnostics && requireNamespace("ggplot2", quietly = TRUE)) {
  tryCatch({
    resid_vec <- stats::residuals(model)
    fit_vec   <- stats::fitted(model)
    # Q-Q plot
    qq_data <- stats::qqnorm(resid_vec, plot.it = FALSE)
    diag_qq_plot <- ggplot2::ggplot(
      data.frame(theoretical = qq_data$x, sample = qq_data$y),
      ggplot2::aes(x = theoretical, y = sample)) +
      ggplot2::geom_point(alpha = 0.5) +
      ggplot2::geom_abline(slope = 1, intercept = 0,
                           linetype = "dashed", colour = "red") +
      ggplot2::theme_classic() +
      ggplot2::labs(title = "Q-Q plot of residuals",
                    x = "Theoretical quantiles", y = "Sample quantiles")
    # Residuals vs Fitted
    diag_resid_fitted_plot <- ggplot2::ggplot(
      data.frame(fitted = fit_vec, residuals = resid_vec),
      ggplot2::aes(x = fitted, y = residuals)) +
      ggplot2::geom_point(alpha = 0.5) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
      ggplot2::geom_smooth(method = "loess", se = FALSE, colour = "blue") +
      ggplot2::theme_classic() +
      ggplot2::labs(title = "Residuals vs Fitted",
                    x = "Fitted values", y = "Residuals")
  }, error = function(e) NULL)
}
# Add to the returned list
```

Then add both fields to the `return(list(...))` block at the end. The dashboard side already works.

**Severity:** medium (one whole UI tab is broken; no incorrect results, just no diagnostics).

### (2) Bayesian residuals / posterior-dist / credible-intervals plots wasted

**Symptom:** `bayesian_tumor_growth()`, `bayesian_body_weight()`, `bayesian_survival()`, etc. compute three plots that are never displayed:

- `residuals_plot` — Bayesian residuals vs day, faceted by group
- `posterior_dist_plot` — `bayesplot::mcmc_areas` for treatment parameters
- `credible_intervals_plot` — forest plot of group EMMs

The dashboard never renders these. They consume compile + render time when `plots = TRUE` (which is the default for BW, Survival, Synergy, DR, TWM) and the result is discarded.

**Fix options:**

a. **Render them.** Add to each module's MCMC Diagnostics tab or a new "Bayesian Residuals" sub-tab. Easy win — backend is already producing them.
b. **Stop computing them.** Pass `plots_extras = FALSE` (would need a new argument) so brms doesn't generate them. Saves a few seconds per Bayesian fit.

Recommend (a) for the first two (they're genuinely useful diagnostics); (b) for `credible_intervals_plot` since the dashboard already builds its own forest plot from the data.

**Severity:** low (no user-visible bug; just inefficient).

### (3) GAM / GAMM has no rendered diagnostics

**Symptom:** GAM model_type produces no actionable diagnostics. Users can't tell if `k` was big enough or if there's concurvity.

**Fix:** in `R/tgs_gam.R::tgs_gam_diagnostics()`, add:

```r
k_check <- tryCatch(mgcv::k.check(gam_obj$gam), error = function(e) NULL)
concurvity <- tryCatch(mgcv::concurvity(gam_obj$gam, full = FALSE),
                       error = function(e) NULL)
```

Surface in the GAM result list and add a dashboard tab "GAM Diagnostics" rendering both as DT tables with a banner when `k.index < 1`.

**Severity:** medium (GAM users get no quality signal at all).

### (4) Schoenfeld plot missing for Cox PH

**Symptom:** the PH Assumption tab shows the `cox.zph` test table but not the Schoenfeld residuals plot.

**Fix:** add a `plotOutput` rendering `plot(rv$surv_results$ph_test)`. Two lines of code in `mod_survival.R`.

**Severity:** medium (the test tells you "PH is violated for term X"; the plot tells you HOW — linear drift, periodic, sudden change — which determines the fix).

### (5) Frequentist DR has no goodness-of-fit

**Symptom:** the user fits a Hill / Emax curve and gets no indication of fit quality. Doesn't know whether to switch to a 3PL or 2PL.

**Fix:** in `R/dose_response_statistics.R`, add:

```r
fit_diagnostics <- list(
  lack_of_fit  = tryCatch(drc::modelFit(model), error = function(e) NULL),
  model_select = tryCatch(drc::mselect(model,
                                       fctList = list(drc::LL.3(), drc::LL.2())),
                          error = function(e) NULL),
  residuals    = stats::residuals(model)
)
```

Plus a residual plot. Surface in a DR Diagnostics tab.

**Severity:** medium-high (regulatory-quality DR analyses need lack-of-fit).

### (6) Frequentist synergy has no diagnostics

The t-test path returns p-values without checking the t-test's assumptions. At minimum, surface per-group AUC normality (QQ) and a boxplot of per-mouse AUCs.

**Severity:** low (Welch is robust; nice-to-have).

### (7) Random-effects QQ missing across all LMM fits

Across `tumor_growth_statistics()`, `analyze_body_weight()`, `repeated_measures_anova()`: the BLUPs distribution is never QQ-checked. This is an LMM model assumption (random effects are normally distributed); easy to compute via `lme4::ranef()`.

**Severity:** low-medium (assumption is rarely problematic in practice but should be checkable).

### (8) `repeated_measures_anova()` has no diagnostics at all

Returns just the model and ANOVA table. Should at least return residual QQ + residuals-vs-fitted.

**Severity:** low (the function is rarely used as a primary analysis).

### (9) ESS adequacy flag missing

The Bayesian `mcmc_diagnostics` has a `Converged` flag (Rhat ≤ 1.01) but no flag for low ESS. Add `ESS_Adequate = Bulk_ESS >= 400 & Tail_ESS >= 400`.

### (10) PPC coverage visualisation

Current dl row: "50%: 48% | 80%: 78% | 95%: 94% (n_obs = 120)". Easy to read but a small bar chart (bars at 50/80/95 with target line) makes under-coverage pop.

### (11) Pareto-k distribution histogram

The dashboard currently shows the table of high-k observations. A small histogram of all Pareto-k values with the 0.7 cut-off line gives context.

### (12) Bayesian rank plots

Modern Stan diagnostic that replaces Rhat for finer convergence detection. `bayesplot::mcmc_rank_overlay` would add a tab.

### (13) Influence diagnostics for LMM

`influence.ME::influence(model)` is the standard. Compute and surface Cook's distance + DFBETAS for high-leverage animals.

---

## 8. Recommended dashboard tab structure (post-fix)

After fixing (1)-(8), the diagnostic surface across modules should look like this:

### Tumor Growth — Diagnostics tab
- Q-Q plot (residuals; LME4/AUC/GAM)
- Residuals vs Fitted (LME4/AUC/GAM)
- Scale-location plot
- Random-effects Q-Q (when `random_effects_specification` is non-trivial)
- For GAM: k.check table + concurvity table

### Tumor Growth — Posterior Check tab (Bayesian)
- PPC density overlay (`pp_check_plot`)
- PPC coverage bar chart (new)

### Tumor Growth — MCMC Diagnostics tab (Bayesian)
- Sampler warnings banner
- `bayes_diagnostics_panel()` (Bayes R² + PPC coverage + LOO)
- MCMC diagnostics table (per-parameter Rhat / ESS / ESS_Adequate flag)
- Pareto-k histogram (new)
- Influential observations table
- MCMC trace plot
- Posterior density areas (new — render `posterior_dist_plot`)
- Bayesian residuals vs day (new — render `residuals_plot`)

### Tumor Growth — Prior vs Posterior tab (Bayesian) — unchanged

### Survival — PH Assumption tab
- cox.zph table (existing)
- Schoenfeld residuals plot (new)

### Survival — Diagnostics tab (new)
- Martingale residuals
- Deviance residuals
- C-index summary (move from current "model_quality_panel" location)

### Toxicity / Body Weight — Diagnostics tab (existing, works correctly)

### Toxicity / Body Weight — Bayesian sub-tabs — same additions as TG

### Dose-Response — Diagnostics tab (new)
- Lack-of-fit test (drc::modelFit)
- Model selection table (drc::mselect)
- Residuals plot

### Drug Synergy — Diagnostics tab (light, new)
- Per-group AUC Q-Q
- Per-group AUC boxplot

---

## 9. Implementation priorities

If you fix in one batch:

1. **Bug fix (1)** — TG Diagnostics tab. Two-helper-function addition. Half-day.
2. **(4) Schoenfeld plot** — three lines in `mod_survival.R`. One hour.
3. **(3) GAM k.check + concurvity** — backend changes in `tgs_gam.R`; one new dashboard tab. Half-day.
4. **(5) DR lack-of-fit** — backend changes in `dose_response_statistics.R`; new dashboard tab. Half-day.

That's a clean ~2-day batch closing the four highest-impact items.

If you want to also close (2) (render the wasted plots), add one more half-day.

The remaining items are quality-of-life rather than load-bearing; consider tackling them when one of the underlying modules is being reworked for another reason.

---

## 10. Reference

For background on the diagnostics themselves (what each tests, thresholds, interpretation):

- Bates et al. 2015 — "Fitting Linear Mixed-Effects Models Using lme4" (JSS 67(1)) — the lme4 diagnostic primer
- Therneau & Grambsch 2000 — *Modeling Survival Data: Extending the Cox Model* — Schoenfeld residuals + cox.zph
- Wood 2017 — *Generalized Additive Models: An Introduction with R* — mgcv diagnostics
- Vehtari et al. 2021 — "Rank-normalization, folding, and localization: An improved Rhat for assessing convergence of MCMC" (Bayesian Analysis 16(2)) — modern Rhat / ESS guidance
- Vehtari et al. 2017 — "Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC" (Statistics and Computing 27(5)) — PSIS-LOO + Pareto-k
- Gelman et al. 2019 — "R-squared for Bayesian regression models" (American Statistician 73(3)) — Bayes R²
- Ritz et al. 2015 — "Dose-Response Analysis Using R" (PLOS ONE 10(12)) — drc diagnostics

---

## See also

- [`METHODS.md`](METHODS.md) — when to use each method
- [`BAYESIAN.md`](BAYESIAN.md) — Bayesian diagnostics in depth
- `R/utils_bayes.R` — Bayesian diagnostic implementations
- `R/survival_statistics.R` — Cox PH + concordance implementation
- `R/analyze_body_weight.R` — the only frequentist path with rendered diagnostic ggplots (template for fixing TG)
- Dashboard `R/helpers_bayes.R::bayes_diagnostics_panel()` — unified Bayesian diagnostic renderer
