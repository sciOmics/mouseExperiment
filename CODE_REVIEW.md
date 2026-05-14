# mouseExperiment Code Review

**Reviewer:** Claude Code  
**Date:** 2026-05-14  
**Version reviewed:** 0.3.3  
**Scope:** All R source files under `R/`

---

## Summary

The package covers a wide, useful analysis surface and the recent refactors (composite-key fixes, Loewe replacement) improved correctness. However, there are several statistical validity problems, a handful of runtime bugs, and recurring code-quality issues that should be addressed before wide use in publications or regulatory contexts.

Issues are categorised **Critical** (wrong answers or crashes), **Major** (misleading output or silent wrong behaviour), and **Minor** (code quality, inconsistency, efficiency).

Issues marked **FIXED** have been resolved in the codebase.

---

## 1. Critical Issues

### 1.1 Silently ignored `handle_cage_effects` and `random_effects_specification` parameters — FIXED
**File:** `R/tumor_growth_statistics.R` — `tgs_fit_lme4_models()`

Both parameters are validated with `match.arg()` and documented, but are **never acted upon**. The model formula is always `Volume ~ Day * Treatment + (1|ID)` or `Volume ~ Day * Treatment + (Day|ID)`, selected by BIC — the user's `handle_cage_effects` and `random_effects_specification` arguments have zero effect on the model. Callers specifying `handle_cage_effects = "always_include"` or `random_effects_specification = "none"` will silently receive a model that ignores their choice.

The same applies to `polynomial_degree` — the parameter is accepted but never inserted into the formula.

**Resolution:** `tgs_fit_lme4_models()` now accepts `random_effects_specification`, `handle_cage_effects`, `polynomial_degree`, and `cage_collinear` as explicit parameters and builds the model formula from them. `random_effects_specification = "none"` uses `lm()` instead of `lmer()`. `polynomial_degree > 1` inserts `I(Day^n)` terms. `handle_cage_effects` adds cage to fixed or random effects according to the chosen option. The BIC auto-selection between two fixed model structures has been removed. The diagnostics section guards `lme4::ranef()` / `VarCorr()` calls with `inherits(model, "lmerMod")` to handle the `"none"` path.

---

### 1.2 Runtime error: `verbose` not in scope in helper functions — FIXED
**Files:** `R/dose_response_statistics.R` — `try_nonlinear_models()` and `analyze_polynomial_trends()`

Both functions reference `verbose` (lines ~341 and ~452) but it is not in their parameter lists. When the main function is called with `verbose = TRUE` (the default), execution will reach these helpers and fail with `Error: object 'verbose' not found`.

```r
# try_nonlinear_models — no verbose parameter, but uses it:
if (verbose) {
  message("Selected dose-response model type: ", model_type)
```

**Resolution:** Added `verbose = TRUE` to both function signatures and passed it through from `perform_statistical_analyses()`.

---

### 1.3 Runtime error: `verbose` not in scope in `fit_survival_model()` — FIXED
**File:** `R/survival_statistics.R` — `fit_survival_model()` logrank branch (~line 454)

```r
if (verbose) message(paste(utils::capture.output(print(surv_diff)), collapse = "\n"))
```

`verbose` is a parameter of the exported `survival_statistics()`, not of the internal `fit_survival_model()`. This crashes when `has_issues = TRUE` and `firth_correction = FALSE`.

**Resolution:** Added `verbose = TRUE` to `fit_survival_model()` signature and passed it from `survival_statistics()`.

---

### 1.4 `method_used` extracted from wrong object type in survival model — FIXED
**File:** `R/survival_statistics.R` — `fit_survival_model()` coxphf branch (~lines 437–445)

After the Firth model block, the code does:
```r
method_used <- results$method_used
```
but `results` at that point is the data.frame returned by the inner `tryCatch`, not a list with a `method_used` field. `results$method_used` will be `NULL`, and the function returns `method_used = NULL` instead of `"coxphf"`.

**Resolution:** Reordered the three assignments so `method_used <- results$method_used` is extracted from the tryCatch list _before_ `results` is overwritten with `results$results`.

---

### 1.5 `tumor_auc_analysis()` composite ID uses `_` as separator — fragile parsing — FIXED
**File:** `R/tumor_auc_analysis.R` (~lines 110, 229)

```r
composite_ids <- paste(df[[id_column]], df[[treatment_column]], df[[cage_column]], sep = "_")
...
id_parts <- strsplit(unique_id, "_")[[1]]
original_id <- id_parts[1]
treatment   <- id_parts[2]   # wrong if treatment name contains "_"
cage        <- id_parts[3]
```

If any ID, treatment name, or cage name contains `_` (e.g. `"Drug_A"`, `"Group_1"`), the split produces more than three elements and `id_parts[2]` / `id_parts[3]` extract the wrong fields silently. Every other composite-ID constructor in the package uses `"|||"` or `"⁠"`, which are safe delimiters.

**Resolution:** Replaced with `make_mouse_key()` / `split_mouse_key()` helpers (see Issue 3.1). All composite-ID construction and parsing across the package now uses the `"|||"` delimiter via these shared utilities.

---

### 1.6 `efficacy_toxicity_bivariate()` — control AUC computed over all mice pooled — FIXED
**File:** `R/efficacy_toxicity_bivariate.R` (~lines 104–108)

```r
ctrl_auc <- trap_auc(ctrl_data$Day, ctrl_data$Volume)
```

`ctrl_data` contains all control observations from all mice. `trap_auc` sorts by time and applies the trapezoidal rule, so every (Day, Volume) pair across all control mice is treated as a single time series. For any day with multiple mice, this interleaves observations and produces a nonsensical integral. The per-mouse AUC for each control subject should be computed first and averaged.

**Resolution:** Compute per-mouse control AUCs using composite keys, then take their mean (`ctrl_mean_auc`). The per-mouse efficacy loop now references `ctrl_mean_auc` instead of calling `trap_auc` on pooled control data.

---

### 1.7 `total_benefit_area()` — baseline computed from first row of group aggregate — FIXED
**File:** `R/total_benefit_area.R` (~line 87)

```r
baseline_grp <- stats::aggregate(Net_Weight ~ Treatment, data = wd, FUN = function(x) x[1])
```

`aggregate()` does not guarantee that `x[1]` is the earliest observation; the row order within groups depends on the data's initial order. The result is that the toxicity AUC is computed relative to an arbitrary weight observation rather than the true baseline. The same pattern appears in `therapeutic_window_metric.R:74` and `weight_corrected_tgi.R:63`.

**Resolution:** All three files now filter to `wd[wd$Day == min(wd$Day), ]` before aggregating, ensuring the baseline is always the first study day regardless of row order.

---

### 1.8 `weight_corrected_tgi()` — composite key not used for ID uniqueness — FIXED
**File:** `R/weight_corrected_tgi.R` (~lines 63, 85)

The function aggregates by `ID` alone:
```r
baseline <- stats::aggregate(Net_Weight ~ ID, data = wd, FUN = function(x) x[1])
```
and then filters out excluded mice by `ID`:
```r
safe_data <- wd[!wd$ID %in% excluded_ids, ]
```

If mouse IDs are shared across treatment groups (the exact problem fixed in `body_weight_auc.R` v0.3.1), all mice sharing an ID in any group will be excluded or have wrong baselines.

**Resolution:** Baseline aggregation now uses `ID + Treatment`. Exclusion uses `make_mouse_key(Treatment, ID)` composite keys so an over-threshold mouse in one group does not incorrectly exclude a same-ID mouse in another group.

---

## 2. Major Statistical Issues

### 2.1 Post-hoc power analysis is statistically invalid as used — FIXED
**File:** `R/post_power_analysis.R` *(deleted)*

Post-hoc ("observed") power analysis is widely criticised and considered uninformative. Computing power from the observed effect size of the same study creates a 1-to-1 mapping: a non-significant result will always produce low post-hoc power; a significant result will produce high post-hoc power. The function provides no information beyond the p-value itself (Hoenig & Heisey, 2001; Gelman & Carlin, 2014).

**Resolution:** `post_power_analysis()` and the entire post-hoc mode were eliminated. `R/post_power_analysis.R` was deleted, its `NAMESPACE` export removed, and the dashboard `power_modules.R` rewritten to expose only the two a priori modes (`apriori_analytic` and `apriori_lmm`). All sample-size and power questions must now be answered prospectively with a user-supplied, design-based effect size.

---

### 2.2 Power simulation on linear scale, analysis on log scale
**File:** `R/apriori_power_simulation.R` (~lines 105–118)

The simulation generates data as:
```r
vol <- (baseline_volume + b0[j]) + (growth_rates[g] + b1[j]) * timepoints + noise
```
This is a **linear** data-generating model. However, `tumor_growth_statistics()` applies a **log transformation** before fitting. Simulating on the linear scale and fitting on the log scale produces anti-conservative power estimates — the residual variance after log-transformation is much smaller than on the raw scale, so the simulation understates noise, inflating apparent power.

**Fix:** Simulate log-volume data (i.e., the linear model on the log scale) and exponentiate, or simulate volumes as log-normal: `vol <- exp(log(baseline_volume) + b0[j] + (growth_rates[g] + b1[j]) * timepoints + noise)`, then fit with log-transformation as in the real analysis pipeline.

---

### 2.3 Cohen's d to f conversion for ANOVA power is scenario-specific
**File:** `R/apriori_power_analysis.R` (~line 90)

```r
f_val <- effect_size / sqrt(2)
```

This comment says "equal group means, two extreme groups." This conversion of Cohen's d to Cohen's f is valid only for the specific scenario where two groups have different means and all others are zero. For a typical oncology study where all treatment groups have reduced growth relative to control, the true Cohen's f would be different (generally smaller). Users specifying a Cohen's d for a pairwise comparison and assuming it drives the ANOVA correctly will over-estimate power.

**Fix:** Accept Cohen's f directly for ANOVA, or document the specific scenario assumed and direct users to compute f from their anticipated group means.

---

### 2.4 Bliss Independence applied to TGI — known limitation not documented
**Files:** `R/analyze_drug_synergy.R`, `R/analyze_drug_synergy_over_time.R`

Bliss Independence was formulated for probability of cell death, not for proportional inhibition (TGI). Applying it to TGI is a common pragmatic choice but carries specific implications: when individual drug effects are large (TGI > 50% each), the Bliss expected combined effect can approach 100%, making it nearly impossible to demonstrate synergy by the Bliss criterion regardless of actual biological interaction. There is no documentation of this limitation.

Similarly, the Loewe Additivity single-dose approximation (`min(FE_A + FE_B, 1) / FE_combo`) is explicitly noted as requiring the "linear dose-response assumption," but the associated limitations (unknown dose-response curvature, absence of IC50 estimates) are not communicated to the user.

**Fix:** Add a `Notes` section to each function's documentation explaining the assumptions and appropriate interpretation range.

---

### 2.5 `tumor_auc_analysis()` LOCF method is not an AUC
**File:** `R/tumor_auc_analysis.R` (~lines 202–217)

```r
} else if (method == "last_observation") {
  ...
  return(list(auc = last_volume + extrapolated_value, extrapolated = is_extrapolated))
}
```

The "last observation" return value is `last_volume` (a single volume measurement in mm³), not an area. Adding `extrapolated_value` (an area in mm³·day) to `last_volume` (mm³) is dimensionally incoherent. The result will be numerically close to the last volume for short extrapolation windows but structurally wrong. LOCF should carry the last observation forward to fill the time-series and then apply the trapezoidal rule over the full period.

---

### 2.6 `dose_response_statistics()` — AIC comparison between `lm` and `drc` models is invalid
**File:** `R/dose_response_statistics.R` (~lines 359–362)

```r
statistics$linear_aic  <- AIC(linear_model)
statistics$nonlinear_aic <- AIC(dr_model)
```

AIC values from `stats::lm` and `drc::drm` are not directly comparable because `drc` uses a different parameterisation of the likelihood (it does not include the `log(2π)` constant term in the same way and estimates separate variance parameters). Comparing these values as if they were on the same scale will mislead model selection.

---

### 2.7 `tumor_growth_statistics()` — emmeans contrast at mean time point only
**File:** `R/tumor_growth_statistics.R` (~lines 955–1003)

```r
lsmeans_obj <- emmeans::emmeans(model, specs = treatment_column,
  at = stats::setNames(list(mean(analysis_df[[time_column]])), time_column))
```

Marginalising over a single time point (the mean day) to summarise treatment effects from a `Day * Treatment` interaction model discards the full interaction structure. If treatments diverge over time (the biologically relevant signal), the marginal mean at a single time point is an incomplete summary. The user should receive estimated marginal means across multiple time points or the interaction terms themselves should be the primary output.

---

### 2.8 Survival analysis: logrank p-value inflated for multi-group comparison
**File:** `R/survival_statistics.R` — `fit_survival_model()` (~line 457–458)

When multiple treatment groups have separation and the logrank fallback is used:
```r
p_value <- 1 - stats::pchisq(chisq, df = length(treatment_groups) - 1)
```
This single omnibus p-value is assigned to all non-reference groups:
```r
results$P_Value[-ref_idx] <- p_value
```

An omnibus test p-value is not a valid per-comparison p-value. Users may interpret these identical p-values as pairwise comparisons with the reference group when they are not.

---

### 2.9 `therapeutic_window_metric()` — max weight loss taken as maximum across group, not mean
**File:** `R/therapeutic_window_metric.R` (~lines 87–89)

```r
group_wl <- stats::aggregate(Pct_Loss ~ Treatment, data = mouse_wl, FUN = max, na.rm = TRUE)
```

The TWM denominator is the single worst-case mouse in each group. This makes TWM hypersensitive to a single outlier and not representative of typical toxicity. A mean or median weight loss would be more interpretable.

---

## 3. Code Quality Issues

### 3.1 Composite ID delimiter inconsistency
**Files:** Multiple

| File | Separator |
|------|-----------|
| `tumor_growth_statistics.R` | `"|||"` |
| `tumor_auc_analysis.R` | `"_"` |
| `body_weight_auc.R` | `"⁠"` (word-joiner) |
| `weight_loss_threshold.R` | `"⁠"` |
| `efficacy_toxicity_bivariate.R` | `"⁠"` |
| `post_power_analysis.R` | `"|||"` |

Using `_` (common in IDs) is actively dangerous. A single shared helper function (e.g., `make_mouse_key()` / `split_mouse_key()`) using `"|||"` should replace all variants.

---

### 3.2 `rbind` growing a data frame inside a loop — O(n²) memory copies — FIXED
**Files:** `R/analyze_drug_synergy_over_time.R` (~line 230), `R/post_power_analysis.R` (many loops)

```r
# analyze_drug_synergy_over_time.R
synergy_summary <- rbind(synergy_summary, tp_row)

# post_power_analysis.R (ppa_calculate_effect_sizes, ppa_estimate_sample_sizes, etc.)
effect_sizes_df <- rbind(effect_sizes_df, data.frame(...))
```

R copies the entire data frame on every `rbind` call. For `post_power_analysis.R`, the pattern is nested across four levels of loops. Pre-allocate a list and call `do.call(rbind, list)` once after the loop.

**Resolution:** `analyze_drug_synergy_over_time.R` now accumulates rows in `synergy_rows` list and combines with `do.call(rbind, ...)` after the loop. In `post_power_analysis.R`, all inner loops in `ppa_calculate_effect_sizes()`, `ppa_run_power_analysis()`, `ppa_estimate_sample_sizes()`, and `create_all_pairwise_effects()` were converted to pre-allocated list accumulation.

---

### 3.3 `calculate_volume()` — `in_place` modifies parent environment via `assign()` — FIXED
**File:** `R/calculate_volume.R` (~lines 126–131)

```r
if (in_place) {
  parent_frame <- parent.frame()
  df_name <- deparse(substitute(df))
  assign(df_name, result, envir = parent_frame)
}
```

`parent.frame()` refers to the immediate caller's environment. This fails silently when `calculate_volume` is called from inside another function (the assignment modifies the inner function's frame, not the user's workspace). The deprecation warning is appropriate; the implementation should simply stop trying to side-effect the parent and just return the result.

**Resolution:** Removed the `parent.frame()` / `deparse(substitute())` / `assign()` block. The function now simply returns `result` in all cases.

---

### 3.4 `generate_summary_statistics()` always prints regardless of `verbose` — FIXED
**File:** `R/dose_response_statistics.R` (~lines 162–164)

```r
message("Summary statistics by dose level:")
message(paste(utils::capture.output(print(summary_stats)), collapse = "\n"))
```

The `verbose` parameter from the top-level call is not passed to this helper, so it always prints. This also applies to `generate_user_report()`, which is correctly gated with `if (verbose)` at the call site but internally has no guard.

**Resolution:** Added `verbose = TRUE` parameter to `generate_summary_statistics()` and gated its `message()` calls with `if (isTRUE(verbose))`. The call site in `dose_response_statistics()` now passes `verbose = verbose`.

---

### 3.5 `post_power_analysis.R` — excessive complexity with duplicated fallback logic — FIXED
**File:** `R/post_power_analysis.R`

`ppa_calculate_effect_sizes()` is ~350 lines with the same "use default effects if X fails" pattern repeated four times. The logic for extracting treatments is attempted in at least six different places in the same function. This depth of defensive coding masks bugs rather than preventing them.

**Recommendation:** Simplify to: validate inputs at entry, compute treatments once from data, compute effect sizes once from data or accept user-supplied values. Remove the cascading fallback logic.

**Resolution:** The three inline "fallback to default effects" branches (when `nrow(auc_summary) < 2`, when pooled SD is invalid, and when `AUC.SD` is absent) all replicated the same loop as `create_all_pairwise_effects()`. Each was replaced with a single `create_all_pairwise_effects(treatments, default_effect_sizes)` call. The parametric/simulation fallback branch was similarly collapsed.

---

### 3.6 Duplicated trapezoidal AUC implementations — FIXED
**Files:** `R/utils_auc.R`, `R/body_weight_auc.R`, `R/efficacy_toxicity_bivariate.R`, `R/total_benefit_area.R`

The body of `calculate_auc()` from `utils_auc.R` is copy-pasted as `trap_auc()` in at least three other files. These local copies are slightly different (e.g., they don't validate input lengths), leading to inconsistent edge-case handling. All callers should use the exported `calculate_auc()`.

**Resolution:** Removed the local `trap_auc()` definition from `body_weight_auc.R`, `efficacy_toxicity_bivariate.R`, and `total_benefit_area.R`. All call sites now use `calculate_auc()` directly.

---

### 3.7 `apriori_power_simulation.R` — p-value extraction logic is inverted — FIXED
**File:** `R/apriori_power_simulation.R` (~lines 143–173)

```r
if (nrow(int_rows) == 0) {
  # Fall back: Satterthwaite via lmerTest if available
  ...
} else {
  # Use likelihood-ratio test between full and reduced model
  ...
}
```

`lme4::lmer` does not compute p-values by default, so `int_rows` (filtered by `":"` in row names) will almost always be empty, because the coefficient table has no `Pr(>|t|)` column. The LRT branch (statistically preferred for LMMs) is reached only in the rare case where coefficients have p-values. The condition should be inverted: prefer LRT, fall back to Satterthwaite if LRT fails.

**Resolution:** Inverted the condition. The LRT (fit vs. fit_red) is now attempted first. Satterthwaite via `lmerTest` is the fallback when `fit_red` cannot be fitted.

---

### 3.8 `plot_synergy_trend()` — `size` aesthetic deprecation warning — FIXED
**File:** `R/analyze_drug_synergy_over_time.R` (~line 346)

```r
ggplot2::geom_line(ggplot2::aes(y = TGI_Combo, color = combo_name), size = 1.2)
```

`size` is deprecated in `geom_line()` since ggplot2 3.4.0; the replacement is `linewidth`. This generates a deprecation warning in every call.

---

### 3.9 Annotation coordinates in `plot_synergy_trend()` break when time starts at 0 — FIXED
**File:** `R/analyze_drug_synergy_over_time.R` (~lines 363–382)

```r
ggplot2::annotate("rect",
  xmin = min(synergy_summary$Time_Point) * 1.05,
  xmax = min(synergy_summary$Time_Point) * 1.15, ...)
```

When `Time_Point` starts at 0 or 1, these expressions evaluate to `0` and `0` (or tiny values near the y-axis), placing annotation boxes at or outside the visible plot area.

**Resolution:** Replaced multiplicative offsets with additive range-based offsets. Pre-computed `.tp_lo`, `.tp_span` (clamped to at least 1), `.tp_x1`, `.tp_x2`, `.tp_x3` before the `ggplot` call; used these in all four `annotate()` calls.

---

### 3.10 `tumor_auc_analysis()` checks `exists("plot_auc", mode = "function")` at runtime — FIXED
**File:** `R/tumor_auc_analysis.R` (~line 339)

```r
if (exists("plot_auc", mode = "function")) {
```

`plot_auc` is an exported function in the same package. It will always exist when `mouseExperiment` is loaded. The check is unnecessary and the fallback `ggplot2` code is dead code.

**Resolution:** Removed the `exists("plot_auc", mode = "function")` guard and the dead fallback branch. The call is now a single `tryCatch({ plot_auc(...) }, error = ...)` with a basic fallback only for genuine errors.

---

### 3.11 `survival_statistics.R` — median survival calculated redundantly in `print_results()` — FIXED
**File:** `R/survival_statistics.R` (~lines 578–697)

The `print_results()` function re-fits `survfit` and re-calculates median survival for groups where median is `NA` — work already done in `survival_statistics()`. The same logic appears twice, diverging in subtle ways. The calculated value is updated in the local copy of `results` inside `print_results()` but not propagated back to the caller, so the returned `result_list$results` may still contain `NA` medians. Move all median calculation into `survival_statistics()` before calling `print_results()`.

**Resolution:** Simplified both the loop body and the summary table section of `print_results()` to display-only: `NA` medians print as "Not reached"; non-NA medians print as formatted days. All median calculation remains in `survival_statistics()` where it is also written back to `result_list$results`.

---

### 3.12 Inconsistent `isTRUE(verbose)` vs `if (verbose)` usage — FIXED
Some functions use `if (isTRUE(verbose))` (safer, handles `NULL`/`NA`) and others use `if (verbose)`. Standardise on `isTRUE(verbose)` throughout.

**Resolution:** Replaced all `if (verbose)` with `if (isTRUE(verbose))` in `tumor_growth_statistics.R`, `dose_response_statistics.R`, and `survival_statistics.R`.

---

### 3.13 `dose_response_statistics.R` — `verbose` not declared in `dose_response_statistics()` signature — FIXED
**File:** `R/dose_response_statistics.R` (~line 52 vs usage)

`verbose = TRUE` appears in the function signature but is not in the `@param` documentation block, and is passed to `perform_statistical_analyses()` but not to `generate_summary_statistics()` or `create_dose_plots()`. The documentation-signature mismatch is a maintenance hazard.

**Resolution:** Added `@param verbose` to the roxygen documentation block. `verbose` is now also passed to `generate_summary_statistics()` (see 3.4).

---

### 3.14 `analyze_drug_synergy_over_time.R` redundant Bliss recalculation — FIXED
**File:** `R/analyze_drug_synergy_over_time.R` (~lines 182–208)

After calling `analyze_drug_synergy()` which already computed Bliss values, the function recomputes them identically and then validates that the two match. Since `analyze_drug_synergy()` is a pure function called with identical inputs, a mismatch is impossible. The validation block adds ~25 lines and a `Validation_Check` column in the output for zero benefit.

**Resolution:** Removed the fe_a/fe_b recalculation and validation block. Bliss expected value is now read directly from `summary_df$TGI_Percent[summary_df$Treatment == "Bliss Expected"]`. The `Validation_Check` column is removed from the output data frame and verbose printing.

---

## 4. API Design Concerns

### 4.1 `dose_response_statistics()` exports too many internal helpers — FIXED
**File:** `R/dose_response_statistics.R`

`prepare_dose_data()`, `generate_summary_statistics()`, and `create_dose_plots()` are exported (they appear in `NAMESPACE`) but are really internal steps of `dose_response_statistics()`. They have no independent use case, and exporting them exposes implementation details that are hard to maintain with a stable API. Mark them `@noRd` / `@keywords internal`.

**Resolution:** Changed `@export` to `@noRd` + `@keywords internal` on `prepare_dose_data()` and `generate_summary_statistics()` (the two that were incorrectly exported; `create_dose_plots()` was already internal). Both lines removed from `NAMESPACE`.

---

### 4.2 `post_power_analysis.R` — private helpers mixed into the same file without `@noRd` — FIXED
**File:** `R/post_power_analysis.R`

`calculate_auc_values()`, `create_all_pairwise_effects()`, and `extract_test_data_treatments()` are defined at the bottom of the file without any roxygen documentation or `@noRd` tags. They are not exported but are not clearly internal either. Add `@noRd` and `@keywords internal`.

**Resolution:** Added `#' @noRd` + `#' @keywords internal` roxygen blocks before all three functions.

---

### 4.3 `tumor_growth_statistics()` returns inconsistent structures for `lme4` vs `auc` modes — FIXED
The `lme4` path returns `pairwise_comparisons` (an `emmeans` contrast object) while the `auc` path returns `posthoc` (a list with a data frame inside). Code downstream (including the dashboard) must branch on `model_type` to access the same logical information. A consistent return structure would reduce coupling.

**Resolution:** Both paths now return both fields. The AUC return gains `pairwise_comparisons = posthoc$pairwise` (the data frame). The lme4 return gains `posthoc = list(method = ..., pairwise = as.data.frame(summary(pairwise_comp)))` so callers can use `$posthoc$pairwise` regardless of model type. Existing code using either field continues to work without modification.

---

### 4.4 `analyze_drug_synergy()` has a hardcoded special-case comment with real names — FIXED
**File:** `R/analyze_drug_synergy.R` (~lines 198–200)

```r
# Clean up combo_name if needed - handle "HDACi + PD1" vs "aPD1"
if (combo_name == "HDACi + PD1" && drug_b_name == "aPD1") {
  # Variables are already correct - no change needed
}
```

This is dead code (the body does nothing) referencing specific drug names from what appears to be a real experiment. The comment and block should be removed.

**Resolution:** Removed the entire dead if-block and its comments.

---

## 5. Quick-Reference Checklist

| # | Issue | Severity | File | Status |
|---|-------|----------|------|--------|
| 1.1 | `handle_cage_effects`, `random_effects_specification`, `polynomial_degree` silently ignored | Critical | `tumor_growth_statistics.R` | ✅ Fixed |
| 1.2 | `verbose` out-of-scope in `try_nonlinear_models` / `analyze_polynomial_trends` | Critical | `dose_response_statistics.R` | ✅ Fixed |
| 1.3 | `verbose` out-of-scope in `fit_survival_model` logrank branch | Critical | `survival_statistics.R` | ✅ Fixed |
| 1.4 | `method_used` extracted from wrong object in Firth branch | Critical | `survival_statistics.R` | ✅ Fixed |
| 1.5 | Composite ID uses `_` separator — breaks on IDs/treatments containing `_` | Critical | `tumor_auc_analysis.R` | ✅ Fixed |
| 1.6 | Control AUC pooled across all mice instead of per-mouse average | Critical | `efficacy_toxicity_bivariate.R` | ✅ Fixed |
| 1.7 | Baseline weight from `aggregate(x[1])` without ordering guarantee | Critical | `total_benefit_area.R`, `therapeutic_window_metric.R` | ✅ Fixed |
| 1.8 | Composite key not used — shared IDs across groups produce wrong TGI | Critical | `weight_corrected_tgi.R` | ✅ Fixed |
| 2.1 | Post-hoc power is statistically invalid when effect sizes are data-derived | Major | `post_power_analysis.R` | ✅ Fixed (deleted) |
| 2.2 | Simulation generates linear-scale data; analysis fits log-scale — anti-conservative power | Major | `apriori_power_simulation.R` | Open |
| 2.3 | d→f conversion for ANOVA is scenario-specific, not general | Major | `apriori_power_analysis.R` | Open |
| 2.4 | Bliss/Loewe assumptions on TGI not documented | Major | `analyze_drug_synergy*.R` | Open |
| 2.5 | LOCF AUC returns volume (mm³), not area (mm³·day) | Major | `tumor_auc_analysis.R` | Open |
| 2.6 | AIC comparison between `lm` and `drc` models is invalid | Major | `dose_response_statistics.R` | Open |
| 2.7 | emmeans at single mean time point discards interaction | Major | `tumor_growth_statistics.R` | Open |
| 2.8 | Omnibus logrank p-value assigned as per-group p-value | Major | `survival_statistics.R` | Open |
| 2.9 | TWM denominator is max of one mouse, not group mean | Major | `therapeutic_window_metric.R` | Open |
| 3.1 | Composite ID separator inconsistent across files | Minor | Multiple | ✅ Fixed |
| 3.2 | `rbind` in loops — O(n²) | Minor | Multiple | ✅ Fixed |
| 3.3 | `in_place` assigns to parent frame via `assign()` | Minor | `calculate_volume.R` | ✅ Fixed |
| 3.4 | `generate_summary_statistics()` ignores `verbose` | Minor | `dose_response_statistics.R` | ✅ Fixed |
| 3.5 | `ppa_calculate_effect_sizes()` — excessive complexity, duplicated fallback | Minor | `post_power_analysis.R` | ✅ Fixed |
| 3.6 | Trapezoidal AUC copy-pasted in 3 files instead of using `calculate_auc()` | Minor | Multiple | ✅ Fixed |
| 3.7 | p-value extraction condition inverted — LRT branch rarely reached | Minor | `apriori_power_simulation.R` | ✅ Fixed |
| 3.8 | `size` deprecated in `geom_line` | Minor | `analyze_drug_synergy_over_time.R` | ✅ Fixed |
| 3.9 | Synergy plot annotation coordinates break at `Time_Point = 0` | Minor | `analyze_drug_synergy_over_time.R` | ✅ Fixed |
| 3.10 | Dead-code check for `plot_auc` existence | Minor | `tumor_auc_analysis.R` | ✅ Fixed |
| 3.11 | Median survival recalculated inside `print_results()` — result not propagated | Minor | `survival_statistics.R` | ✅ Fixed |
| 3.12 | Inconsistent `isTRUE(verbose)` vs `if (verbose)` | Minor | Multiple | ✅ Fixed |
| 3.13 | `verbose` parameter undocumented in `dose_response_statistics()` | Minor | `dose_response_statistics.R` | ✅ Fixed |
| 3.14 | Redundant Bliss recalculation in `analyze_drug_synergy_over_time()` | Minor | `analyze_drug_synergy_over_time.R` | ✅ Fixed |
| 4.1 | Internal helper functions incorrectly exported | Minor | `dose_response_statistics.R` | ✅ Fixed |
| 4.2 | Private helpers in `post_power_analysis.R` missing `@noRd` | Minor | `post_power_analysis.R` | ✅ Fixed |
| 4.3 | Inconsistent return structure between `lme4` and `auc` model types | Minor | `tumor_growth_statistics.R` | ✅ Fixed |
| 4.4 | Dead code referencing specific drug names | Minor | `analyze_drug_synergy.R` | ✅ Fixed |
