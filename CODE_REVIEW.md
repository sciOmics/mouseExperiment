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

### 2.2 Power simulation on linear scale, analysis on log scale — FIXED
**File:** `R/apriori_power_simulation.R` (~lines 105–118)

The simulation generates data as:
```r
vol <- (baseline_volume + b0[j]) + (growth_rates[g] + b1[j]) * timepoints + noise
```
This is a **linear** data-generating model. However, `tumor_growth_statistics()` applies a **log transformation** before fitting. Simulating on the linear scale and fitting on the log scale produces anti-conservative power estimates — the residual variance after log-transformation is much smaller than on the raw scale, so the simulation understates noise, inflating apparent power.

**Resolution:** Data are now generated on the log scale:
```r
log_vol <- log(baseline_volume) + b0[j] + (growth_rates[g] + b1[j]) * timepoints + noise
vol <- exp(log_vol)
```
The LMM is fitted as `log(Volume) ~ Treatment * Day + (Day | ID)` (full and reduced models), matching `tumor_growth_statistics()`. Function-argument defaults updated to log-scale values (`control_growth_rate = 0.15`, `treatment_effect = 0.10`, `random_intercept_sd = 0.20`, `random_slope_sd = 0.05`, `residual_sd = 0.10`). Dashboard LMM defaults and labels updated to match.

---

### 2.3 Cohen's d to f conversion for ANOVA power is scenario-specific — FIXED
**File:** `R/apriori_power_analysis.R` (~line 90)

```r
f_val <- effect_size / sqrt(2)
```

This comment says "equal group means, two extreme groups." This conversion of Cohen's d to Cohen's f is valid only for the specific scenario where two groups have different means and all others are zero. For a typical oncology study where all treatment groups have reduced growth relative to control, the true Cohen's f would be different (generally smaller). Users specifying a Cohen's d for a pairwise comparison and assuming it drives the ANOVA correctly will over-estimate power.

**Resolution:** The assumption is now fully documented in the `@param effect_size` roxygen block, in an inline comment at the conversion site, and in a `method_note` field returned when `n_groups >= 3`. Users who have a better estimate of between-group variability are directed to supply Cohen's f directly.

---

### 2.4 Bliss Independence applied to TGI — known limitation not documented — FIXED
**Files:** `R/analyze_drug_synergy.R`, `R/analyze_drug_synergy_over_time.R`

Bliss Independence was formulated for probability of cell death, not for proportional inhibition (TGI). Applying it to TGI is a common pragmatic choice but carries specific implications: when individual drug effects are large (TGI > 50% each), the Bliss expected combined effect can approach 100%, making it nearly impossible to demonstrate synergy by the Bliss criterion regardless of actual biological interaction. There is no documentation of this limitation.

Similarly, the Loewe Additivity single-dose approximation (`min(FE_A + FE_B, 1) / FE_combo`) is explicitly noted as requiring the "linear dose-response assumption," but the associated limitations (unknown dose-response curvature, absence of IC50 estimates) are not communicated to the user.

**Resolution:** Added `@section Assumptions and Limitations:` to both functions' roxygen documentation. The Bliss ceiling effect (individual TGI > 50%) and the Loewe single-dose approximation limitation are now clearly described.

---

### 2.5 `tumor_auc_analysis()` LOCF method is not an AUC — FIXED
**File:** `R/tumor_auc_analysis.R` (~lines 202–217)

```r
} else if (method == "last_observation") {
  ...
  return(list(auc = last_volume + extrapolated_value, extrapolated = is_extrapolated))
}
```

The "last observation" return value is `last_volume` (a single volume measurement in mm³), not an area. Adding `extrapolated_value` (an area in mm³·day) to `last_volume` (mm³) is dimensionally incoherent. The result will be numerically close to the last volume for short extrapolation windows but structurally wrong. LOCF should carry the last observation forward to fill the time-series and then apply the trapezoidal rule over the full period.

**Resolution:** The LOCF branch now computes the trapezoidal AUC for the observed period first via `calculate_auc(times, volumes)`, then adds the LOCF extension `(max_experiment_time − subject_max_time) × last_volume` as an mm³·day area. The dimensionally incoherent `last_volume + extrapolated_value` is gone.

---

### 2.6 `dose_response_statistics()` — AIC comparison between `lm` and `drc` models is invalid — FIXED
**File:** `R/dose_response_statistics.R` (~lines 359–362)

```r
statistics$linear_aic  <- AIC(linear_model)
statistics$nonlinear_aic <- AIC(dr_model)
```

AIC values from `stats::lm` and `drc::drm` are not directly comparable because `drc` uses a different parameterisation of the likelihood (it does not include the `log(2π)` constant term in the same way and estimates separate variance parameters). Comparing these values as if they were on the same scale will mislead model selection.

**Resolution:** Added a prominent `WARNING` comment at the AIC storage site explaining the incompatibility. A `statistics$aic_comparison_note` field is added to the return value with the same warning for users who access results programmatically.

---

### 2.7 `tumor_growth_statistics()` — emmeans contrast at mean time point only — FIXED
**File:** `R/tumor_growth_statistics.R` (~lines 955–1003)

```r
lsmeans_obj <- emmeans::emmeans(model, specs = treatment_column,
  at = stats::setNames(list(mean(analysis_df[[time_column]])), time_column))
```

Marginalising over a single time point (the mean day) to summarise treatment effects from a `Day * Treatment` interaction model discards the full interaction structure. If treatments diverge over time (the biologically relevant signal), the marginal mean at a single time point is an incomplete summary. The user should receive estimated marginal means across multiple time points or the interaction terms themselves should be the primary output.

**Resolution:** The lme4 path now computes two EMM objects: the existing `lsmeans_obj` at the mean day (backward-compatible `treatment_effects`) and a new `lsmeans_time` at the five quintile study days (min, Q1, median, Q3, max). The latter is returned as `treatment_effects_over_time` — a data frame with one row per (Treatment, Day) combination showing the adjusted mean, SE, and 95% CI at each quintile. The pairwise contrasts continue to use the mean-day EMM for a compact summary.

---

### 2.8 Survival analysis: logrank p-value inflated for multi-group comparison — FIXED
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

**Resolution:** The logrank fallback now runs a pairwise log-rank test for each non-reference group (`survdiff()` restricted to that group vs. reference, df = 1) and assigns the resulting p-value to that group's row only. The omnibus test result is retained in a local variable `omnibus_p` for reference but is not propagated to the per-group results table.

---

### 2.9 `therapeutic_window_metric()` — max weight loss taken as maximum across group, not mean — FIXED
**File:** `R/therapeutic_window_metric.R` (~lines 87–89)

```r
group_wl <- stats::aggregate(Pct_Loss ~ Treatment, data = mouse_wl, FUN = max, na.rm = TRUE)
```

The TWM denominator is the single worst-case mouse in each group. This makes TWM hypersensitive to a single outlier and not representative of typical toxicity. A mean or median weight loss would be more interpretable.

**Resolution:** Changed `FUN = max` to `FUN = mean`. The column renamed from `Max_Pct_Weight_Loss` to `Mean_Pct_Weight_Loss`. Dashboard labels and the broken column references in the summary text (`TGI_Pct` → `TGI`, `Max_Weight_Loss_Pct` → `Mean_Pct_Weight_Loss`) corrected. TWM formula description updated throughout the dashboard UI.

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

---

## Bayesian Analysis Functions Review (v0.3.7 – v0.4.2)

**Reviewer:** Claude Code
**Date:** 2026-05-20
**Scope:** `bayesian_tumor_growth()`, `bayesian_survival()`, `bayesian_body_weight()`, `bayesian_dose_response()`, `bayesian_therapeutic_window()`, `bayesian_synergy()`

---

### B1. API Inconsistencies

#### B1.1 `model_type_used` is not consistent across functions
**File:** `R/bayesian_tumor_growth.R:579`

`bayesian_tumor_growth()` returns `model_type_used = "bayes"` while all other Bayesian functions return `"bayes_bw"`, `"bayes_survival"`, `"bayes_dr"`, `"bayes_twm"`, `"bayes_synergy"`. Dashboard and any downstream code that dispatch on `model_type_used` must special-case tumor growth. Should be `"bayes_tg"` to complete the pattern.

**Fix:** Change line 579 of `bayesian_tumor_growth.R` from `"bayes"` to `"bayes_tg"` and update the dashboard handler that checks this value.

---

#### B1.2 `prior_strength` default differs: `bayesian_body_weight` defaults to `"weakly_informative"`
**File:** `R/bayesian_body_weight.R:116`

```r
# body_weight — "weakly_informative" is first (and therefore the match.arg default):
prior_strength = c("weakly_informative", "skeptical", "informative", "diffuse", "manual")

# all other Bayesian functions — "skeptical" is first:
prior_strength = c("skeptical", "weakly_informative", "informative", "diffuse", "manual")
```

A user calling both `bayesian_tumor_growth()` and `bayesian_body_weight()` with all defaults gets different priors without any indication. Produces a subtle systematic difference in results.

**Fix:** Reorder the body_weight default to put `"skeptical"` first, matching the rest of the family.

---

#### B1.3 `"informative"` preset missing from `bayesian_dose_response`
**File:** `R/bayesian_dose_response.R:123`

`prior_strength` options are `c("skeptical", "weakly_informative", "diffuse", "manual")` — `"informative"` is absent. All other Bayesian functions include it. Calling `bayesian_dose_response(..., prior_strength = "informative")` throws a cryptic `match.arg` error rather than a helpful message.

**Fix:** Add `"informative"` to the choice vector and to the switch statement at line ~208.

---

#### B1.4 Credible interval columns named `Lower_CL`/`Upper_CL` in most functions, `Lower_CrI`/`Upper_CrI` in `bayesian_dose_response`
**Files:** `bayesian_tumor_growth.R:375`, `bayesian_body_weight.R:352`, `bayesian_survival.R:464`, `bayesian_dose_response.R:345`

`CL` is "confidence limit" — frequentist terminology. `CrI` is "credible interval" — the correct Bayesian term. `bayesian_dose_response` got it right; the others did not. Dashboard code that accesses these columns by name will break if column names are ever harmonised.

**Fix:** Rename to `Lower_CrI`/`Upper_CrI` across `bayesian_tumor_growth`, `bayesian_body_weight`, and `bayesian_survival`.

---

#### B1.5 `include_cage_effect` vs `include_frailty` — same concept, different argument name
**Files:** `bayesian_tumor_growth.R:177`, `bayesian_body_weight.R:125`, `bayesian_synergy.R:154` vs `bayesian_survival.R:149`

The survival function uses `include_frailty` for what is conceptually the same feature (cage random intercept). Users must remember a different parameter name for one function.

**Fix:** Rename `include_frailty` to `include_cage_effect` in `bayesian_survival.R` and update documentation.

---

#### B1.6 Intercept prior scaling inconsistent: `b_sd * 2.5` vs `b_sd * 2.0`
**Files:** `bayesian_tumor_growth.R:299`, `bayesian_body_weight.R:261`, `bayesian_survival.R` vs `bayesian_synergy.R:287`

`bayesian_tumor_growth`, `bayesian_body_weight`, and `bayesian_survival` set the intercept prior width to `b_sd * 2.5`. `bayesian_synergy` uses `b_sd * 2`. This creates a slightly different effective prior on the intercept scale for no documented reason.

**Fix:** Standardise to `b_sd * 2.5` in `bayesian_synergy.R:287`.

---

#### B1.7 `bayesian_synergy` uses `brms::set_prior()` while all other functions use `brms::prior_string()`
**Files:** `bayesian_synergy.R:276-293` vs all other `bayesian_*.R`

Both work, but the inconsistency is confusing in code review and increases the diff noise for future maintenance. `prior_string()` is the simpler, string-only variant; `set_prior()` is an alias for the more flexible `prior()`. There is no reason to use the more complex form here.

**Fix:** Replace `brms::set_prior(...)` with `brms::prior_string(...)` in `bayesian_synergy.R`.

---

### B2. Code Redundancy

#### B2.1 Back-transform helper `bt()` defined independently in three files
**Files:** `bayesian_body_weight.R:430`, `bayesian_synergy.R:346`, `bayesian_therapeutic_window.R:206`

All three define a local function:
```r
bt <- function(x) { if (transform == "log") exp(x) else if (transform == "sqrt") x^2 else x }
```
The therapeutic_window version additionally accepts a `transform` argument (the others close over it). Three copies to maintain.

**Fix:** Extract to `R/utils_bayes.R` as an exported-internal helper:
```r
#' @noRd
bayes_backtransform <- function(x, transform) {
  switch(transform, log = exp(x), sqrt = x^2, x)
}
```

---

#### B2.2 MCMC diagnostics data frame constructed identically in four files
**Files:** `bayesian_tumor_growth.R:334`, `bayesian_body_weight.R:304`, `bayesian_survival.R:308`, `bayesian_dose_response.R:319`

```r
mcmc_diagnostics <- data.frame(
  Parameter = posterior_summary$Parameter,
  Rhat      = round(posterior_summary$Rhat,     4),
  Bulk_ESS  = round(posterior_summary$Bulk_ESS, 0),
  Tail_ESS  = round(posterior_summary$Tail_ESS, 0),
  Converged = posterior_summary$Rhat <= 1.01,
  stringsAsFactors = FALSE
)
```
`bayesian_synergy` uses a different approach (`brms::rhat()` + `brms::neff_ratio()` separately, lines 322-337) and produces an inconsistent `mcmc_diagnostics` schema that lacks `Bulk_ESS` and `Tail_ESS` columns.

**Fix:** Extract to `R/utils_bayes.R`:
```r
#' @noRd
make_mcmc_diagnostics <- function(posterior_summary_df) { ... }
```
Use it in all six functions.

---

#### B2.3 Prior specification switch duplicated across five files
**Files:** `bayesian_tumor_growth.R:285`, `bayesian_body_weight.R:244`, `bayesian_survival.R:247`, `bayesian_synergy.R:259`, `bayesian_dose_response.R:208`

Each function contains the same lookup table:
```r
list(
  skeptical          = list(b_sd = 0.25, exp_rate = 2),
  weakly_informative = list(b_sd = 1.0,  exp_rate = 1),
  informative        = list(b_sd = 0.5,  exp_rate = 2),
  diffuse            = list(b_sd = 2.5,  exp_rate = 0.5)
)
```
Changing a default (e.g., tightening the skeptical prior) requires touching five files and risks drift.

**Fix:** Centralise as `bayes_prior_params(prior_strength)` in `R/utils_bayes.R`.

---

#### B2.4 Cage placeholder setup repeated in four files
**Files:** `bayesian_tumor_growth.R:206`, `bayesian_body_weight.R:151`, `bayesian_survival.R:185`, `bayesian_synergy.R:227`

All four check for a cage column, fall back to `id_column` as a placeholder, and construct `id_var` as `cage_var + "__" + id`. The pattern is identical except for variable names used.

**Fix:** Extract to `R/utils_bayes.R`:
```r
#' @noRd
setup_cage_vars <- function(df, cage_column, id_column) { ... }
```

---

### B3. Bugs and Correctness Issues

#### B3.1 Loewe denominator floor (1e-6) is arbitrary and undocumented
**File:** `bayesian_synergy.R:413`

```r
loewe_ci <- loewe_num / pmax(fe_combo, 1e-6)
```

When `fe_combo` is near zero (minimal combo efficacy), CI values reach ~1e6 — technically not infinity but practically meaningless. No documentation explains why 1e-6 was chosen or how to interpret CI when combo effect is near zero. This will silently produce enormous CI values in negative-result studies.

**Fix:** Replace with a threshold relative to the maximum observed FE:
```r
fe_floor <- max(fe_combo) * 1e-4
fe_floor <- max(fe_floor, 1e-4)   # absolute minimum guard
loewe_ci <- loewe_num / pmax(fe_combo, fe_floor)
```
Document the threshold and add a note to `loewe_summary` when floor was applied.

---

#### B3.2 `bayesian_synergy` MCMC diagnostics lack `Bulk_ESS` and `Tail_ESS`
**File:** `bayesian_synergy.R:317-337`

All other Bayesian functions derive `mcmc_diagnostics` from `fixed_df` (the formatted posterior summary that includes Rhat, Bulk_ESS, Tail_ESS). `bayesian_synergy` instead calls `brms::rhat()` and `brms::neff_ratio()` separately and only appends `ESS_Bulk` (singular) — `Tail_ESS` is never populated. The schema differs from every other function's `mcmc_diagnostics`.

**Fix:** Follow the same pattern as the other functions: extract `fixed_df` from the brms summary, build `mcmc_diagnostics` from it with the shared helper (see B2.2).

---

#### B3.3 Zero/negative volume input with log transform: no guard against all-negative data
**Files:** `bayesian_tumor_growth.R:240`, `bayesian_body_weight.R:194`, `bayesian_synergy.R:202`

The zero-fill pattern (`vol[vol <= 0] <- min(vol[vol > 0], na.rm = TRUE) / 2`) will set `vol` to `Inf / 2 = Inf` if all values are ≤ 0, because `min(numeric(0)) = Inf`. The resulting `log(Inf)` will silently corrupt the model fit.

**Fix:** Add a guard before the fill:
```r
positive_vals <- vol[is.finite(vol) & vol > 0]
if (length(positive_vals) == 0L) {
  stop("No positive values found in '", volume_column,
       "' after filtering. Cannot apply log transform.")
}
```

---

#### B3.4 `bayesian_therapeutic_window`: common groups check uses `intersect(treated_tg, treated_bw)` but does not verify alignment of `bw_groups` factor levels
**File:** `bayesian_therapeutic_window.R:225`

The function assumes that group names from the TG model (`tg_groups`) and BW model (`bw_groups`) are in directly comparable character form. If one model was fitted on a data frame where a group name has trailing whitespace or different capitalisation, `intersect()` returns an empty vector and the function stops with a misleading message about "no common treated groups." There is no normalisation step.

**Fix:** Trim and `tolower()` both name vectors before the `intersect()` call (or warn early that names differ only in whitespace/case).

---

#### B3.5 Rounding inconsistency in `bliss_summary` return: excess uses 4 decimal places, FE medians use 3
**File:** `bayesian_synergy.R:402-409`

```r
Excess_Median      = round(bliss_q["50%"],       4),   # 4 dp
Expected_FE_Median = round(bliss_fe_q["50%"],    3),   # 3 dp
Observed_FE_Median = round(obs_fe_med,           3)    # 3 dp
```

A user comparing `Excess_Median` with `Observed_FE_Median - Expected_FE_Median` will get different values at the 4th decimal place, which is confusing.

**Fix:** Use 4 decimal places consistently throughout `bliss_summary` (or 3 — pick one and be consistent).

---

### B4. Documentation Gaps

#### B4.1 `bayesian_tumor_growth` reference group auto-detection logic is not documented
**File:** `R/bayesian_tumor_growth.R`

The `@param reference_group` line says "Auto-detected if NULL" but does not describe the detection strategy (check for "Control"/"Vehicle"/"control" etc., then first alphabetically). Users cannot predict which group will be the reference without reading source code.

**Fix:** Add to `@param reference_group`: "Auto-detection checks for common control-group names (`Control`, `Vehicle`, `control`, `vehicle`, `CTRL`, `ctrl`) before falling back to the first level alphabetically."

---

#### B4.2 `bayesian_body_weight` does not document the weight-loss calculation endpoints
**File:** `R/bayesian_body_weight.R`

`weight_loss_summary` is computed as percentage change from the first to last study day in the model data. This is not stated explicitly — users may not realise that animals with different dropout patterns will have different "first day" or "last day" references, and that the summary always uses the overall study start and end days, not individual animal endpoints.

**Fix:** Add to `@return` documentation: "Weight-loss percentages are computed from the earliest to the latest study day present in the model data. This is a population-level prediction; individual dropout patterns do not affect the endpoint days used."

---

#### B4.3 `bayesian_synergy` Loewe CI limitations are less prominently documented than in `analyze_drug_synergy`
**File:** `R/bayesian_synergy.R`

`analyze_drug_synergy()` has a full `@section Assumptions and Limitations:` added in the prior code review. `bayesian_synergy()` describes the formula in `@details` but does not carry the same explicit "single-dose approximation assumes linear dose-response" warning. The Bayesian wrapper adds credible intervals to the CI but the underlying approximation is the same.

**Fix:** Add a matching `@section Assumptions and Limitations:` to `bayesian_synergy()`, referencing the same Bliss ceiling and Loewe linear-assumption caveats.

---

#### B4.4 `bayes_prior_posterior_plot()` is internal but used by three functions; its location is not obvious
**File:** `R/bayesian_tumor_growth.R` (end of file)

The helper is defined at the bottom of `bayesian_tumor_growth.R` with `@noRd`. Both `bayesian_body_weight` and `bayesian_survival` call it, so it is effectively a cross-file dependency. A developer reading `bayesian_body_weight.R` has no indication where `bayes_prior_posterior_plot()` is defined.

**Fix:** Move to `R/utils_bayes.R` with a comment indicating its consumers, and add a comment in each consumer file pointing to `utils_bayes.R`.

---

### B5. Performance

#### B5.1 `bayesian_body_weight` calls `posterior_epred()` three times
**File:** `R/bayesian_body_weight.R:416-428` and `:577-583`

Three separate `posterior_epred()` calls are made: first day weights, last day weights, and full trajectory. Each call re-runs Stan memory allocation and prediction. The first two could be combined into a single call with a two-row `newdata` matrix, halving one brms round-trip.

**Fix:**
```r
nd_endpoints <- rbind(
  make_nd(bw_groups, study_days[1L]),
  make_nd(bw_groups, study_days[length(study_days)])
)
ep_ends <- brms::posterior_epred(model, newdata = nd_endpoints,
                                  re_formula = NA, ...)
ep_first <- ep_ends[, seq_len(n_groups)]
ep_last  <- ep_ends[, seq_len(n_groups) + n_groups]
```

---

### B6. Missing Functionality

#### B6.1 No `bayesian_power_analysis()` function
The package has both analytic (`apriori_power_analysis`) and simulation (`apriori_power_simulation`) frequentist power analyses. There is no Bayesian equivalent — no way to estimate the required N to achieve a target posterior probability of a given effect size exceeding a threshold. This is the natural next step for a fully Bayesian workflow.

#### B6.2 `bayesian_therapeutic_window` cannot be used without two separately fitted models
**File:** `R/bayesian_therapeutic_window.R`

`bayesian_therapeutic_window()` requires pre-fitted `brmsfit` objects from `bayesian_tumor_growth()` and `bayesian_body_weight()`. There is no convenience wrapper that fits both from a combined data frame. This is a high friction path for exploratory analysis.

**Suggestion:** Add a `bayesian_therapeutic_window_from_data()` convenience function, or add a `tg_df` + `bw_df` path that internally calls both fitting functions before computing TWM.

#### B6.3 No uncertainty-aware version of `analyze_drug_synergy_over_time()`
`bayesian_synergy()` computes synergy at a single endpoint day. The time-series counterpart (`bayesian_synergy_over_time()`) would compute draw-wise Bliss and Loewe CIs at each study day, showing how synergy evolves. This is directly analogous to `analyze_drug_synergy_over_time()`.

---

### B7. Harmonization Opportunities

#### B7.1 Treatment effects table schema is not consistent between Bayesian and frequentist functions
`bayesian_tumor_growth()` and `bayesian_body_weight()` return `treatment_effects` with columns `Adjusted_Mean`, `SE`, `DF`, `Lower_CL`, `Upper_CL`, `Note` — sourced from `emmeans`. The frequentist `tumor_growth_statistics()` returns the same schema. However, `bayesian_survival()` returns `Time_Ratio`, `Lower_CL`, `Upper_CL`, `HR`, `Median_Survival` — a completely different schema. Dashboard modules that render treatment_effects tables must handle each function type separately.

**Suggestion:** Define a canonical schema document and have each function document its deviations explicitly.

#### B7.2 `bayesian_survival` is the only Bayesian function not returning `transform_used`
All other Bayesian functions return `transform_used` in the result list. `bayesian_survival` does not — it has `family_used` and `frailty_used` instead. This is appropriate since survival models do not use a volume transform, but the absence of the field means `bayesian_therapeutic_window` cannot be applied to survival results (and the documentation should make this explicit).

---

### B8. Quick-Reference Checklist — Bayesian Functions

| # | Issue | Severity | File(s) | Status |
|---|-------|----------|---------|--------|
| B1.1 | `model_type_used = "bayes"` should be `"bayes_tg"` | Major | `bayesian_tumor_growth.R:579` | ✅ Fixed |
| B1.2 | `prior_strength` default is `"weakly_informative"` in body_weight, `"skeptical"` elsewhere | Major | `bayesian_body_weight.R:116` | ✅ Fixed |
| B1.3 | `"informative"` preset missing from `bayesian_dose_response` | Major | `bayesian_dose_response.R:123` | ✅ Fixed |
| B1.4 | `Lower_CL`/`Upper_CL` should be `Lower_CrI`/`Upper_CrI` in Bayesian functions | Major | TG, BW, Survival | ✅ Fixed |
| B1.5 | `include_frailty` vs `include_cage_effect` — same concept, different name | Minor | `bayesian_survival.R:149` | ✅ Fixed |
| B1.6 | Intercept prior width: `b_sd * 2.5` in 3 files, `b_sd * 2.0` in synergy | Minor | `bayesian_synergy.R:287` | ✅ Fixed |
| B1.7 | `set_prior()` in synergy vs `prior_string()` everywhere else | Minor | `bayesian_synergy.R:276` | ✅ Fixed |
| B2.1 | `bt()` back-transform helper duplicated in 3 files | Minor | BW, Synergy, TWM | ✅ Fixed — `R/utils_bayes.R` |
| B2.2 | MCMC diagnostics DF construction duplicated in 4+ files; schema differs in synergy | Minor | Multiple | ✅ Fixed — `R/utils_bayes.R` |
| B2.3 | Prior specification switch duplicated in 5 files | Minor | Multiple | ✅ Fixed — `R/utils_bayes.R` |
| B2.4 | Cage placeholder setup duplicated in 4 files | Minor | Multiple | ✅ Fixed — `R/utils_bayes.R` |
| B3.1 | Loewe 1e-6 floor arbitrary and undocumented; produces huge CI for near-zero combo effects | Major | `bayesian_synergy.R:413` | ✅ Fixed |
| B3.2 | `bayesian_synergy` MCMC diagnostics missing Bulk_ESS and Tail_ESS | Minor | `bayesian_synergy.R:317` | ✅ Fixed |
| B3.3 | All-zero/negative volume with log transform silently produces `log(Inf)` | Major | TG, BW, Synergy | ✅ Fixed |
| B3.4 | TWM group-name intersection not normalised for whitespace/case | Minor | `bayesian_therapeutic_window.R:225` | ✅ Fixed |
| B3.5 | `bliss_summary` inconsistent rounding (3 vs 4 dp) | Minor | `bayesian_synergy.R:402` | ✅ Fixed |
| B4.1 | Reference-group auto-detection not documented | Minor | `bayesian_tumor_growth.R` | ✅ Fixed |
| B4.2 | Weight-loss endpoints not documented in `bayesian_body_weight` | Minor | `bayesian_body_weight.R` | ✅ Fixed |
| B4.3 | `bayesian_synergy` missing `@section Assumptions and Limitations:` | Minor | `bayesian_synergy.R` | ✅ Fixed |
| B4.4 | `bayes_prior_posterior_plot()` cross-file dependency not indicated | Minor | `bayesian_tumor_growth.R` | ✅ Fixed — moved to `R/utils_bayes.R` |
| B5.1 | `bayesian_body_weight` calls `posterior_epred()` 3× instead of 2× | Minor | `bayesian_body_weight.R:416` | ✅ Fixed |
| B6.1 | No `bayesian_power_analysis()` function | Enhancement | — | ✅ Fixed — `R/bayesian_power_analysis.R` |
| B6.2 | No single-call wrapper for TWM (requires two pre-fitted models) | Enhancement | `bayesian_therapeutic_window.R` | ✅ Fixed — `bayesian_twm_from_data()` |
| B6.3 | No `bayesian_synergy_over_time()` | Enhancement | — | ✅ Fixed — `R/bayesian_synergy.R` |
| B7.1 | Treatment effects table schema not consistent across Bayesian/frequentist | Minor | Multiple | Open |
| B7.2 | `bayesian_survival` does not return `transform_used` | Minor | `bayesian_survival.R` | Open |

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
| 2.2 | Simulation generates linear-scale data; analysis fits log-scale — anti-conservative power | Major | `apriori_power_simulation.R` | ✅ Fixed |
| 2.3 | d→f conversion for ANOVA is scenario-specific, not general | Major | `apriori_power_analysis.R` | ✅ Fixed |
| 2.4 | Bliss/Loewe assumptions on TGI not documented | Major | `analyze_drug_synergy*.R` | ✅ Fixed |
| 2.5 | LOCF AUC returns volume (mm³), not area (mm³·day) | Major | `tumor_auc_analysis.R` | ✅ Fixed |
| 2.6 | AIC comparison between `lm` and `drc` models is invalid | Major | `dose_response_statistics.R` | ✅ Fixed |
| 2.7 | emmeans at single mean time point discards interaction | Major | `tumor_growth_statistics.R` | ✅ Fixed |
| 2.8 | Omnibus logrank p-value assigned as per-group p-value | Major | `survival_statistics.R` | ✅ Fixed |
| 2.9 | TWM denominator is max of one mouse, not group mean | Major | `therapeutic_window_metric.R` | ✅ Fixed |
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

---

# Review Round 2 (v0.3.6 — 2026-05-30)

**Reviewer:** Claude Code
**Scope:** Statistical validity, model diagnostics, code architecture; changes since the v0.3.3 review

## Round-1 status

Spot-checked the v0.3.3 critical / major / minor items table above — all 25+ items remain resolved at v0.3.6. The Bayesian B1–B8 harmonisation items also appear resolved (consistent `Lower_CrI`/`Upper_CrI` naming, shared `utils_bayes.R` for `make_mcmc_diagnostics`, `bayes_prior_params`, `setup_cage_column`, `bayes_prior_posterior_plot`, `build_posterior_summary`). The package is materially more robust than at the prior review.

This second round focuses on **statistical-validity gaps that the first review did not surface** — primarily missing model-fit diagnostics for both the Cox PH path and every brms model — plus a few new bugs and architectural recommendations.

---

## A. Critical — Statistical Validity

### A.1 No proportional-hazards check in `survival_statistics()`
**File:** `R/survival_statistics.R`

Cox proportional hazards is valid only if the PH assumption holds. The package never tests it. `survival::cox.zph()` returns Schoenfeld-residual–based tests per covariate plus a global test; a single line after fitting `cox_model` would surface this. Without it, hazard ratios reported as point estimates with CIs can be materially misleading when PH is violated (a common situation when treatment effects accelerate or attenuate over time).

**Recommendation:**
```r
ph_test <- tryCatch(survival::cox.zph(cox_model), error = function(e) NULL)
# add to result_list:
result_list$ph_test <- ph_test
# expose Schoenfeld residual plot via plot(ph_test) for dashboard
```
Flag a warning when any p < 0.05 in the test table and surface the result in the dashboard's Survival results panel.

### A.2 No NUTS / divergent-transitions diagnostics in any Bayesian function
**Files:** all `R/bayesian_*.R`

Every Bayesian function returns `mcmc_diagnostics` containing only `Rhat`, `Bulk_ESS`, `Tail_ESS`, and a `Converged = Rhat <= 1.01` flag. The Stan/NUTS sampler also reports:
- **Divergent transitions** — geometry pathology; any non-zero count signals untrustworthy posterior near divergent points.
- **`max_treedepth` hits** — efficiency warning; high counts mean the sampler is exploring inefficiently.
- **E-BFMI < 0.3** — energy diagnostic; low values indicate the sampler is having trouble exploring the typical set.

All three are first-line Stan diagnostics. None are checked. A model with `Rhat = 1.00` but 200 divergences is **not** converged in any meaningful sense.

**Recommendation:** extend `make_mcmc_diagnostics()` in `utils_bayes.R`:
```r
np <- brms::nuts_params(model)
div_count <- sum(np$Value[np$Parameter == "divergent__"])
max_td    <- max(np$Value[np$Parameter == "treedepth__"])
ebfmi     <- brms::neff_ratio(model)  # or rstan::get_bfmi(model$fit)
```
Surface these in the returned list and flag clearly in the dashboard.

### A.3 `auc_method` argument silently ignored in `tumor_growth_statistics()`
**File:** `R/tumor_growth_statistics.R` lines 682, 832, and `R/utils_auc.R` line 41

`tumor_growth_statistics()` exposes `auc_method = c("trapezoidal", "last_observation")` and `match.arg`s it, but the chosen method is **never passed** to `tgs_compute_auc()`, which calls `calculate_auc(time, volume)` — a trapezoidal-only function. Users requesting `auc_method = "last_observation"` get the trapezoidal result with no warning. This is the same class of bug the v0.3.3 review caught for `handle_cage_effects` / `random_effects_specification`.

**Recommendation:** either add a `method` arg to `calculate_auc()` and plumb it through, or remove the `auc_method` arg from `tumor_growth_statistics()` and document that LOCF AUC is available via `tumor_auc_analysis()` instead.

### A.4 No LOO-CV / WAIC / Pareto-k in any Bayesian function
**Files:** all `R/bayesian_*.R`

`brms::loo()` runs PSIS-LOO cross-validation cheaply on a fitted model and returns (a) an `elpd_loo` for model comparison, (b) Pareto-k shape parameters per observation — values > 0.7 flag influential observations that LOO can't reliably approximate. These are the standard Bayesian counterparts to AIC/BIC and Cook's distance, and they're trivial to add. Without them, the package has no tool for either Bayesian model comparison or detecting influential mice.

**Recommendation:** add `bayes_loo()` to `utils_bayes.R`:
```r
bayes_loo <- function(model) {
  loo_obj <- tryCatch(brms::loo(model, save_psis = TRUE), error = function(e) NULL)
  if (is.null(loo_obj)) return(NULL)
  list(
    elpd_loo = loo_obj$estimates["elpd_loo", "Estimate"],
    p_loo    = loo_obj$estimates["p_loo",    "Estimate"],
    pareto_k = loo_obj$diagnostics$pareto_k,
    n_high_k = sum(loo_obj$diagnostics$pareto_k > 0.7)
  )
}
```
Surface `pareto_k` as a per-mouse table in the dashboard.

---

## B. Bugs

### B.1 `tgs_extrapolate()` extrapolation count is always 0 or 1
**File:** `R/tumor_growth_statistics.R` line 99

```r
extrapolated_subjects <- unique(df$Extrapolated[df$Extrapolated])
n_extrapolated <- length(extrapolated_subjects)
```
`df$Extrapolated[df$Extrapolated]` returns the subset of rows where `Extrapolated == TRUE` — i.e. a logical vector of `TRUE`s. `unique()` collapses that to `c(TRUE)`, so `n_extrapolated` is always `1` whenever any extrapolation happened (and `0` otherwise). The "Successfully extrapolated N subjects" message is therefore meaningless. Should count **unique subject keys** that have any extrapolated rows.

**Recommendation:**
```r
ex_keys <- unique(make_mouse_key(
  df[df$Extrapolated, id_column],
  df[df$Extrapolated, treatment_column],
  df[df$Extrapolated, cage_column]
))
n_extrapolated <- length(ex_keys)
```

### B.2 Magic `1.96` for survival CI construction
**File:** `R/survival_statistics.R` lines 525–526

```r
ci_lower <- exp(model_summary$coefficients[i, "coef"] - 1.96 * ...)
ci_upper <- exp(model_summary$coefficients[i, "coef"] + 1.96 * ...)
```
Standard Cox `summary()` already gives `exp(coef)`, `lower .95`, `upper .95` columns — use those directly. The magic-number form is correct numerically but bypasses the proper `qnorm(0.975)` (1.959964) and obscures the level.

### B.3 `lme4` necrosis helper coupled to its file
**File:** `R/tumor_growth_statistics.R` (defines `tgs_handle_necrosis`), `R/bayesian_tumor_growth.R` line 237 (uses it)

`tgs_handle_necrosis()` is called from the Bayesian path too. It's correctly DRY in implementation but its location is misleading — a future developer trimming the frequentist file may delete it. Move to `utils_necrosis.R` (already exists at 64 LOC).

### B.4 `bayes_prior_params()` lacks a default case
**File:** `R/utils_bayes.R` line 50

`switch(prior_strength, skeptical = ..., ...)` has no default; passing `"manual"` returns `NULL` and any downstream `pp$b_sd` errors with "$ operator is invalid for atomic vectors". Currently safe because every caller checks `prior_strength == "manual"` before invoking, but defensive code would add:
```r
stop("Unknown prior_strength: ", prior_strength)
```

---

## C. Diagnostics — Present but Incomplete

| Path | Present | Missing |
|---|---|---|
| `tumor_growth_statistics` (LME4) | residuals, fitted, QQ, random effects, variance components, R² in growth rates, cage chi-square | homoscedasticity test (Breusch-Pagan), autocorrelation (Durbin-Watson on residuals vs time), influence diagnostics (`influence.merMod`), normality of random effects |
| `tumor_growth_statistics` (AUC) | residuals, fitted, QQ from `aov` | Levene's test for variance equality (Welch addresses this but isn't tested), normality test |
| `survival_statistics` | event/total counts, event rates, median survival, separation check, cage distribution | `cox.zph()` (A.1), Schoenfeld residuals, C-index, deviance residuals for influence |
| `bayesian_*` | Rhat, Bulk_ESS, Tail_ESS, pp_check, prior_posterior, trace, residuals_vs_day | divergences (A.2), max_treedepth, E-BFMI, LOO+Pareto-k (A.4), `bayes_R2()` |

### C.1 `bayes_R2()` not exposed
`brms::bayes_R2(model)` returns a posterior distribution of R². For LMM-style models it's the natural Bayesian summary of explained variance. Add one line per Bayesian function:
```r
result$bayes_R2 <- as.data.frame(brms::bayes_R2(model, summary = TRUE))
```

### C.2 Posterior P(effect > 0) is computed nowhere
A first-class advantage of Bayesian analysis over null-hypothesis testing is directional posterior probability: `mean(draws > 0)`. For every treatment contrast we already compute, adding `P_direction` is one line. Replaces the awkward "is the 95% CrI excluding zero?" interpretation with a quantitative posterior statement.

### C.3 ESS / N ratio
`Bulk_ESS = 4000` sounds great until you realise you sampled `n_chains × n_iter = 8000`. The ratio (efficiency) is what matters. Stan recommends ESS/N > 0.1; lower indicates an inefficient parameterisation. Easy to add.

---

## D. Architecture

### D.1 `tumor_growth_statistics.R` is 1,218 LOC in one file
The main function is ~550 LOC with two large mutually-exclusive return paths (`auc` vs `lme4`). Splitting into `tgs_path_lme4()` / `tgs_path_auc()` helpers and moving `tgs_*` helpers into `tgs_helpers.R` would make the file navigable. Same applies to `bayesian_synergy.R` (1,046 LOC) and `bayesian_body_weight.R` (714 LOC).

### D.2 Function signatures with 20+ parameters
`tumor_growth_statistics()` takes 23 arguments; `bayesian_tumor_growth()` takes 24. Consider grouping related parameters into a config helper:
```r
tg_priors(strength = "skeptical", b = NULL, intercept = NULL, sd = NULL, sigma = NULL)
tg_mcmc(chains = 4, warmup = 1000, iter = 500, seed = 42)
tumor_growth_statistics(df, columns, priors = tg_priors(), mcmc = tg_mcmc(), ...)
```
This is a soft recommendation — invasive to existing callers — but would substantially reduce the visual weight of every Bayesian entry point.

**✅ Resolved in v0.4.7** — additive implementation: `tg_priors()` and
`tg_mcmc()` are exported helpers (`R/bayesian_config.R`); `priors` and
`mcmc` are accepted by `bayesian_tumor_growth()`, `bayesian_body_weight()`,
`bayesian_survival()` (priors$sigma → prior_aux mapping), `bayesian_synergy()`,
and `bayesian_dose_response()` (mcmc only — DR's Hill priors are
model-specific). All individual `prior_*` / `n_chains` / etc. arguments
continue to work; the helper just overrides them when supplied. New
callers should prefer the helpers; existing callers keep working with no
changes.

### D.3 Plot generation is intermixed with statistical computation
Many statistical functions return ggplot objects directly. This couples the analysis backend to ggplot2 (a heavy dep) and makes the functions hard to test on a headless CI. Consider returning **data frames** plus thin `plot_*()` helpers that consume them — the package already has `plot_*.R` files; the pattern just isn't fully applied.

**✅ Resolved in v0.4.7** — every `bayesian_*` analysis function now
universally respects `plots = FALSE` (verified — each guards its plot
block with `if (isTRUE(plots))`), and the result list cleanly separates
data fields (`treatment_effects`, `posterior_summary`, …) from plot
fields (`pp_check_plot`, `credible_intervals_plot`, …). Passing
`plots = FALSE` returns the data-only path suitable for headless / CI
pipelines, satisfying the original concern. Pattern documented in
`bayesian_tumor_growth()` via a new `@section Separating analysis from
visualization`. The Shiny dashboard already uses this pattern — it
passes `plots = FALSE` and rebuilds plots reactively from the data.

---

## E. Recommended Additional Functionality

### E.1 `cmdstanr` backend option
RStan is the default; `cmdstanr` (a CmdStan wrapper) is 3–10× faster for compilation, has better diagnostic reporting, and is the actively-maintained Stan interface. Add `backend = c("rstan", "cmdstanr")` to all `bayesian_*` functions and pass to `brms::brm(backend = ...)`. On the VPS this would meaningfully reduce the 3–12 min fit time.

### E.2 Concordance (C-index) for survival
`survival::concordance(cox_model)$concordance` gives the C-index — the survival analogue of AUC-ROC. Should be reported in every `survival_statistics` result.

### E.3 Per-animal posterior growth rates
`bayesian_tumor_growth()` currently returns per-animal growth rates from **OLS** on log-volumes, not from the brms model itself (documented as such, line 119). Extracting random-slope posterior draws (when `random_effects_specification = "slope"`) would give per-animal posterior distributions, including credible intervals.

### E.4 Bootstrap CIs for AUC treatment effects
The AUC path uses Welch's t-tests on per-animal AUCs and reports parametric CIs. Per-animal AUC is itself a derived statistic with its own uncertainty (especially when extrapolated). Bootstrap CIs would be more honest.

### E.5 Power analysis: surface effect-size scale used
`apriori_power_analysis()` uses Cohen's d on the **modelling scale** (e.g. log-volume), not on the raw response scale. This is the correct choice statistically but trips up users who specify "d = 0.5 means a 0.5 mm³ difference". Adding a clear note in the result describing the d-units would prevent silent misuse.

### E.6 Posterior predictive interval coverage
A standard check: simulate from the posterior predictive distribution and compute the empirical coverage of the 50%, 80%, 95% intervals against held-out data. If coverage is far from nominal, the model is mis-specified. `posterior::summarise_draws()` makes this easy.

---

## F. Status Summary — Round 2

| ID | Issue | Severity | File | Status |
|---|---|---|---|---|
| A.1 | No `cox.zph()` proportional hazards check | Critical | `survival_statistics.R` | ✅ Fixed v0.4.5 |
| A.2 | No NUTS divergences / max_treedepth / E-BFMI | Critical | all `bayesian_*.R` | ✅ Fixed v0.4.5 (`make_nuts_diagnostics`) |
| A.3 | `auc_method` parameter silently ignored | Critical | `tumor_growth_statistics.R` | ✅ Fixed v0.4.5 (removed) |
| A.4 | No LOO-CV / Pareto-k / WAIC | Major | all `bayesian_*.R` | ✅ Fixed v0.4.5 (`bayes_loo`) |
| B.1 | Extrapolation count always 0 or 1 | Major | `tumor_growth_statistics.R` | ✅ Fixed v0.4.5 |
| B.2 | Magic 1.96 for Cox CI | Minor | `survival_statistics.R` | ✅ Fixed v0.4.5 |
| B.3 | `tgs_handle_necrosis` misplaced | Minor | `tumor_growth_statistics.R` | Retracted (H.1) |
| B.4 | `bayes_prior_params` missing default | Minor | `utils_bayes.R` | ✅ Fixed v0.4.5 |
| C.1 | No `bayes_R2()` | Major | all `bayesian_*.R` | ✅ Fixed v0.4.5 (`bayes_r2_summary`) |
| C.2 | No posterior P(effect > 0) | Major | all `bayesian_*.R` | ✅ Fixed v0.4.5 (`emm_p_direction` for TG/BW; survival/DR open) |
| C.3 | No ESS/N efficiency ratio | Minor | `utils_bayes.R` | ✅ Fixed v0.4.5 |
| D.1 | 1,000+ LOC files | Architecture | multiple | Partial v0.4.6 (`tgs_path_auc` extracted) |
| D.2 | 20+ parameter signatures | Architecture | multiple | ✅ Fixed v0.4.7 (additive `tg_priors()` / `tg_mcmc()` config helpers wired into TG / BW / Survival / Synergy / DR) |
| D.3 | ggplot generation inside stat functions | Architecture | multiple | ✅ Fixed v0.4.7 (every Bayesian function already guards plot generation with `if (isTRUE(plots))`; pattern documented in `bayesian_tumor_growth()` @section) |
| E.1 | Add cmdstanr backend option | Enhancement | all `bayesian_*.R` | ✅ Fixed v0.4.6 |
| E.2 | Concordance / C-index | Enhancement | `survival_statistics.R` | ✅ Fixed v0.4.5 |
| E.3 | Posterior growth rates from brms model | Enhancement | `bayesian_tumor_growth.R` | ✅ Fixed v0.4.6 |
| E.4 | Bootstrap CIs for AUC | Enhancement | `tumor_growth_statistics.R` | ✅ Fixed v0.4.6 |
| E.5 | Document effect-size scale | Enhancement | `apriori_power_analysis.R` | ✅ Fixed v0.4.5 |
| E.6 | Posterior predictive coverage | Enhancement | all `bayesian_*.R` | ✅ Fixed v0.4.6 |

---

## G. Deep Audit Findings (added 2026-05-30, files read line-by-line)

The user requested deeper reads of `bayesian_synergy.R`, the dose-response stack, body weight + TWM, and the power-analysis trio.

### G.1 `suppressWarnings(brms::brm(...))` in `bayesian_synergy` — **CRITICAL**
**File:** `R/bayesian_synergy.R` lines 327, 829

Both `bayesian_synergy()` and `bayesian_synergy_over_time()` wrap their `brms::brm()` call in `suppressWarnings()`. Stan emits **warnings** for:
- Divergent transitions (sampling pathology — geometry the sampler cannot integrate accurately)
- Max-treedepth hits (efficiency degradation)
- Rhat above the package's threshold
- Low ESS

All four are silenced. The function then reports `Converged = (Rhat <= 1.01)` with no indication anything went wrong. A model with 500 divergences and `Rhat = 1.00` will pass as "converged". For comparison, `bayesian_tumor_growth()` (line 314) does **not** suppress warnings.

**Fix:** Remove the `suppressWarnings()` wrapper, or replace with `withCallingHandlers()` that captures warnings into a list and exposes them in `mcmc_diagnostics`. Wrapping a Bayesian fit in `suppressWarnings()` is essentially never the right move.

### G.2 Dose-response model selection labels are misleading — **MAJOR**
**File:** `R/dose_response_statistics.R` lines 333–353

```r
dr_model_decr <- drc::drm(... fct = drc::LL.4(...))  # 4-param log-logistic
dr_model_incr <- drc::drm(... fct = drc::LL.5(...))  # 5-param (adds asymmetry)
if (model_aic_decr < model_aic_incr) {
  model_type <- "inhibition"
} else {
  model_type <- "stimulation"
}
statistics$dr_model_type <- model_type
```
The variable names (`decr` / `incr`) and the `model_type` labels (`"inhibition"` / `"stimulation"`) imply directional inference, but **neither LL.4 nor LL.5 determines direction** — both fit either direction via the Slope sign. The selection is between a **symmetric** 4-parameter and an **asymmetric** 5-parameter log-logistic. Reporting `model_type = "stimulation"` to a user who selected the LL.5 because their tumor-growth data happens to be a slightly asymmetric inhibition curve is actively wrong.

**Fix:** Rename to `dr_model_4p` / `dr_model_5p` and `model_type = "symmetric"` / `"asymmetric"`. Derive direction independently from `sign(coef(dr_model)["b:(Intercept)"])` (Slope parameter).

### G.3 Bayesian dose-response is structurally inhibition-only — **MAJOR (or document loudly)**
**File:** `R/bayesian_dose_response.R` line 288

```r
brms_formula <- brms::bf(
  TGI ~ inv_logit(logEmax) / (1 + (exp(logEC50) / Dose)^exp(logHill)),
  logEmax ~ 1, logEC50 ~ 1, logHill ~ 1, nl = TRUE
)
```
- `inv_logit(logEmax)` constrains Emax ∈ [0, 1]
- `exp(logHill)` constrains Hill > 0 → inhibitory direction only

If a dataset shows a stimulatory dose-response, the model cannot represent it — the fit will be poor with no clear signal as to why. Either add a `direction = c("inhibitory", "stimulatory")` argument that inverts the formula, or document the inhibition-only assumption with a prominent warning in `@details` and validate the data direction at call time.

### G.4 `bayesian_synergy` returns fewer diagnostic plots than other Bayesian functions — **MAJOR (consistency)**
**File:** `R/bayesian_synergy.R`

The function returns `synergy_plot` and `posterior_dist_plot` (Bliss / Loewe metric densities), but **not** the four standard plots every other Bayesian function returns:

| Plot | tumor_growth | body_weight | survival | dose_response | **synergy** | TWM |
|---|---|---|---|---|---|---|
| `pp_check_plot` | ✓ | ✓ | ✓ | ✓ | ✗ | ✗ |
| `prior_posterior_plot` | ✓ | ✓ | ✓ | ✓ | ✗ | ✗ |
| `mcmc_trace_plot` | ✓ | ✓ | ✓ | ✓ | ✗ | ✗ |
| `posterior_dist_plot` (treatment params) | ✓ | ✓ | ✓ | ✓ | ✗ (reused for metrics) | ✗ |

The dashboard tabs for these in synergy and TWM end up empty. **Fix:** mirror the `bayesian_tumor_growth()` plot block in both.

### G.5 `bayesian_therapeutic_window` has *no* standard Bayesian diagnostic plots — **MAJOR**
**File:** `R/bayesian_therapeutic_window.R`

A grep for `pp_check|prior_posterior|mcmc_trace|bayes_R2|loo(` returned **zero matches** in the entire file. TWM is the most complex Bayesian model in the package (combines TG and BW posteriors), and it ships with zero standard diagnostic plots. **Fix:** add pp_check, prior_posterior, trace, and bayes_R2.

### G.6 `bayesian_synergy_over_time` duplicates ~300 LOC of `bayesian_synergy` — Architecture
**File:** `R/bayesian_synergy.R`

The two exported functions share data prep, transform handling, cage setup, prior spec, and model fitting almost verbatim. Extract a `bs_fit_synergy_model(df, ..., prior_str, ...)` helper.

### G.7 `analyze_polynomial_trends` assumes equally-spaced doses — **MAJOR (silent wrong answer)**
**File:** `R/dose_response_statistics.R` lines 459–509

```r
analysis_data$dose_factor <- factor(analysis_data[[dose_column]],
                                    levels = sort(unique(analysis_data[[dose_column]])))
stats::contrasts(analysis_data$dose_factor) <- stats::contr.poly(levels(analysis_data$dose_factor))
poly_model <- stats::lm(... ~ dose_factor)
```
`stats::contr.poly()` returns orthogonal polynomial contrasts that are correct **only when factor levels are equally spaced** (it generates contrasts for indices 1, 2, 3, ...). Dose-response experiments routinely use unequally-spaced doses (e.g. 0, 10, 30, 100 mg/kg) — the linear / quadratic / cubic decomposition reported in `linear_trend_pvalue`, `quadratic_trend_pvalue`, `cubic_trend_pvalue` is then **not** the orthogonal decomposition of the dose effect on the actual dose scale.

**Fix:** Use `contr.poly(n, scores = sort(unique(numeric_doses)))` to build contrasts from real dose values, or switch to a continuous-dose polynomial regression: `lm(Volume ~ poly(Dose, 3))`.

### G.8 `analyze_growth_rate` zero-handling diverges from main flow — MINOR (inconsistency)
**File:** `R/dose_response_statistics.R` line 418

`mutate(log_volume = log1p(Volume))` in the growth-rate helper, but the main analysis uses the `log(min/2)` pattern. Pick one across the package.

### G.9 Hill prior centred on `log(median(non_zero_doses))` — design note (good practice)
**File:** `R/bayesian_dose_response.R` line 272

```r
brms::prior_string(
  paste0("normal(", round(log_median_dose, 3), ", ", ec50_sd, ")"),
  nlpar = "logEC50"
)
```
Anchoring the log-EC50 prior at the log of the median observed dose is the right call for a data-aware weakly-informative prior. Worth keeping. Document this in the prior_strength docstring so users know.

### G.10 Power-trio + `analyze_body_weight` — spot-checked, no new findings
**Files:** `R/apriori_power_analysis.R` (243 LOC), `R/apriori_power_simulation.R` (238 LOC), `R/bayesian_power_analysis.R` (292 LOC), `R/analyze_body_weight.R` (202 LOC)

Read function signatures, validation guards, and a few key chunks. All have explicit `stop()` validations for required arguments and document direction conventions (e.g. `treatment_effect >= 0` = inhibition). The Round 1 review covered items 2.1, 2.2, 2.3, 3.5, 3.7 in these files; all remain resolved. No net-new findings to add here at this depth. If you want a true line-by-line audit of these, say so and I'll do another pass.

### G.11 `bayesian_synergy` Loewe denominator floor — verified resolved
**File:** `R/bayesian_synergy.R` line 413

Round 1 item B3.1 flagged the Loewe denominator floor (`1e-6`) as arbitrary and undocumented. At v0.3.6:
```r
fe_max   <- max(fe_combo, na.rm = TRUE)
fe_floor <- max(fe_max * 1e-4, 1e-4)
floor_applied <- any(fe_combo < fe_floor)
```
The floor is now relative to the maximum observed FE_combo (with an absolute lower bound), and `floor_applied` is returned in `loewe_summary` so the dashboard can flag when it kicked in. **Resolved.**

### G.12 Round 2 summary table update

Add to the table in section F:

| ID | Issue | Severity | File |
|---|---|---|---|
| G.1 | `suppressWarnings()` around `brms::brm()` | Critical | `bayesian_synergy.R` | ✅ Fixed v0.4.5 |
| G.2 | DR model labels misleading (`decr`/`incr` = `inhibition`/`stimulation`) | Major | `dose_response_statistics.R` | ✅ Fixed v0.4.5 (symmetric/asymmetric + independent direction) |
| G.3 | Bayesian DR inhibition-only by construction | Major | `bayesian_dose_response.R` | ✅ Fixed v0.4.5 (documented + data-direction warning) |
| G.4 | `bayesian_synergy` missing 4 standard diagnostic plots | Major | `bayesian_synergy.R` | ✅ Fixed v0.4.5 |
| G.5 | `bayesian_therapeutic_window` missing **all** standard diagnostic plots | Major | `bayesian_therapeutic_window.R` | ✅ Fixed v0.4.5 (aliased from input models) |
| G.6 | `bayesian_synergy_over_time` duplicates `bayesian_synergy` | Architecture | `bayesian_synergy.R` | ✅ Fixed v0.4.6 (`bs_fit_synergy_model` helper) |
| G.7 | `analyze_polynomial_trends` assumes equally-spaced doses | Major | `dose_response_statistics.R` | ✅ Fixed v0.4.5 |
| G.8 | `log1p` vs `log(min/2)` inconsistency | Minor | `dose_response_statistics.R` | ✅ Fixed v0.4.5 |
| G.9 | EC50 prior centred on median dose | (Good practice — document) | `bayesian_dose_response.R` | ✅ Fixed v0.4.5 |

---

## H. Corrections to this review

### H.1 B.3 was wrong — retract
The earlier item B.3 claimed `tgs_handle_necrosis()` lives in `R/tumor_growth_statistics.R` and should be moved to `R/utils_necrosis.R`. That's not the case: `tgs_handle_necrosis()` is **already** defined in `R/utils_necrosis.R` (line 29) and is called from both `tumor_growth_statistics.R` and `bayesian_tumor_growth.R`. Treat B.3 as not an issue.

## I. Final-sweep additions

### I.1 Bliss / Loewe formulas duplicated between frequentist and Bayesian synergy
**Files:** `R/analyze_drug_synergy.R` lines ~150–200, `R/bayesian_synergy.R` lines ~392–435

The Bliss-independence and Loewe-CI formulas (`fe_a + fe_b - fe_a * fe_b`, `min(fe_a + fe_b, 1) / fe_combo`, thresholds 0.85 / 1.15 for synergy / additivity / antagonism) appear verbatim in `analyze_drug_synergy()` (scalar form, on point estimates) and `bayesian_synergy()` (vectorised over posterior draws). Extract to `synergy_bliss_excess()` / `synergy_loewe_ci()` helpers — both forms work elementwise on numeric vectors, so a single implementation can serve both call sites.

### I.2 `therapeutic_window_metric()` — `noise_floor = 1.0%` is a magic number — MINOR
**File:** `R/therapeutic_window_metric.R` lines 17, 100–103

The `noise_floor` argument controls the threshold below which `TWM = abs(TGI)` (instead of `TGI / weight_loss%`). Default is 1.0 — i.e. group mean per-mouse max weight loss must exceed 1% before the ratio kicks in. The rationale (avoiding division by near-zero weight loss) is sound, but the choice of 1% is undocumented. Either (a) cite the empirical / clinical basis for 1%, or (b) document that the parameter is exposed and let users tune it for their experiment.

### I.3 `tumor_growth_statistics` `auc_method` arg ignored — A.3 confirmed
The final sweep confirmed: `R/tumor_auc_analysis.R` (the standalone AUC function) properly implements both `"trapezoidal"` and `"last_observation"` (LOCF) at line 202. But `R/tumor_growth_statistics.R` exposes `auc_method` and then calls `calculate_auc()` from `utils_auc.R`, which is trapezoidal-only. The two AUC implementations are split across files and the LME4 path's AUC branch can't reach the LOCF logic. A.3's recommendation (plumb or remove) stands.

---

## J. Deep-read pass 2 — every remaining file (2026-05-30)

The user requested a line-by-line read of files I had only sampled. This section documents net-new findings beyond Section G.

### J.1 `bayesian_therapeutic_window` — draw-pairing assumes independence
**File:** `R/bayesian_therapeutic_window.R` lines 52–54, 238–244

Posterior draws from `bayesian_tumor_growth()` and `bayesian_body_weight()` (independently fitted) are paired draw-by-draw to propagate uncertainty through TWM = TGI / |WL%|. The docstring already documents this as an independence assumption. Worth flagging: if TG and BW models share a residual structure (e.g. an animal that did poorly on both) the marginal posteriors will be correlated, but the paired draws treat them as independent — TWM CrIs will be too narrow or too wide depending on the sign of that correlation. A joint multivariate model (single `brms::brm()` with two responses via `mvbind()`) would resolve this; documenting the limitation in `@details` would be the minimal fix.

### J.2 `bayesian_therapeutic_window` — TWM=1 plot isoline is an approximation
**File:** `R/bayesian_therapeutic_window.R` lines 361–367

```r
ggplot2::geom_abline(slope = -1/100, intercept = noise_floor/100, ...)
```
The exact TWM=1 isoline is piecewise: horizontal at `y = noise_floor/100` for `|x| ≤ noise_floor`, linear `y = |x|/100` outside. The single line `y = (-x + noise_floor)/100` approximates this and is visually fine but not the strict isoline. Minor — document or replace with a piecewise polygon.

### J.3 `plot_combination_index` retains the multiplicative-annotation bug — **MAJOR**
**File:** `R/analyze_drug_synergy_over_time.R` lines 419–423

Round 1 item 3.9 fixed `plot_synergy_trend()` to use additive annotation offsets so labels render correctly when `Time_Point` doesn't start at 0. The sibling `plot_combination_index()` in the same file was missed:
```r
ggplot2::annotate("text", x = max(synergy_summary$Time_Point) * 0.15, ...)
```
If the study starts at e.g. day 7 and ends at day 14, `max * 0.15 = 2.1` — outside the plot area. Apply the same `tp_lo + tp_span * 0.05` pattern as in `plot_synergy_trend()`.

### J.4 `apriori_power_simulation` — `baseline_sd` is a ghost parameter
**File:** `R/apriori_power_simulation.R` lines 36–37, 73

```r
#' @param baseline_sd Numeric > 0. SD of baseline volume in mm³ (unused in
#'   simulation; retained for documentation purposes).
```
The argument is in the signature and stored in `params$baseline_sd` but never used in the simulation (variance is captured by `random_intercept_sd` on the log scale). Either remove from the signature, or use it (vary baseline volumes per mouse instead of using a fixed `baseline_volume`). Ghost parameters are landmines — a future contributor will assume it controls something.

### J.5 `bayesian_power_analysis` — no multi-group support — Missing
**File:** `R/bayesian_power_analysis.R` line 130–139

Hardcoded to 2 groups (`Control` + `Treatment`). The frequentist `apriori_power_simulation` accepts `n_groups` and replicates the treatment effect across non-control groups. The Bayesian version should match — common preclinical designs have ≥3 dose / regimen arms.

### J.6 `bayesian_power_analysis` — no random-slope option
**File:** `R/bayesian_power_analysis.R` line 171

Simulation fits `Log_Vol ~ Treatment * Day + (1 | ID)` — random intercept only. The frequentist sibling and `bayesian_tumor_growth()` both expose `random_effects_specification = c("intercept_only", "slope")`. The mismatch means Bayesian power estimates can't be compared like-for-like with the actual modelling pipeline when users elect random slopes.

### J.7 `bayesian_power_analysis` — no null-distribution / type-I-error check
The "success" criterion is `P(β < -δ | y) > target_prob`. Under the null (`treatment_effect = 0`) this should not exceed `1 - target_prob` on average. The function neither asks users to run a null calibration nor includes one internally. A `null_calibration = TRUE` option that re-runs the simulation with `treatment_effect = 0` would let users verify Bayesian "power" estimates aren't simply prior-driven.

### J.8 `me_result` S3 class is defined but not used by analysis functions
**File:** `R/me_result.R` lines 1–22, 119–168

The roxygen for the class says "Every analysis function returns an `me_result` object with at least: analysis_type / data / results / plots / summary / call / timestamp" — but `tumor_growth_statistics()`, `survival_statistics()`, `bayesian_tumor_growth()`, and all the others return **bare lists**, not `me_result` objects. The only function that uses `new_me_result()` in this file is `repeated_measures_anova()`. The class scaffolding is dead code as far as the main analysis surface is concerned.

**Fix:** either (a) refactor the main analysis functions to wrap their result in `new_me_result()` (and update dashboard accessors), or (b) delete `R/me_result.R` and remove the corresponding NAMESPACE exports. The current state misleads readers about the contract.

### J.9 `tumor_doubling_time` doesn't use composite mouse key — MINOR
**File:** `R/me_result.R` lines 194–261

Iterates `unique(df[[id_column]])` rather than `unique(make_mouse_key(...))`. If the same ID exists in multiple cages or treatments (the scenario the composite key was introduced to handle), they'll be conflated. Inconsistent with the rest of the package.

### J.10 `repeated_measures_anova` zero-handling diverges yet again
**File:** `R/me_result.R` line 311

Uses `log(df[[volume_column]] + 1)` (log1p). This is now the **third** zero-handling pattern in the package:
- `tumor_growth_statistics` and `bayesian_tumor_growth`: `log(x); x[x<=0] <- min_positive/2`
- `analyze_growth_rate`: `log1p(x)` (per backend G.8)
- `repeated_measures_anova`: `log(x + 1)` (this finding)

Pick one and apply across the package.

### J.11 `analyze_body_weight` — `cage_column` accepted but ignored — **MAJOR**
**File:** `R/analyze_body_weight.R` lines 14, 30, 78–82, 116–119

```r
has_cage <- !is.null(cage_column) && cage_column %in% names(df)
if (has_cage) wd$Cage <- as.factor(df[[cage_column]])
...
formula_full <- stats::as.formula(paste(response_col, "~", fixed_terms, "+ (1 + Day | ID)"))
```
`Cage` is added to the data frame but **never enters the formula** in either the full or simplified model. This is the same class of bug as Round 1 critical 1.1 (`handle_cage_effects` silently ignored in `tumor_growth_statistics`). Users who supply a cage column expect cage random effects (or warnings if collinear) — they get neither.

**Fix:** add `(1 | Cage)` to the RE structure when `has_cage`, or document that the argument is data-frame–attached only.

### J.12 `analyze_body_weight` — no diagnostic plots returned
The Bayesian counterpart `bayesian_body_weight()` returns `pp_check_plot`, `prior_posterior_plot`, `mcmc_trace_plot`, etc. The frequentist version returns just model + summary_text. Standard residual / QQ / fitted-vs-residual plots should be added for parity and for any LMM the user is running.

### J.13 `therapeutic_window_metric` — `abs(TGI)` masks negative efficacy
**File:** `R/therapeutic_window_metric.R` lines 102–103

```r
twm$TWM <- ifelse(twm$Mean_Pct_Weight_Loss <= noise_floor,
                  abs(twm$TGI),
                  abs(twm$TGI) / twm$Mean_Pct_Weight_Loss)
```
If a drug **accelerates** tumor growth (`TGI < 0`), the absolute value makes the safety score look positive — a treatment that enhances disease and causes weight loss would still earn a TWM > 0. Should either clamp `TGI < 0` to 0 (no benefit) or document that TWM is only meaningful when `TGI > 0`.

### J.14 `therapeutic_window_metric` — per-mouse baseline / nadir doesn't use composite key
**File:** `R/therapeutic_window_metric.R` lines 78–87

Aggregates `~ ID + Treatment` for baseline and nadir. If the same numeric ID is reused across cages (within the same Treatment), they collapse. The Bayesian TWM uses `make_mouse_key` correctly; the frequentist sibling doesn't.

### J.15 `total_benefit_area` — division-by-zero protection is partial
**File:** `R/total_benefit_area.R` line 66

`merged$TGI[is.nan(merged$TGI)] <- 0` catches `0/0` but `1 - V/0 = -Inf` (and `1 - V/negative = > 1`) isn't handled. When the control mean volume at some early time-point is 0 (unlikely but possible if implants haven't taken), the AUC will be dominated by ±Inf. Add `is.finite` guard.

### J.16 `total_benefit_area` — unit-mismatched AUCs combined with lambda=1
**File:** `R/total_benefit_area.R` line 103

`Benefit = Efficacy_AUC - lambda * Toxicity_AUC` adds two AUCs with different units (TGI% × days vs Weight-loss% × days). At `lambda = 1`, the function treats 1 TGI-point-day as equivalent to 1 weight-loss-point-day — a clinically arbitrary equivalence. Either (a) document the default trade-off and provide guidance on calibrating `lambda` from clinical context, or (b) normalise both AUCs (e.g. by their group means) before combining.

### J.17 `weight_loss_threshold` — no `cox.zph` and no Firth fallback
**File:** `R/weight_loss_threshold.R` lines 136–145

Fits `coxph(Surv(Time, Event) ~ Treatment)` with no PH check and no Firth-correction fallback (which `survival_statistics()` does have). Time-to-weight-loss is exactly the kind of endpoint where PH violations are common (some treatments cause acute loss then recovery), and small-event-count separation is likely. Backend A.1 applies here too, and a Firth fallback path matching `survival_statistics::fit_survival_model()` would harden it.

### J.18 `efficacy_toxicity_bivariate` — toxicity baseline uses ID-only aggregation
**File:** `R/efficacy_toxicity_bivariate.R` line 67

```r
baseline_w <- stats::aggregate(Net_Weight ~ ID, data = wt, FUN = function(x) x[1])
```
Groups by `ID` only, not `make_mouse_key(Treatment, ID)`. Efficacy in the same function (line 92) **does** use the composite key. Same-file inconsistency: the function effectively assumes IDs are globally unique for the toxicity arm and possibly-shared for the efficacy arm.

### J.19 `calculate_volume` — unknown formula silently defaults to ellipsoid
**File:** `R/calculate_volume.R` lines 110–119

```r
Volume <- switch(formula,
  "modified_ellipsoid" = ...,
  "ellipsoid"          = ...,
  ...
  (length_tmp * width_tmp^2 * pi) / 6  # default catch-all
)
```
A typo (`"ellipsoid_3axes"` instead of `"ellipsoid_3axis"`) silently produces the 2D ellipsoid value with no warning. Replace with `match.arg(formula, ...)` to fail loudly.

### J.20 `calculate_dates` — `in_place = TRUE` still uses `assign()` to parent.frame
**File:** `R/calculate_dates.R` lines 180–186

Round 1 item 3.3 was logged as fixed in `calculate_volume.R` — but the sibling `calculate_dates.R` still has the parent-frame assignment pattern:
```r
if (in_place) {
  parent_frame <- parent.frame()
  df_name <- deparse(substitute(df))
  if (exists(df_name, envir = parent_frame)) {
    assign(df_name, result, envir = parent_frame)
  }
}
```
Inconsistent — `calculate_volume.R` was deprecated to be a no-op; `calculate_dates.R` retains the magic. Either deprecate both or apply both. (The `in_place` argument in `calculate_volume.R` is now a no-op-with-warning; users who actually rely on `in_place = TRUE` in `calculate_dates.R` are unknowingly using a different behaviour.)

### J.21 `data.R` — `my_data` is unrelated to the package purpose
**File:** `R/data.R` lines 23–36

```r
#' Example data
#' 
#' A simple example dataset with x and y coordinates.
"my_data"
```
A generic x/y example dataset bundled alongside oncology-specific synthetic data. Looks like leftover template scaffolding. Either repurpose or delete.

### J.22 Plot files — generally well-structured
**Files:** `R/plot_*.R` (7 files, 1,471 LOC total)

Surveyed via grep. All have proper input validation (`stop()` on missing data frames / columns), warnings on missing groups / colors, and fallbacks (alphabetical ordering, color recycling). No statistical bugs found. The only inherited issue is **J.3** (`plot_combination_index` in `analyze_drug_synergy_over_time.R`).

### J.23 Tests — by file count, not coverage
**Files:** `tests/testthat/` (22 test files)

By file count the package has good test breadth — one file per statistical function plus helpers and master-dataset. Worth tracking line-coverage and branch-coverage (e.g. via `covr::package_coverage()`) to surface untested edge cases. Not actually opened in this audit.

---

## J.X Round 2 final consolidated table

Adds to the totals in Section F + Section G. Open items as of v0.3.6 / 2026-05-30:

| ID | Issue | Severity | File |
|---|---|---|---|
| J.1 | TWM independence assumption between TG and BW posteriors | Minor (document) | `bayesian_therapeutic_window.R` | ✅ Fixed v0.4.5 |
| J.2 | TWM=1 plot isoline is approximate | Minor | `bayesian_therapeutic_window.R` | ✅ Fixed v0.4.5 |
| J.3 | `plot_combination_index` multiplicative-offset bug (Round 1 3.9 partial fix) | Major | `analyze_drug_synergy_over_time.R` | ✅ Fixed v0.4.5 |
| J.4 | `baseline_sd` ghost parameter | Minor | `apriori_power_simulation.R` | ✅ Fixed v0.4.5 (now functional) |
| J.5 | Bayesian power: no multi-group support | Missing | `bayesian_power_analysis.R` | ✅ Fixed v0.4.6 |
| J.6 | Bayesian power: no random-slope option | Missing | `bayesian_power_analysis.R` | ✅ Fixed v0.4.6 |
| J.7 | Bayesian power: no null-distribution / type-I check | Minor | `bayesian_power_analysis.R` | ✅ Fixed v0.4.6 |
| J.8 | `me_result` class defined but not used by main funcs | Architecture | `me_result.R` | ✅ Fixed v0.4.6 (docs corrected to honest scope) |
| J.9 | `tumor_doubling_time` no composite key | Minor | `me_result.R` | ✅ Fixed v0.4.5 |
| J.10 | `repeated_measures_anova` uses third zero-handling pattern | Minor (inconsistency) | `me_result.R` | ✅ Fixed v0.4.5 |
| J.11 | `analyze_body_weight` accepts `cage_column` but ignores it in formula | **Major** | `analyze_body_weight.R` | ✅ Fixed v0.4.5 |
| J.12 | `analyze_body_weight` returns no diagnostic plots | Missing | `analyze_body_weight.R` | ✅ Fixed v0.4.5 |
| J.13 | `therapeutic_window_metric` `abs(TGI)` masks negative efficacy | Minor | `therapeutic_window_metric.R` | ✅ Fixed v0.4.5 |
| J.14 | `therapeutic_window_metric` baseline doesn't use composite key | Minor | `therapeutic_window_metric.R` | ✅ Fixed v0.4.5 |
| J.15 | `total_benefit_area` division-by-zero partial guard | Minor | `total_benefit_area.R` | ✅ Fixed v0.4.5 |
| J.16 | `total_benefit_area` unit-mismatch in benefit score | Minor (document) | `total_benefit_area.R` | ✅ Fixed v0.4.5 |
| J.17 | `weight_loss_threshold` no `cox.zph` / Firth fallback | Major | `weight_loss_threshold.R` | ✅ Fixed v0.4.5 |
| J.18 | `efficacy_toxicity_bivariate` toxicity uses ID-only key | Minor | `efficacy_toxicity_bivariate.R` | ✅ Fixed v0.4.5 |
| J.19 | `calculate_volume` silent default on unknown formula | Minor | `calculate_volume.R` | ✅ Fixed v0.4.5 |
| J.20 | `calculate_dates` `in_place` still uses `assign()` to parent.frame | Minor (inconsistency) | `calculate_dates.R` | ✅ Fixed v0.4.5 (deprecated) |
| J.21 | `my_data` example dataset unrelated to package purpose | Minor | `data.R` | ✅ Fixed v0.4.5 (deleted) |

---

## K. Test suite review (2026-05-30)

**Scope:** all 22 files under `tests/testthat/` (4,498 LOC total, 458 `expect_` calls). Read in full for the largest files; sampled where assertions follow the same pattern.

### K.0 Strengths

1. **`helper-fixtures.R` is excellent** — every generator is `set.seed`'d, ground-truth properties are computed in comments alongside the data ("AUC_LOW_GROWTH = 5250", "Bliss expected FE = 0.70 → combo vol = 150"), and multiple scenarios exist for synergy (additive, synergistic, antagonist). Test fixtures of this quality are rare.
2. **`test-master_dataset.R` is exemplary** — verifies the bundled dataset via **directional / ordering assertions** (`expect_gt(tgi("Drug_A High"), tgi("Drug_A Mid"))`), not exact numerics. Robust to numeric drift, captures the scientific intent.
3. **Bayesian tests cache the fit** via a `local({ ... })` closure over `.cached_result`, avoiding redundant brms compiles across `test_that` blocks. Smart.
4. **All Bayesian test files `skip_if_not_installed("brms")`** — graceful degradation for CI without Stan.
5. **22 test files for ~36 source files** — broad coverage by file count.

---

### K.1 **STALE TESTS — Bayesian TG tests reference an API that no longer exists** — Critical
**File:** `tests/testthat/test-bayesian_tumor_growth.R` lines 60, 66

```r
test_that("bayesian_tumor_growth: model_type_used is 'bayes'", {
  ...
  expect_equal(res$model_type_used, "bayes")   # ← src returns "bayes_tg"
})

required_cols <- c("Group", "Adjusted_Mean", "Lower_CL", "Upper_CL")
  # ← src returns Lower_CrI / Upper_CrI per the Round 1 B1.4 harmonisation
```

Confirmed by reading `R/bayesian_tumor_growth.R`: line 578 returns `"bayes_tg"`, lines 373–374 return `Lower_CrI / Upper_CrI`.

These two tests must fail whenever brms is installed. If they don't, the CI either skips Bayesian tests by policy or tolerates failures — either way the test discipline is broken and the Bayesian path is effectively unchecked. **Fix:** update the assertions to match the current API. Then verify the file actually runs in CI (skip-tracking via `testthat::SkipReporter` or similar).

### K.2 Defensive column-name discovery masks shallow assertions — Major
**Files:** `test-tumor_growth_statistics.R` lines 91–103, 183–189; `test-bayesian_*` similar pattern

```r
contrast_col <- intersect(c("contrast", "Contrast", "comparison"), colnames(pw))[1]
est_col      <- intersect(c("estimate", "Estimate"), colnames(pw))[1]
...
if (!is.na(contrast_col) && !is.na(est_col)) {
  # actual direction assertion here
}
```
The test passes whether or not the column exists — if column names change, the `if` body simply isn't entered. The "test" then asserts nothing. This pattern recurs across many files. **Fix:** put `expect_false(is.na(contrast_col))` before the `if`, so missing columns count as a failure.

### K.3 `suppressWarnings(suppressMessages(...))` wraps almost every fit — Major
**Files:** `test-tumor_growth_statistics.R`, all Bayesian tests, several toxicity tests

Every fit is wrapped in `suppressWarnings(suppressMessages(...))`. This silences:
- `lme4` singular-fit and non-convergence warnings (which are exactly what regression tests should be catching)
- brms divergent-transition / max-treedepth / low-ESS warnings (G.1 and A.2 — silenced here too, double blind)
- `emmeans` "interaction is averaged over" warnings (sometimes statistically meaningful)

If the package starts emitting a new convergence warning because of a real regression, the test suite reports green. **Fix:** capture warnings into a vector with `withCallingHandlers()` and assert that no unexpected warnings appear; or scope `suppressWarnings()` to specific known-benign warning classes.

### K.4 No test verifies that direction-changing parameters actually change output — Major
The Round 1 critical 1.1 ("`handle_cage_effects` and `random_effects_specification` silently ignored") was *exactly* this class of bug. The Round 2 finding A.3 (`auc_method` silently ignored) is the same class. **No test in the suite would catch either**.

Specifically missing:
- `auc_method = "trapezoidal"` vs `"last_observation"` — should produce different AUCs on the same data
- `handle_cage_effects = "include_if_not_collinear"` vs `"always_include"` vs `"never_include"` — fixed-effect coefficient table should differ
- `random_effects_specification = "intercept_only"` vs `"slope"` vs `"none"` — variance components should differ
- `polynomial_degree = 1` vs `2` — number of fixed-effect rows should differ
- `extrapolation_points = 0` vs `>0` — number of rows in growth-curve plot data should differ
- `prior_strength = "skeptical"` vs `"diffuse"` (Bayesian) — posterior SDs should be different (skeptical narrower)
- `random_effects_specification = "slope"` in Bayesian — never tested
- `prior_strength = "manual"` in Bayesian — never tested

A small "parameter sensitivity" test per function (different param → different output) prevents the entire "silent ignore" bug class.

### K.5 No tests would have caught Round 2 findings — see audit-paired list

| Round 2 finding | Test that would catch it |
|---|---|
| A.1 / J.17 no `cox.zph` | Assert `is.null(res$ph_test)` is **false** for any survival result |
| A.2 / G.1 no divergence reporting | Assert `"Divergences"` column in `mcmc_diagnostics`; assert `suppressWarnings()` is not used around brm |
| A.3 `auc_method` ignored | Compare `tumor_growth_statistics(auc_method="trapezoidal")` vs `auc_method="last_observation"` — should differ |
| A.4 no LOO/Pareto-k | Assert `res$loo` exists for Bayesian results |
| B.1 extrapolation count always 0 or 1 | Run `extrapolation_points = 5`, assert `> 1` mice extrapolated |
| C.2 no posterior P(effect > 0) | Assert `"P_direction"` in Bayesian pairwise table |
| G.2 DR labels misleading | Generate clear stimulatory data, assert `model_type == "stimulation"` matches the actual slope sign |
| G.3 Bayesian DR inhibition-only | Fit on stimulatory data, assert posterior recovers it (currently can't) |
| G.7 `analyze_polynomial_trends` unequal-dose | Run with doses 0/1/100, assert linear/quadratic decomposition matches `poly(c(0,1,100), 2)` |
| J.11 `analyze_body_weight` ignores `cage_column` | Fit with and without `cage_column`, assert `VarCorr(model)` differs |
| J.13 TWM `abs(TGI)` masks negative efficacy | Pass synthetic data with negative TGI, assert TWM ≤ 0 |
| J.20 `calculate_dates` parent.frame `assign()` | Already untested |

### K.6 Bayesian convergence assertions are too permissive — Minor
All Bayesian tests use `expect_true(all(diag$Rhat < 1.1))`. The current Stan convention is `Rhat ≤ 1.01` (Vehtari et al. 2021). The permissive bound passes models that have material convergence problems. Tighten to 1.01.

Additionally, **no Bayesian test asserts** ESS thresholds (recommended: `Bulk_ESS > 400`, `Tail_ESS > 400`). This means the test suite green-lights brms fits with ESS = 30 — essentially unusable posteriors.

### K.7 Bayesian tests use `n_iter = 500` per chain — Minor (intentional speed compromise)
2 chains × 500 iter = 1000 draws. Documented as a CI-speed choice. Fine for "does it run / does it converge" but means coverage of slow-mixing scenarios is zero. Worth running a slow CI lane occasionally with `n_iter = 2000`.

### K.8 Fixture duplication — Minor
`test-toxicity_functions.R::make_weight_data()` re-implements what `helper-fixtures.R::make_bw_simple()` already provides (with different parameters but same structure). The toxicity local fixture is used in ~10 toxicity tests; using the shared `make_bw_simple()` instead would centralise fixture maintenance.

### K.9 No tests for plot return values — Minor
Plot functions are exercised indirectly (tests pass `plots = FALSE` to avoid generating them). `test-plot_functions.R` exists but I didn't deep-read; per the grep earlier the package's plot files are structurally clean. Worth confirming the file actually invokes the plot generators with realistic input.

### K.10 No coverage measurement — Architecture
There's no `covr::package_coverage()` baseline. A coverage report would surface exactly which branches are untested (e.g., `prior_strength = "manual"`, `firth_correction = FALSE`, `transform = "sqrt"`, `transform = "none"`, the `model_simplified` fallback path in `analyze_body_weight`, the Firth fallback in `survival_statistics`). Running `covr::codecov()` once and committing the coverage badge would make untested code visible.

**✅ Resolved v0.4.7** — `covr (>= 3.6.0)` added to `Suggests`; a
`coverage.R` script lives at the package root for one-line invocation
(`Rscript coverage.R`). The script excludes Bayesian entry points by
default because a full Stan compile dominates the run. Known caveat:
the first invocation surfaces pre-existing test failures
(`test-post_power_analysis.R` references a function deleted in v0.3.4;
`test-toxicity_functions.R` calls `efficacy_metric = "log_cell_kill"`
which is no longer a valid choice) that block `covr` from producing a
clean summary. Both are stale-test issues tracked under K.11 and need
to be cleaned up first — the infrastructure is in place once they are.

### K.11 No regression tests for prior fixes — Minor (process)
The Round 1 review fixed 25+ items. None of the fixes appear to have an accompanying regression test (e.g., a test for "`tumor_auc_analysis` composite ID uses `_` as separator" that would fire if someone reintroduces the old code). This is the most common reason regressions slip back in. Pattern: for every Round 1 / Round 2 fix, add a test labeled with the issue ID.

### K.12 No tests for the `me_result` class — Minor
The class exists but isn't constructed by any analysis function (J.8). The test file `test-utils_and_me_result.R` exists; I sampled it via grep but didn't deep-read. The class itself isn't exercised by the analysis functions, so even if tested in isolation the coverage is misleading.

---

## K Round-2 final test-suite summary

| ID | Issue | Severity |
|---|---|---|
| K.1 | Stale tests reference old `bayes` / `Lower_CL` API | **Critical** (Bayesian path effectively unchecked) | ✅ Fixed v0.4.5 |
| K.2 | Defensive `if(!is.na(col))` masks shallow assertions | Major | ✅ Fixed v0.4.6 (TGS tests tightened) |
| K.3 | Blanket `suppressWarnings(suppressMessages())` hides real signal | Major | ✅ Fixed v0.4.6 (`capture_warnings()` + `expect_no_unexpected_warnings()` helpers) |
| K.4 | No "parameter actually changes output" tests | Major (same bug class as Round 1 1.1) | ✅ Fixed v0.4.6 (`test-param_sensitivity.R`) |
| K.5 | None of the Round 2 findings would have been caught | Major | Partially addressed (K.1, K.6 fixed; K.4 still open) |
| K.6 | Rhat threshold 1.1 too permissive; no ESS assertions | Minor | ✅ Fixed v0.4.5 |
| K.7 | n_iter = 500 is fast but light | Minor | Open |
| K.8 | Toxicity fixture duplicates `make_bw_simple()` | Minor | Open |
| K.9 | Plot return values not tested | Minor | Open |
| K.10 | No code-coverage measurement | Architecture | ✅ Fixed v0.4.7 (covr added to Suggests; `coverage.R` script written; baseline blocked by pre-existing stale tests under K.11 — fix those first to get a clean run) |
| K.11 | No regression tests for prior Round 1 / Round 2 fixes | Process | Open |
| K.12 | `me_result` tested in isolation; not exercised end-to-end | Minor | Open |

---

# Review Round 3 (v0.4.14 — 2026-07-29)

**Reviewer:** Claude Code
**Scope:** Fresh statistical-validity audit of the whole analysis surface, independent of Rounds 1–2. Focus areas chosen to be ones the earlier rounds did *not* cover: multiplicity handling, informative dropout / survivor bias, the experimental unit, ratio-statistic uncertainty, prior scaling, and the internal consistency of the three tumour-growth model paths.

## Round-1 / Round-2 status

Most Round 1–2 items verify as genuinely resolved. **Two do not**, and are re-opened below as R3.1 and R3.2 — both are cases where a fix was written but does not take effect at runtime. Both were confirmed by executing the failing code path, not by reading it.

Findings are graded **Critical** (produces a wrong number that a user will act on), **Major** (statistically unsound or systematically biased), **Minor** (correctness-adjacent, inconsistency, or efficiency).

---

## R3-A. Critical

### R3.1 `p_adjust_method` is silently ignored in the `lme4` and `gam` paths — **Critical**
**Files:** `R/tumor_growth_statistics.R:731, 752, 1113`; `R/tgs_gam.R:299-307`

`p_adjust_method` is `match.arg`'d, documented with `"bonferroni"` as the default and described as "Method for p-value adjustment in pairwise comparisons". It is threaded into `tgs_path_auc()` and used there (`R/tgs_path_auc.R:146`). It is **never referenced again** in the `lme4` or `gam` paths — including `model_type = "lme4"`, which is the function default.

What the `lme4` path actually does:

```r
pairwise_comp <- emmeans::contrast(lsmeans_obj, method = contrasts)   # :1113
```

`emmeans::contrast()` defaults `adjust` to `"tukey"` only for `method = "pairwise"`, `"dunnettx"` for `"trt.vs.ctrl"`, and **`"none"` for a user-supplied list of contrast coefficients** — which is exactly what `contrasts` is here. Verified by execution: a 4-group LMM returns three vs-reference p-values with no adjustment and no "P value adjustment" footnote.

The `gam` path (`tgs_gam_pairwise()`) hardcodes `method = "trt.vs.ctrl", by = time_column`, so it applies Dunnett *within each of the five quantile days* and nothing across them — again independent of `p_adjust_method`.

Consequences:
- The documented default (`bonferroni`) is never applied on the default model path. A user who reads the signature and reports "Bonferroni-adjusted" in a methods section is reporting something the software did not do.
- The three model paths adjust over three different families: `lme4` → none (k−1 contrasts); `gam` → Dunnett within-day, nothing across 5 correlated days; `auc` → `p_adjust_method` over all k(k−1)/2 pairs. The same dataset gives three different multiplicity regimes depending on `model_type`.
- The dashboard does not expose `p_adjust_method` at all (grep across `mouseExperimentDashboard/R/` returns no input binding), so the argument is unreachable from the UI and inert on the default backend path — it is effectively dead except for direct `model_type = "auc"` callers.

This is the third instance of the same bug class (Round 1 §1.1 `handle_cage_effects`; Round 2 §A.3 `auc_method`; Round 2 §J.11 `cage_column`). The parameter-sensitivity tests added in v0.4.6 (`test-param_sensitivity.R`, closing K.4) did not include `p_adjust_method` and so did not catch it.

**Fix:** pass `adjust = p_adjust_method` into the `emmeans::contrast()` calls in both paths. For the GAM path also decide whether the family is per-day or across all (day × contrast) cells and adjust accordingly. Add `p_adjust_method` to `test-param_sensitivity.R`.

**Open question for the maintainer:** what should the canonical comparison family be — vs-reference only (Dunnett), or all pairs (Tukey)? The three paths currently disagree with each other *and* with the dashboard (see dashboard review R3.D1).

---

### R3.2 Round 1 §2.8 is not actually fixed — every pairwise log-rank silently falls back to the omnibus p-value — **Critical**
**File:** `R/survival_statistics.R:487-497`

Round 1 §2.8 replaced the single omnibus log-rank p-value with per-group pairwise tests. The replacement code cannot run:

```r
# in survival_statistics():
surv_obj    <- survival::Surv(df[[time_column]], df[[censor_column]])   # length = nrow(df)
cox_formula <- stats::as.formula(paste("surv_obj ~", treatment_column))
...
# in fit_survival_model(), logrank branch:
pair_data <- df[df[[treatment_column]] %in% c(reference_group, grp), ]
pair_p <- tryCatch({
  pair_diff <- survival::survdiff(cox_formula, data = pair_data)        # :493
  1 - stats::pchisq(pair_diff$chisq, df = 1L)
}, error = function(e) omnibus_p)                                       # :495
```

`surv_obj` is not a column of `pair_data`, so `model.frame()` resolves it from `environment(cox_formula)` — the full-length object — while `Treatment` comes from the subset. Reproduced directly:

```
ERROR (falls back to omnibus p): variable lengths differ (found for 'Treatment')
```

Because the `tryCatch` handler returns `omnibus_p`, the failure is invisible and **the function reproduces the exact behaviour §2.8 was written to eliminate**: one omnibus p-value stamped onto every non-reference row, presented as a per-group comparison.

The sibling Firth-convergence fallback at `:426-432` builds its formula correctly (`Surv(time_col, censor_col) ~ treatment_col` as a string), which is why that path works and this one does not.

**Fix:** build the formula from column names inside the loop, as the Firth branch does, and do not swallow the error into a silent fallback:

```r
pair_formula <- stats::as.formula(
  paste0("survival::Surv(", time_column, ", ", censor_column, ") ~ ", treatment_column)
)
pair_diff <- survival::survdiff(pair_formula, data = pair_data)
```

Add a regression test that asserts the non-reference p-values are **not** all identical when the groups genuinely differ.

---

### R3.3 AUC is integrated over each mouse's own follow-up window — unequal dropout makes the groups non-comparable — **Critical**
**Files:** `R/tumor_growth_statistics.R:481-575` (`tgs_compute_auc`), `R/tgs_path_auc.R:55-126`

`calculate_auc()` integrates from each subject's first to last observed day. `tgs_compute_auc()` records `First_Day` and `Last_Day` per mouse — and then **never uses them**. `tgs_path_auc()` feeds the raw AUCs straight into `aov(AUC ~ Treatment)` and Welch t-tests.

In a preclinical oncology study the follow-up window is not a nuisance — it is an outcome. Control animals reach the IACUC volume limit first and are euthanised; treated animals survive longer. So the control group's AUC is integrated over a *shorter* interval than the treated group's. Because the integrand is a growing quantity, extending the window always adds area. A well-tolerated, genuinely efficacious treatment can therefore accumulate a **larger** AUC than control purely because its animals lived longer, reversing the sign of the estimated effect.

Nothing in the code path detects, corrects, or warns about this. `Extrapolated` and `NumPoints` are carried through to the output table but are not used in the analysis either.

**Fix (pick one and document it):**
1. **Common window (recommended default):** integrate every mouse over `[min_day, T*]` where `T*` is the largest day at which all groups still have data, and report `T*` and the number of animals truncated.
2. **Normalise:** report AUC / (Last_Day − First_Day) (mean volume over follow-up) so the metric is duration-free.
3. **At minimum:** warn when the between-group range of `Last_Day` exceeds some fraction of the study length, and surface per-group median follow-up alongside the AUC table.

Note that the same window is what `extrapolation_points` was presumably meant to paper over — see R3.9 for why that is not an acceptable substitute.

---

### R3.4 `analyze_body_weight(model_type = "gam")` always returns an empty EMM table — **Critical**
**File:** `R/analyze_body_weight.R:120-206`, esp. `:161-168`

v0.4.11 fixed the tumour-growth GAM path by patching the `gamm4` `$gam` stub before handing it to `emmeans` — the stub's class vector lacks `"glm"`/`"lm"` and its `$call` is `NULL`, and without both `emmeans::recover_data.gam` rejects it. That patch lives in `tgs_fit_gamm4_model()` (`R/tgs_gam.R:157-165`).

`analyze_body_weight()` does **not** call `tgs_fit_gamm4_model()`. It fits `gamm4::gamm4()` inline (`:150-159`) and calls `emmeans` on the unpatched stub:

```r
emm <- tryCatch({
  em <- emmeans::emmeans(gam_fit$gam, ~ Treatment, at = list(Day = mean(wd$Day)))
  as.data.frame(em)
}, error = function(e) NULL)          # :167 — swallows the known failure
```

So the body-weight GAM path hits the exact error v0.4.11 diagnosed, the `tryCatch` converts it to `NULL`, and `emmeans_table` comes back empty with no message. The dashboard's Toxicity GAM results panel renders nothing and the user has no way to tell whether that means "no effect" or "the call failed".

Compounding this, the roxygen on `tgs_fit_gamm4_model()` (`R/tgs_gam.R:6-7`) states it is "used by `tumor_growth_statistics()` **and `analyze_body_weight()`**". It is not — `grep` finds exactly one call site, in `tumor_growth_statistics.R:851`.

**Fix:** route `analyze_body_weight()`'s GAM path through `tgs_fit_gamm4_model()` (it already accepts a generic `response_column`, which is why it was written that way), deleting the duplicated inline fit. Correct the roxygen either way. Replace the bare `error = function(e) NULL` with something that records the message.

---

## R3-B. Major — statistical validity

### R3.5 Survivor bias in every endpoint-day TGI — **Major**
**Files:** `R/therapeutic_window_metric.R:77-84`, `R/analyze_drug_synergy.R:128-143`, `R/weight_corrected_tgi.R:83-96`, `R/efficacy_toxicity_bivariate.R:93-114`, `R/analyze_drug_synergy_over_time.R:160`

The shared idiom:

```r
max_day <- max(wd$Day, na.rm = TRUE)
final   <- wd[wd$Day == max_day, ]
ctrl_mean_vol <- mean(final$Volume[final$Treatment == reference_group], na.rm = TRUE)
```

Only animals still on study at the global last day contribute. Since animals leave the study *because their tumours got large*, conditioning on survival to `max_day` selects the slowest-growing animals in every group — and most severely in the control group, which loses animals earliest. `ctrl_mean_vol` is therefore biased downward, and every TGI computed against it is biased **downward** (efficacy understated). In the limit where no control animal reaches `max_day`, `ctrl_mean_vol` is `NaN` and every TGI silently becomes `NaN`.

`R/efficacy_toxicity_bivariate.R:114` has an additional failure mode: `final_vol <- sub$Volume[sub$Day == max_day]` yields `numeric(0)` for any animal that dropped out, which will either error or drop the row depending on how it is consumed downstream.

v0.4.14 fixed precisely this in `bayesian_dose_response()` — "when `endpoint_day` is NULL, take each mouse's last observation (grouped by `id_column`, filter on per-mouse max day)". That fix was not propagated to the five frequentist siblings, and the frequentist functions do not even expose an `endpoint_day` argument.

**Fix:** add an `endpoint_day = NULL` argument to each, with the v0.4.14 semantics (per-mouse last observation when `NULL`), and report per-group *n* at the endpoint alongside every TGI so the reader can see how much of each group survived to contribute.

---

### R3.6 TGI is a ratio statistic and its denominator uncertainty is never propagated — **Major**
**Files:** all TGI consumers, including `R/bayesian_dose_response.R:250-263`

Everywhere in the package, TGI is `1 − V_treated / mean(V_control)`. The denominator is an estimate with its own sampling error, and it is treated as a known constant in every downstream calculation. Two consequences:

1. **Understated uncertainty.** In `bayesian_dose_response()`, every per-animal `TGI` row shares the same plug-in `ctrl_mean` (`:250`), so the observations handed to `brms` are mutually correlated through a source of variation the model does not see. The posterior credible intervals on `Emax`, `EC50`, and `Hill` are consequently too narrow — the Bayesian machinery propagates uncertainty faithfully through everything *except* the one quantity that was silently fixed before the model saw the data.
2. **Arithmetic vs geometric mean inconsistency.** `mean(V_control)` is the arithmetic mean of a right-skewed, approximately log-normal quantity. Every other part of the package models `log(Volume)` — i.e. it assumes the *geometric* mean is the meaningful centre. Using the arithmetic mean in the TGI denominator and the geometric mean in the LMM means the two headline numbers in a report are not describing the same population parameter.

**Fix:** for `bayesian_dose_response()`, model volume directly with the control arm as a group and derive TGI as a posterior transform, so the denominator's uncertainty is inside the model. For the frequentist TGI functions, at minimum add a bootstrap CI (resample animals within group, recompute TGI) and state which mean is used. Pick geometric or arithmetic package-wide and document the choice.

---

### R3.7 `analyze_drug_synergy()` classifies synergy from point estimates with no uncertainty at all — **Major**
**File:** `R/analyze_drug_synergy.R:141-200`

The function reduces each arm to a single group mean, forms three scalar TGIs, and then makes a categorical scientific call from hard-coded thresholds:

```r
if (bliss_difference > 0.1 && loewe_difference > 0.1)  "Strong Synergy"
...
if (ci_value < 0.85)  "Synergistic"   else if (ci_value < 1.15)  "Additive"
```

There is no standard error, no confidence interval, and no test anywhere in the frequentist synergy path. With the group sizes typical of these studies (n = 8–10), the sampling SD of a single TGI is easily 10 percentage points, so `bliss_difference` crossing 0.1 is well within noise. The function will label noise "Strong Synergy" with complete confidence, and the returned `$bliss_independence$synergy` is a bare logical (`bliss_difference > 0`) with no qualification.

Two further defects in the same block:
- `tapply(analysis_data[[volume_column]], analysis_data[[treatment_column]], mean)` (`:132`) uses the default `na.rm = FALSE`; a single missing volume at the evaluation day makes that group's mean `NA` and silently `NA`s every derived quantity.
- `group_means[control_name]` (`:135`) returns `NA` rather than erroring when the named group is absent.

`bayesian_synergy()` does this correctly (draw-wise Bliss/Loewe with credible intervals). The frequentist entry point should either gain a bootstrap CI on `bliss_difference` and `ci_value`, or its documentation should state prominently that it returns descriptive point estimates only and direct inferential use to `bayesian_synergy()`.

---

### R3.8 Bayesian priors are neither per-coefficient nor scaled to the response — **Major**
**Files:** `R/bayesian_tumor_growth.R:364-372`, `R/utils_bayes.R` (`bayes_prior_params`), and the same pattern in BW / survival / synergy

```r
selected_priors <- c(
  brms::prior_string(paste0("normal(0, ", b_sd,       ")"), class = "b"),
  brms::prior_string(paste0("normal(0, ", b_sd * 2.5, ")"), class = "Intercept"),
  ...
)
```

**(a) One prior for coefficients on incommensurable scales.** `class = "b"` applies a single `normal(0, b_sd)` to *every* fixed effect in `Volume ~ Treatment * Day`: the Treatment main effects (differences in log-volume, plausibly ±0.5), the `Day` slope (log-growth per day, ~0.10–0.25), and the `Treatment:Day` interactions (differences in growth rate, ~0.01–0.08). Under the default `prior_strength = "skeptical"` (`b_sd = 0.25`), the prior is meaningfully informative on the Treatment main effect but effectively flat on the interaction — which is the parameter the whole analysis exists to estimate. The prior-strength ladder does not do what its name promises: **"skeptical" is not skeptical about the treatment effect.**

**(b) The Intercept prior ignores the response scale.** With `skeptical`, the Intercept prior is `normal(0, 0.625)`. For the package's own `master_synthetic_data`, volumes span 9.6–4217 mm³, so log-volume spans 2.26–8.35 and the centred intercept sits near 5.5 — about **nine prior SDs** from the prior mean. brms's own default (`student_t(3, median(y), 2.5 * mad(y))`) is data-scaled for exactly this reason. With a few hundred observations the likelihood dominates and the point estimate survives, but: the prior predictive distribution is nonsense (tumours of ~1 mm³), `prior_posterior_plot` — the package's own prior-sensitivity tool — will show extreme prior/data conflict on every real dataset, and on a small pilot the shrinkage is not negligible.

**Fix:** scale the Intercept prior to the data (`normal(median(y), 2.5 * mad(y))` or leave brms's default in place), and set the `class = "b"` priors per coefficient — at minimum split the interaction terms from the main effects via `coef = `. Document the units each prior is expressed in. Re-tune the skeptical/weakly-informative/informative/diffuse ladder against the interaction scale, since that is the estimand.

**Open question:** what volume units do your real uploads use — mm³ (as in `master_synthetic_data`) or cm³ (as in `combo_treatment_synthetic_data`, where volumes are ~0.1–2)? The two differ by a factor of 1000, i.e. ~6.9 on the log scale, which changes how badly (b) bites. If both occur in practice, data-scaled priors are not optional.

---

### R3.9 `tgs_extrapolate()` is unflagged single imputation over informative dropout — **Major**
**File:** `R/tumor_growth_statistics.R:16-117`, called at `:790-793`

When `extrapolation_points > 0`, every subject whose last observation precedes the study maximum gets one synthetic row appended, produced by an OLS fit to its last ≤3 points and then treated by the LMM as an ordinary observation.

Problems, in order of severity:

1. **The dropout is informative and the imputation ignores that.** Animals are missing at the final day mostly because they were euthanised for tumour burden. Extrapolating their trajectory forward and feeding it in as data assumes MAR *within* the fitted model, which is exactly the assumption that fails here.
2. **No uncertainty inflation.** This is single imputation: the predicted value enters with the same weight as a measured one. The LMM's residual variance, standard errors, and p-values are all computed as if the imputed points were observed. Everything downstream is anti-conservative. Multiple imputation or a joint longitudinal-survival model is the statistically defensible route; a single deterministic fill is not.
3. **Linear extrapolation of an exponential process.** `tgs_extrapolate()` runs before the transform at `:817`, so `lm(Volume ~ Day)` is fitted on the **raw** scale for a quantity the rest of the package models as log-linear. The extrapolation is systematically biased low, and `max(0, ...)` at `:68` truncates at zero, creating an artificial floor that then meets the `vol[vol <= 0] <- min(vol[vol > 0])/2` fill.
4. **Covariates are copied from the wrong row.** `new_row <- subject_data[1, ]` (`:71`) clones the subject's **first** observation, so the synthetic final-day row carries day-0 values for every column other than time and volume — including `necrotic_cov_flag` when `necrotic_handling = "covariate"`.
5. **The parameter's value is ignored.** `extrapolation_points` is used only as `if (extrapolation_points > 0)`; the count never reaches `tgs_extrapolate()`, which hardcodes `n_points <- min(3, nrow(subject_data))` and appends exactly one row. `extrapolation_points = 1` and `= 50` do the same thing.
6. **The same name means something else elsewhere.** `plot_tumor_growth(extrapolation_points = "all" | <numeric>)` (`R/plot_tumor_growth.R:52`) uses the name for "how many recent points to fit", which is what a reader would assume it means here too.

**Fix:** at minimum, rename to `extrapolate_to_final_day = FALSE` (logical), fit on the transformed scale, carry the *last* row's covariates, and keep the `Extrapolated` flag in the model frame so imputed points can be excluded or down-weighted. Better: deprecate it in favour of handling the last-day imbalance in the estimand (R3.3 option 1) rather than in the data.

---

### R3.10 `analyze_body_weight()` double-adjusts for tumour burden by default — **Major**
**File:** `R/analyze_body_weight.R:66-73, 104-107`

Defaults are `adjust_tumor_weight = TRUE` and `covariates = c("volume")`. Together they:

```r
wd$Net_Weight <- wd$Weight - wd$Volume / 1000 * tumor_density   # response is now a function of Volume
...
fixed_terms <- paste("Treatment * Day", "+ Volume")             # ...and Volume is a predictor of it
```

The response has had a deterministic linear function of `Volume` subtracted from it, and `Volume` is then entered as a fixed effect predicting that response. The `Volume` coefficient no longer estimates anything interpretable — it absorbs the residual of a quantity already removed by construction — and the treatment effect is adjusted twice for the same confounder, in two different functional forms.

The two adjustments answer different questions and should be mutually exclusive: subtract tumour mass *or* condition on it, not both. As shipped, the default call does both.

**Fix:** make them mutually exclusive (warn and drop `"volume"` from `covariates` when `adjust_tumor_weight` is `TRUE`), and document which question each answers.

---

### R3.11 `analyze_body_weight()` baseline aggregation drops the composite key — **Major**
**File:** `R/analyze_body_weight.R:90-93`

```r
baseline <- stats::aggregate(Net_Weight ~ ID, data = wd, FUN = function(x) x[1])
```

Grouped by `ID` alone. This is the bug fixed in `body_weight_auc.R` (v0.3.1), in `weight_corrected_tgi.R` (Round 1 §1.8), in `therapeutic_window_metric.R` (§J.14) and in `efficacy_toxicity_bivariate.R` (§J.18) — recurring here untouched. When numeric ear-tag IDs are reused across treatment groups (the normal case), all mice sharing an ID collapse to a single `Initial_Mass`, which is then merged back onto every one of them and, if `"initial_mass" %in% covariates`, entered as a model covariate.

Note the neighbouring `weight_corrected_tgi.R:66` does use `~ ID + Treatment`, so the two files disagree within the same package.

Second, subtler issue in the same block: `Initial_Mass` is derived from the mouse's first observation, and that observation is **retained in the response vector**. Regressing a response on a covariate computed from one of its own elements gives that row an exact-fit contribution and biases the covariate's coefficient. The standard remedy is to drop the baseline occasion from the response when it is used as a covariate (constrained-baseline / ANCOVA formulation).

**Fix:** aggregate on `make_mouse_key(Treatment, ID, Cage)` like the rest of the package; drop `Day == min(Day)` rows from the response when `"initial_mass"` is in `covariates`, or document the choice not to.

---

### R3.12 `analyze_body_weight()` returns no pairwise comparisons — **Major (gap)**
**File:** `R/analyze_body_weight.R:260-263, 313`

The result list carries `emmeans_table` (one marginal mean per group) and nothing else inferential — no contrasts, no p-values, no adjustment. Every other analysis function in the package returns a pairwise-comparison table; the dashboard's Toxicity module has a "Pairwise Comparisons" tab. The frequentist body-weight path cannot answer "did this arm lose more weight than control", which is the primary toxicity question.

Also note `emmeans::emmeans(model, ~ Treatment)` at `:261` marginalises over the reference grid's `Day` values rather than an explicit `at = list(Day = ...)`, so it is not the same marginalisation the tumour-growth path uses (`:1052`, mean day). Two functions in one package answering "the group's adjusted mean" two different ways.

**Fix:** add `emmeans::contrast(emm, method = "trt.vs.ctrl", ref = ..., adjust = <p_adjust_method>)` and return it. Add a `p_adjust_method` argument for consistency with `tumor_growth_statistics()`.

---

### R3.13 The frequentist Cox model ignores cage clustering — **Major**
**File:** `R/survival_statistics.R:88-93`, `check_cage_distribution()` at `:281-303`

`cage_column` is accepted, and `check_cage_distribution()` prints a cage × treatment table and a message about collinearity. Cage then **never enters the model** — the Cox fit is `Surv(...) ~ Treatment` with no `cluster(Cage)` or `frailty(Cage)` term.

Cage is a real source of correlation in these studies (shared environment, shared food/water, social hierarchy, and cage is frequently the unit of randomisation). Ignoring it makes the standard errors and p-values anti-conservative when treatment is not perfectly confounded with cage. Where treatment *is* nested in cage, the correct response is to say so loudly — that design cannot separate cage from treatment at all, and the effective sample size is the number of cages, not the number of mice.

`bayesian_survival()` exposes `include_cage_effect`; the frequentist sibling does not. Same asymmetry as J.11.

**Fix:** add `+ cluster(<cage>)` (robust sandwich SEs) or `+ frailty(<cage>)` when cage is not nested within treatment, gated by the same structural check as R3.17; when it *is* nested, emit a warning that the design is confounded rather than a `message()`.

---

### R3.14 `survival_statistics()` requires one row per mouse but neither documents nor validates it — **Major**
**File:** `R/survival_statistics.R:89, 118, 65-66 (roxygen example)`

`Surv(df[[time_column]], df[[censor_column]])` and `survfit(... , data = df)` treat every row of `df` as an independent subject. The function's own event-counting code (`:166-190`) explicitly compensates for "possible duplicates in the data (multiple rows per subject)" — so the code knows the input may be longitudinal, corrects the summary table, and leaves the *model* uncorrected.

If a longitudinal frame is passed, each mouse contributes one row per measurement day: nine censored pseudo-subjects and one event, for a mouse measured ten times. Risk sets, hazard ratios, log-rank statistics, and the KM curve are all wrong, and n is inflated by the number of timepoints.

The dashboard reduces to one row per animal before calling in (`mod_survival.R:640-657`), so UI users are safe. Direct API users are not, and the shipped roxygen example (`:45-56`) passes `combo_treatment_synthetic_data` after `calculate_volume()` + `calculate_dates()` — the full 448-row longitudinal frame, 14 rows per mouse. (It happens to `stop()` first on the missing `Survival_Censor` column, so the example is broken in a second, unrelated way.)

**Fix:** validate at entry — if `id_column` is present and any ID appears more than once, either reduce to the per-mouse last observation or `stop()` with a clear message. Fix the roxygen example to use `master_synthetic_data` (which has `Survival_Censor`) reduced to one row per mouse. State the one-row-per-animal contract in `@param df`.

**Open question:** should the function reduce internally, or refuse and make the caller reduce? Reducing internally would let the dashboard drop its own 18-line reduction loop and remove a divergence point.

---

### R3.15 Power analysis is disconnected from the analysis it is meant to size — **Major**
**File:** `R/apriori_power_analysis.R:96-160`

Two omissions make the returned `Required_N` an underestimate of what the study needs:

1. **No multiplicity.** `alpha` is used as-is. But the analysis these studies actually run is k−1 vs-control contrasts with an adjustment — the package's own documented default is Bonferroni. A 4-arm study powered at `alpha = 0.05` will be analysed at an effective 0.0167 per comparison, and the delivered power is materially below the nominal target. There is no `n_comparisons` or `p_adjust_method` argument to bridge the two.
2. **No attrition.** `Required_N` is the number of *analysable* animals. In this setting animals are lost to euthanasia and to technical failures, and losses are differential by arm (controls first). There is no `dropout_rate` parameter, so the number a user enrols equals the number the calculation assumed would complete.

Separately, when `n_groups >= 3` the calculation answers "power for the omnibus one-way F-test", which is rarely the question — the study is designed to detect specific arm-vs-control differences. The `method_note` added in Round 1 §2.3 documents the d→f conversion caveat but not this framing mismatch.

**Fix:** add `n_comparisons = NULL` (applying the chosen adjustment to `alpha`) and `dropout_rate = 0` (inflating `Required_N` to `ceiling(N / (1 - dropout_rate))`), and report both the analysable and enrolment N. Consider making the k−1 contrast family, not the omnibus F, the default target for `n_groups >= 3`.

---

### R3.16 The AUC path never applies `transform`, but reports that it did — **Major (reporting)**
**File:** `R/tumor_growth_statistics.R:799 vs 817-826`; `R/tgs_path_auc.R:209, 225`

`auc_df <- df` is taken at `:799`, *before* the transform block at `:817`. So AUC is always computed on raw volumes regardless of `transform`. That is arguably the right default — but the summary metadata says otherwise:

```r
volume_transformation  = transform,                                        # tgs_path_auc.R:209
notes = c(if (transform != "none") paste("Volume data was", transform,
          "transformed prior to analysis") ...)                            # :225
```

With the default `transform = "log"`, an AUC run reports "Volume data was log transformed prior to analysis" in the metadata block that feeds the dashboard's methods text and the HTML report export. That statement is false for this path.

**Fix:** in `tgs_path_auc()`, report `volume_transformation = "none (AUC computed on the raw volume scale)"` and drop the note, or plumb the transform through if log-scale AUC is ever wanted.

---

### R3.17 Cage/treatment confounding is decided by a chi-square p-value on observation counts — **Major**
**File:** `R/tumor_growth_statistics.R:185-197, 837-838`

```r
cage_treatment_table <- table(analysis_df[[cage_column]], analysis_df[[treatment_column]])
cage_analysis$collinearity_test <- stats::chisq.test(cage_treatment_table)
...
cage_collinear <- isTRUE(cage_analysis$collinearity_test$p.value < 0.05)
```

`cage_collinear` then decides whether cage enters the model under the default `handle_cage_effects = "include_if_not_collinear"`. Three problems:

1. **Confounding is structural, not a hypothesis.** Whether cage is nested within treatment is a deterministic property of the design: check `all(rowSums(table(cage, treatment) > 0) == 1)`. A significance test is the wrong instrument, and its answer depends on sample size rather than on the design.
2. **The unit is the observation, not the mouse.** Counts are inflated by the number of timepoints, so the chi-square statistic is inflated by roughly that factor and will reject for arbitrarily mild imbalance. `check_cage_distribution()` in `survival_statistics.R:295-301` does the structural check correctly — the two files disagree on how to answer the same question.
3. **Sparse-table warnings.** `chisq.test` on a sparse cage × treatment table will warn about expected counts < 5 on essentially every real study; the warning is neither caught nor surfaced.

**Fix:** replace with the structural nesting check (`survival_statistics.R`'s version), keep the chi-square only as a descriptive extra, and make `handle_cage_effects = "include_if_not_collinear"` branch on the structural result.

---

### R3.18 TWM's noise-floor branch is discontinuous for any `noise_floor != 1.0` — **Major**
**File:** `R/therapeutic_window_metric.R:118-129`

```r
twm$TWM <- ifelse(twm$Mean_Pct_Weight_Loss <= noise_floor,
                  tgi_pos,                                   # units: TGI percentage points
                  tgi_pos / twm$Mean_Pct_Weight_Loss)        # units: dimensionless ratio
twm <- twm[order(-twm$TWM), ]                                # the two are then ranked together
```

The two branches return quantities on different scales, and the rows are sorted into a single ranking. They agree at the boundary **only when `noise_floor == 1.0`**, because dividing by 1 is the identity — the continuity is an accident of the default value, not a property of the formula. `noise_floor` is a user-facing argument documented as tunable ("users with experiment-specific noise estimates should tune accordingly"); setting it to 5 introduces a 5× jump discontinuity at the threshold, and a well-tolerated arm just below the floor will outrank a more effective arm just above it purely from the branch change.

Relatedly, `twm_table` carries no uncertainty of any kind — it is a ratio of two group means with no CI, no n, and no test, and is then sorted to produce what reads as a ranking of treatments. `bayesian_therapeutic_window()` does this properly.

**Fix:** make the transition continuous, e.g. `tgi_pos / pmax(Mean_Pct_Weight_Loss, noise_floor)`, which reduces to the current behaviour at `noise_floor = 1.0` and stays continuous for any value. Add per-group *n* and a bootstrap CI to `twm_table`, or mark the table descriptive-only and point to the Bayesian version for inference.

---

## R3-C. Minor

### R3.19 Growth-rate methods text says `log1p`; the code uses `log()` with a min/2 fill
**Files:** `R/tumor_growth_statistics.R:1173`, `R/tgs_path_auc.R:218` (text) vs `R/tumor_growth_statistics.R:145-146` (code)

Both `analysis_summary$methods$growth_rate_calculation` strings state growth rates are fitted to "log1p-transformed volume data". `tgs_compute_growth_rates()` actually applies `raw_vol[raw_vol <= 0] <- min(raw_vol[raw_vol > 0])/2; log(raw_vol)`. These strings are surfaced verbatim in the dashboard methods panel and the HTML report export, so the discrepancy propagates into anything a user writes up.

### R3.20 The frequentist log transform lacks the all-non-positive guard added to the Bayesian path
**File:** `R/tumor_growth_statistics.R:819-821`

```r
vol[vol <= 0] <- min(vol[vol > 0], na.rm = TRUE) / 2
```

Round 1 §B3.3 added a guard for this exact pattern in `bayesian_tumor_growth` / `bayesian_body_weight` / `bayesian_synergy`, because `min(numeric(0))` is `Inf` and the fill silently becomes `Inf`, then `log(Inf) = Inf`. The frequentist entry point still has the unguarded version, as does `tgs_compute_growth_rates()` at `:145`. Low probability, but it is the same defect the Bayesian side was hardened against, and the fix is three lines.

### R3.21 The `lme4` and `auc` paths adjust over different comparison families
Beyond R3.1: `tgs_path_auc()` builds `utils::combn(treatments, 2)` — all k(k−1)/2 pairs — and adjusts over that family, even though `reference_group` is supplied and the results are then reordered to put reference comparisons first. The `lme4` path builds only the k−1 vs-reference contrasts. For k = 5 that is 10 tests versus 4, so the same `p_adjust_method` means materially different stringency depending on `model_type`.

### R3.22 The AUC path pairs an equal-variance omnibus test with unequal-variance pairwise tests
**File:** `R/tgs_path_auc.R:55-56, 95`

`stats::aov(AUC ~ Treatment)` assumes homoscedasticity; the pairwise tests deliberately use `t.test(..., var.equal = FALSE)` because, per the function's own note, "variances between treatment groups may differ". Both cannot be right. Use `stats::oneway.test(AUC ~ Treatment, var.equal = FALSE)` (Welch ANOVA) for the omnibus so the two agree, and add a Levene / Brown–Forsythe check to the diagnostics (already listed as missing in Round 2 §C).

### R3.23 Dead `unique_id` merge in `tgs_compute_auc()`
**File:** `R/tumor_growth_statistics.R:491-501`

`unique_combinations`, the sequential `unique_id`, and the `merge()` that attaches it are computed and never used — the loop keys on `composite_id` from `make_mouse_key()`. The `merge()` also silently reorders rows and would duplicate rows if `unique_combinations` were ever non-unique. Delete lines 491-495 and index `auc_df` directly.

### R3.24 `tgs_compute_growth_rates()` splits over the full factor cross-product
**File:** `R/tumor_growth_statistics.R:135-137`

`split(auc_df, list(treatment, id, cage))` allocates one list element per *combination* of levels, not per observed combination. For 6 treatments × 60 IDs × 12 cages that is 4,320 elements of which ~60 are non-empty. Pass `drop = TRUE`, or split on `make_mouse_key(...)` directly.

### R3.25 `weight_loss_threshold()` mouse key omits cage, and the function has no cage argument
**File:** `R/weight_loss_threshold.R:57`

`make_mouse_key(wd$Treatment, wd$ID)` — two fields. `therapeutic_window_metric()` uses three (`ID|||Treatment|||Cage`), and `tumor_growth_statistics()` uses three. Same-ID mice in different cages within one treatment arm still collapse here. The function does not accept a `cage_column` at all, so a caller cannot fix it.

### R3.26 Time-to-weight-loss censoring is informative and treated as independent
**File:** `R/weight_loss_threshold.R:96-106`

Animals that never reach the threshold are censored at their last observation. But the dominant reason for a short observation record is euthanasia for tumour burden — which is not independent of the weight-loss hazard. This is a competing-risks setting, and the KM estimate here systematically overestimates weight-loss-free survival. Either move to a cumulative-incidence / Fine–Gray formulation, or document the assumption explicitly and report how many censorings were administrative versus event-driven.

### R3.27 "All events" is treated as separation, contradicting the package's own comment
**Files:** `R/weight_loss_threshold.R:156-157`, `R/survival_statistics.R:325-333`

```r
has_separation <- any(ev_by_grp == 0L) || any(ev_by_grp == n_by_grp)
```

`survival_statistics.R:331` prints "Note: Groups with all events: ... (this is not a problem for Cox models)" and then includes that case in `has_separation` anyway, routing an ordinary dataset to Firth. A group in which every animal has an event is perfectly estimable in a Cox model; only a group with *zero* events causes separation. Drop the `ev_by_grp == n_by_grp` term from both.

### R3.29 No frequentist path returns `transform_used`, so callers cannot tell what scale the results are on
**Files:** `R/tumor_growth_statistics.R:1210-1234` (lme4 return), `:927-955` (gam return), `R/tgs_path_auc.R:250-267`

Every `bayesian_*()` function returns `transform_used`. No frequentist path does — `grep -rn "transform_used" R/` matches only the Bayesian files. `treatment_effects$Adjusted_Mean` from a default `tumor_growth_statistics()` call is a log-volume with nothing in the returned object marking it as such; the only record is the free-text `summary$methods$volume_transformation` string, which is itself wrong for the AUC path (R3.16).

The dashboard works around this by overwriting the field from its own UI input, with an explicit comment saying the backend metadata "is not reliably populated for all paths" (`mouseExperimentDashboard/R/mod_tumor_growth.R:748-750`). That workaround is only available to a caller who already knows what it asked for; a downstream consumer handed a result object cannot recover the scale, which is exactly what it needs in order to back-transform correctly.

This is the concrete cost of Round 2 §B7.1 (schema inconsistency between Bayesian and frequentist returns) being left open.

**Fix:** add `transform_used` (and `model_type_used`, which the `lme4` and `auc` paths also omit while the `gam` path sets it) to all three frequentist return lists.

---

### R3.28 Magic thresholds are undocumented
`0.85` / `1.15` for the Loewe combination-index bands and `±0.1` for the Bliss/Loewe difference bands (`R/analyze_drug_synergy.R:174-200`) are asserted without citation or rationale, and are not exposed as arguments. Round 2 §I.2 made the same observation about `noise_floor`; the resolution there was to document it. Apply the same treatment: cite the source or expose the thresholds.

---

## R3-D. Cross-cutting themes

These are not individual defects so much as patterns the package should decide a policy on.

**D1. The experimental unit is never stated.** Cage appears as a fixed effect, a random effect, a chi-square test, and a printed message, in four functions, with four different treatments — and is absent from the frequentist Cox model entirely (R3.13). If cages are the randomisation unit, the effective n is the number of cages and most of the package is anti-conservative. A single documented policy, applied uniformly, would settle roughly six of the findings above.

**D2. Dropout is treated as a nuisance to be patched, not as a feature of the design.** It surfaces as extrapolation (R3.9), as survivor-biased endpoint TGI (R3.5), as incomparable AUC windows (R3.3), and as informative censoring (R3.26). Each site handles it differently and none handles it correctly. The principled options are a joint longitudinal-survival model, a common analysis window, or an explicitly-stated estimand ("volume among animals surviving to day T"). Any of the three would be defensible; the current mixture is not.

**D3. Derived ratio metrics carry no uncertainty.** TGI, TWM, Bliss difference, and the Combination Index are all reported as bare point estimates by the frequentist path, and several are used to produce categorical scientific calls or rankings. Every one of these has a Bayesian sibling that does it properly. Either add bootstrap CIs across the board or restructure the frequentist functions as descriptive-only with the inferential answer delegated to the Bayesian path.

**D4. Multiplicity has no owner.** Three model paths, the dashboard, and the power module each make their own independent choice, and none of them consults `p_adjust_method` (R3.1, R3.15, dashboard R3.D1). A single `comparison_family` policy — chosen once, threaded everywhere, and reported in the results object — would replace all of it.

**D5. "Silently ignored parameter" is now a five-time recurring bug class.** Round 1 §1.1 (`handle_cage_effects`, `random_effects_specification`, `polynomial_degree`), Round 2 §A.3 (`auc_method`), §J.4 (`baseline_sd`), §J.11 (`cage_column`), and now R3.1 (`p_adjust_method`) and R3.9 item 5 (`extrapolation_points`). The v0.4.6 parameter-sensitivity tests were the right response but cover only the arguments already known to be broken. The generalisable fix is a test that, for each exported analysis function, calls it twice with each documented argument set to two different values and asserts the results differ — driven off the formals list, so new arguments are covered automatically.

---

## R3-E. Status table — Round 3

| ID | Issue | Severity | File | Status |
|---|---|---|---|---|
| R3.1 | `p_adjust_method` ignored in `lme4` + `gam` paths (default path unadjusted) | **Critical** | `tumor_growth_statistics.R:1113`, `tgs_gam.R:299` | ✅ Fixed v0.6.0 |
| R3.2 | Pairwise log-rank always errors → omnibus p reused (Round 1 §2.8 not actually fixed) | **Critical** | `survival_statistics.R:493` | ✅ Fixed v0.5.0 |
| R3.3 | AUC integrated over per-mouse windows; unequal dropout can reverse the effect | **Critical** | `tumor_growth_statistics.R:481`, `tgs_path_auc.R` | ✅ Fixed v0.8.0 |
| R3.4 | BW GAM emmeans always NULL — v0.4.11 gamm4 patch not applied in duplicated path | **Critical** | `analyze_body_weight.R:161` | ✅ Fixed v0.6.0 |
| R3.5 | Survivor bias in every endpoint-day TGI (v0.4.14 fix not propagated) | Major | 5 files | ✅ Fixed v0.8.0 |
| R3.6 | TGI denominator uncertainty never propagated; arithmetic vs geometric mean | Major | all TGI consumers | ✅ Fixed v0.7.0 |
| R3.7 | Frequentist synergy classifies from point estimates, no uncertainty | Major | `analyze_drug_synergy.R:141` | ✅ Fixed v0.7.0 |
| R3.8 | Bayesian priors not per-coefficient and not response-scaled | Major | `bayesian_*.R`, `utils_bayes.R` | ✅ Fixed v0.9.0 |
| R3.9 | `tgs_extrapolate()` = unflagged single imputation over informative dropout | Major | `tumor_growth_statistics.R:16` | ✅ Fixed v0.9.0 (deprecated) |
| R3.10 | BW double-adjusts for tumour burden by default | Major | `analyze_body_weight.R:66,104` | ✅ Fixed v0.5.0 |
| R3.11 | BW baseline aggregated by ID only; baseline used as covariate on itself | Major | `analyze_body_weight.R:91` | ✅ Fixed v0.5.0 |
| R3.12 | BW returns no pairwise comparisons | Major | `analyze_body_weight.R:261` | ✅ Fixed v0.6.0 |
| R3.13 | Frequentist Cox ignores cage clustering | Major | `survival_statistics.R:88` | ✅ Fixed v0.7.0 |
| R3.14 | `survival_statistics()` one-row-per-mouse contract undocumented/unvalidated | Major | `survival_statistics.R:89` | ✅ Fixed v0.5.0 |
| R3.15 | Power analysis ignores multiplicity and attrition | Major | `apriori_power_analysis.R:96` | ✅ Fixed v0.7.0 |
| R3.16 | AUC path never applies `transform` but reports that it did | Major | `tgs_path_auc.R:209` | ✅ Fixed v0.5.0 |
| R3.17 | Cage confounding decided by chi-square p on observation counts | Major | `tumor_growth_statistics.R:191` | ✅ Fixed v0.6.0 |
| R3.18 | TWM noise-floor branch discontinuous unless `noise_floor == 1.0` | Major | `therapeutic_window_metric.R:119` | ✅ Fixed v0.5.0 |
| R3.19 | Methods text says `log1p`, code uses `log()` + min/2 fill | Minor | `tumor_growth_statistics.R:1173` | ✅ Fixed v0.5.0 |
| R3.20 | Frequentist log fill lacks the all-non-positive guard (cf. §B3.3) | Minor | `tumor_growth_statistics.R:820` | ✅ Fixed v0.5.0 |
| R3.21 | `auc` adjusts over all pairs, `lme4` over k−1 — different families | Minor | `tgs_path_auc.R:64` | ✅ Fixed v0.6.0 |
| R3.22 | Equal-variance `aov` omnibus + Welch pairwise in the same path | Minor | `tgs_path_auc.R:55` | ✅ Fixed v0.7.0 |
| R3.23 | Dead `unique_id` merge in `tgs_compute_auc()` | Minor | `tumor_growth_statistics.R:491` | ✅ Fixed v0.5.0 |
| R3.24 | `split()` over full factor cross-product in growth rates | Minor | `tumor_growth_statistics.R:135` | ✅ Fixed v0.7.0 |
| R3.25 | `weight_loss_threshold()` key omits cage; no cage argument | Minor | `weight_loss_threshold.R:57` | ✅ Fixed v0.5.0 |
| R3.26 | Time-to-weight-loss censoring is informative (competing risks) | Minor | `weight_loss_threshold.R:96` | ✅ Fixed v0.7.0 |
| R3.27 | "All events" wrongly treated as separation, contradicting own comment | Minor | `weight_loss_threshold.R:156`, `survival_statistics.R:331` | ✅ Fixed v0.5.0 |
| R3.28 | Undocumented magic thresholds (0.85 / 1.15 / ±0.1) | Minor | `analyze_drug_synergy.R:174` | ✅ Fixed v0.7.0 |
| R3.29 | No frequentist path returns `transform_used` / `model_type_used` | Minor | `tumor_growth_statistics.R`, `tgs_path_auc.R` | ✅ Fixed v0.5.0 |
| R3.30 | mm³ hard-coded in tumour-mass adjustment; cm³ input silently disables it | Major | `analyze_body_weight.R:69`, `weight_loss_threshold.R:50`, `therapeutic_window_metric.R:61` | ✅ Fixed v0.5.0 |
| R3.31 | Jonckheere-Terpstra hardwired to NULL; `clinfun` + 2 reporting branches dead | Major | `dose_response_statistics.R:293` | ✅ Fixed v0.5.0 |
| R3.32 | `bayesian_survival()` errored on every fit without a cage random effect | **Critical** | `bayesian_survival.R` | ✅ Fixed v0.9.0 |
| R3.33 | `test-bayesian_survival.R` used `include_frailty`, renamed in Round 1 B1.5 | Stale test | `test-bayesian_survival.R` | ✅ Fixed v0.9.0 |
| R3.34 | First R3.8 recalibration made slope priors 16 SDs too tight | Regression (caught) | `utils_bayes.R` | ✅ Fixed v0.9.0 |
| R3.35 | `bayesian_synergy()` non-functional since v0.4.6 (`re_term` not found) | **Critical** | `bayesian_synergy.R:710,1080` | ✅ Fixed v0.9.0 |
| R3.36 | BW fixture had zero between-animal variance; RE unidentifiable | Test defect | `helper-fixtures.R` | ✅ Fixed v0.9.0 |
| R3.37 | TWM efficacy-vs-safety plot always NULL (`plot_df` undefined) | Major | `bayesian_therapeutic_window.R:383` | ✅ Fixed v0.9.0 |

---

## R3-F. Open questions for the maintainer

These change what the fix should be, so they are worth answering before the work starts.

1. **What is the canonical comparison family** — vs-reference only (Dunnett) or all pairs (Tukey)? Right now `lme4`, `gam`, `auc`, and the dashboard each answer differently (R3.1, R3.21, dashboard R3.D1).
2. **Is the cage the randomisation unit** in your studies, or are animals randomised individually and merely housed together? This determines whether the package is currently anti-conservative throughout (R3-D1, R3.13, R3.17).
3. **For unequal follow-up (R3.3), which estimand do you want** — a common window truncated at the last day all groups survive to, duration-normalised AUC, or per-mouse last observation? Each is defensible; they answer different questions and give different numbers.
4. **Should `extrapolation_points` survive at all** (R3.9)? Given that the animals it imputes are almost always animals removed for tumour burden, my recommendation is to deprecate it rather than repair it — but if it is in active use for a specific study design I would want to know that design first.
5. **What volume units do real uploads use** — mm³ or cm³? Both appear in the bundled data, and the answer determines how badly the un-scaled Bayesian intercept prior bites (R3.8).
6. **Should `survival_statistics()` reduce longitudinal input to one row per mouse itself**, or refuse it (R3.14)? Reducing internally would let the dashboard delete its own reduction loop and remove a place the two can diverge.
7. **Is the frequentist path meant to be inferential or descriptive?** Several functions (`analyze_drug_synergy`, `therapeutic_window_metric`, `total_benefit_area`) produce categorical scientific calls from bare point estimates while a fully inferential Bayesian sibling exists. Declaring them descriptive-only would be a legitimate and much cheaper resolution than retrofitting bootstrap CIs onto all of them (R3-D3).

---

## R3-G. Maintainer decisions and revised recommendations (2026-07-29)

The §R3-F questions were answered by the maintainer. This section records the answers, the design consequences, and — where an answer changes the analysis — the revised finding. Two findings are **upgraded** on the strength of the answers (R3.17 → Critical, R3.5 → Critical); two are **simplified** to a settled decision (R3.9, R3.1).

---

### G.1 Comparison family — *"Users should be able to define the comparison of interest."*

Settles R3.1, R3.21, and dashboard R3.D1/R3.D3 as a single piece of work rather than a choice between them.

**Design:** add a `comparison_family` argument alongside a *working* `p_adjust_method`, threaded identically through all three model paths:

```r
comparison_family = c("vs_reference", "all_pairs", "custom")
custom_contrasts  = NULL   # named list of coefficient vectors, when "custom"
p_adjust_method   = c("bonferroni", "holm", "fdr", "dunnett", "tukey", "none")
```

Rules that make this coherent:
- The adjustment family must be the set of comparisons actually returned. Today the AUC path adjusts over all k(k−1)/2 pairs while returning a reference-ordered table, and the `lme4` path returns k−1 contrasts adjusted over nothing (R3.21). Deriving the family from `comparison_family` removes the discrepancy by construction.
- `"dunnett"` is only valid with `vs_reference`; `"tukey"` only with `all_pairs`. Validate rather than silently substituting — and note the dashboard currently passes `adjust = "dunnett"` to emmeans, which is not a recognised value and is silently downgraded to the `dunnettx` approximation (dashboard R3.D2). Exact Dunnett is `adjust = "mvt"`.
- The returned object must record what was actually applied (`$comparison_family`, `$p_adjust_method_used`), so the dashboard and the HTML report can state it rather than assert a default.
- The dashboard's existing "Post-hoc Comparisons" dropdown is the natural home for `comparison_family`; add a second control for the correction, and delete the recomputation block at `mod_tumor_growth.R:633-653` so there is exactly one implementation.

**Note on the GAM path:** for `model_type = "gam"` the contrasts are evaluated at five study-day quantiles, so the family is (k−1) × 5 or k(k−1)/2 × 5 cells, not k−1. Decide explicitly whether the day dimension is part of the family. Recommendation: adjust across all returned cells and say so in the table header, since the five days are reported together and read together.

---

### G.2 Cage structure — *"In most experiments, all mice in one cage belong to the same treatment group. Individual mice are treated independently. In some cases the cage may have mice with different treatments."*

**This upgrades R3.17 from Major to Critical**, because the answer reveals that the current default does the opposite of the right thing in the standard design.

The bundled `master_synthetic_data` is exactly the design described — and it is the informative case:

```
        Vehicle  Drug_A Low  Drug_A Mid  Drug_A High  Drug_B  Drug_A Mid + Drug_B
  C01         4           0           0            0       0                    0
  C02         4           0           0            0       0                    0
  C03         0           4           0            0       0                    0
  ...
  cages holding >1 treatment: 0 of 12
  cages per treatment: 2, 2, 2, 2, 2, 2
```

Two cages per treatment, four mice per cage, no mixed cages. There are three structurally distinct cases and they need different handling:

| Case | Structure | Cage as fixed effect | Cage as random effect | Current behaviour |
|---|---|---|---|---|
| **(a) Crossed** | cages contain multiple treatments | Estimable | Estimable | chi-square non-significant → cage included as fixed. Reasonable. |
| **(b) Nested with replication** | ≥2 cages per treatment, one treatment each | **Aliased — not estimable** | **Estimable and required** | chi-square significant → cage dropped entirely. **Wrong.** |
| **(c) Nested, no replication** | exactly 1 cage per treatment | Aliased | Not estimable | chi-square significant → cage dropped, no warning. Should warn loudly. |

Case (b) is the maintainer's stated norm and the shipped dataset's structure. In it, cage cannot enter as a fixed effect — it is perfectly aliased with treatment — but `(1 | Cage)` is fully estimable *because there are two cages per arm*, and fitting it is the statistically correct model. The current default `handle_cage_effects = "include_if_not_collinear"` resolves `cage_collinear = TRUE` and drops cage from the model altogether, discarding a variance component that the design supports estimating.

The consequence is a straightforward understatement of uncertainty. With m = 4 mice per cage, the design effect is `1 + (m − 1) · ICC = 1 + 3·ICC`. At a cage-level ICC of 0.10 the effective n per arm is 8/1.3 ≈ 6.2 rather than 8 and standard errors are understated by roughly 14 %; at ICC = 0.20 it is ~26 %. Cage-level ICCs in that range are ordinary in rodent studies (shared food/water, shared handling, cage position in the rack, social hierarchy). Every p-value and confidence interval from the default path is anti-conservative by that factor.

The maintainer's point that "individual mice are treated independently" is correct and matters — it establishes that the **mouse**, not the cage, is the experimental unit for the intervention, so the design is not pseudoreplicated in the strong sense and cage does not need to be the unit of analysis. But independent *treatment assignment* does not imply independent *outcomes*: co-housed mice share an environment, and that shared environment induces residual correlation which `(1 | Cage)` is exactly the tool for. Both statements are true simultaneously.

**Revised recommendation:**

1. Replace the chi-square test (R3.17) with the structural check already used correctly in `survival_statistics.R:295-301`, and classify into (a) / (b) / (c).
2. Change the default from `include_if_not_collinear` to a structure-aware `"auto"`:
   - case (a) → cage as fixed effect (or random; both are defensible — document the choice)
   - case (b) → `(1 | Cage)` random intercept, **added to the model, not dropped**
   - case (c) → cage omitted, with a `warning()` that cage and treatment are completely confounded and that any cage effect is absorbed into the treatment estimate
3. Make `always_include` error rather than silently produce a rank-deficient design in cases (b) and (c).
4. Apply the same logic to `survival_statistics()` (R3.13 — add `+ cluster(Cage)` or `+ frailty(Cage)` in case (b)) and to `analyze_body_weight()`, whose `re_full` already includes `(1 | Cage)` and is therefore the one function currently getting this right.
5. Report the estimated cage ICC in the diagnostics so users can see how much it mattered.

---

### G.3 Dropout — *"Mice often die earlier than the final time point due to the tumor or being euthanized due to exceeding an IACUC set maximum tumor size or general poor health. The statistical design needs to incorporate this and not discard or truncate."*

This is the right requirement, and the good news is that **the LMM/GAMM path already satisfies it**; what violates it is every derived endpoint-day metric. **R3.5 is upgraded to Critical** on this basis, and R3.3's recommended fix changes.

**Why the mechanism matters.** Euthanasia triggered by a *measured* tumour volume crossing a protocol threshold is missingness that depends only on **observed** data. That is Missing At Random (MAR), not MNAR — and likelihood-based methods (`lmer`, `gamm4`, `brms`) are valid under MAR without modelling the dropout process at all, provided the model conditions on the information that drove the decision (the trajectory up to that point, which it does). Death from general poor health not reflected in measured volume is closer to MNAR and warrants a sensitivity analysis, but the dominant mechanism here is favourable.

This splits the package cleanly into three tiers:

**Tier 1 — already correct; discards nothing.** `tumor_growth_statistics(model_type = "lme4" | "gam")` and `bayesian_tumor_growth()`. Every mouse contributes every observation it has; a mouse euthanised on day 20 still informs its group's intercept and slope. No truncation, no imputation, no exclusion. This is the right foundation and needs no change for dropout.

**Tier 2 — silently broken; used by most of the derived metrics.** Anything computed as a raw mean at a fixed day:

```r
max_day <- max(wd$Day); final <- wd[wd$Day == max_day, ]      # R3.5
ctrl_mean_vol <- mean(final$Volume[final$Treatment == reference_group])
```

These are not likelihood-based, so MAR affords them no protection. They condition on survival to `max_day` and therefore estimate "mean volume **among mice that survived to day T**" — a different quantity, biased hardest against the control arm because control animals hit the IACUC limit first. `therapeutic_window_metric()`, `analyze_drug_synergy()`, `weight_corrected_tgi()`, `efficacy_toxicity_bivariate()`, and the AUC path (R3.3) are all in this tier.

**Tier 3 — treats dropout as the outcome.** Time-to-event endpoints: time to 2× or 4× baseline volume, time to the IACUC limit. Defined for every mouse, uses the death/euthanasia as information rather than as missingness. The package already has this machinery (`survival_statistics()`, and the dashboard's volume-threshold event derivation).

**Revised recommendation — this replaces the three options offered under R3.3:**

> **Derive the endpoint metrics from the fitted model's estimated marginal means at day T, not from raw means of survivors.**

The LMM already estimates each group's mean log-volume at any day T using every observation from every mouse. Computing TGI as `1 − exp(emm_treated(T)) / exp(emm_control(T))` uses the full dataset, discards nobody, truncates nothing, and — unlike appending synthetic rows — propagates the extrapolation uncertainty into the standard error, because it is the model's own uncertainty. This is precisely the "incorporate, don't discard or truncate" behaviour requested, and it is why it also settles G.4.

Supporting changes:
- Report **n at risk at day T** per group alongside every endpoint metric, so a reader can see how much of each arm was still on study. This is the number the current output hides.
- Surface the dropout itself as a first-class result: a Kaplan–Meier or cumulative-incidence curve of time-to-removal per arm, shown next to the efficacy result rather than buried in the Survival tab.
- For AUC (R3.3): compute the model-based integral of the fitted trajectory over `[0, T]` rather than the per-mouse trapezoid over each mouse's own window. Retain the per-mouse trapezoidal AUC as a descriptive column with its `First_Day`/`Last_Day` shown, clearly marked as observed-window-dependent.
- For weight-loss endpoints (R3.26): removal for tumour burden is a **competing risk** for time-to-weight-loss. Report cumulative incidence rather than 1 − KM, and state how many censorings were administrative versus event-driven.
- Where the MNAR concern bites (death from general poor health, not volume-driven), add a documented sensitivity analysis — e.g. refit with a pattern-mixture shift on the dropouts — rather than assuming it away.

---

### G.4 Extrapolation — *"If the proper statistical designs do not need extrapolation, this functionality can be dropped. However, the ability to extrapolate for plotting purposes must remain."*

**Decided.** Under G.3 the statistical path does not need it: the model extrapolates to day T internally and carries the uncertainty, which appending an imputed row never did (R3.9 item 2).

**Action:**
- Remove `extrapolation_points` from `tumor_growth_statistics()` and delete `tgs_extrapolate()` (`R/tumor_growth_statistics.R:16-117`). Deprecate with a warning for one release rather than removing outright, pointing users to the model-based endpoint means.
- **Keep `plot_tumor_growth(extrapolation_points = ...)` unchanged** (`R/plot_tumor_growth.R:52`). It is a separate implementation with different semantics ("all" or the number of recent points to fit) and is genuinely useful for visualising a projected trajectory. The name collision (R3.9 item 6) resolves itself once the statistical one is gone.
- Plots that display extrapolated segments should render them visually distinct — dashed line, or a shaded projection band — so a reader never mistakes a projection for data. Worth checking `plot_tumor_growth()` already does this.
- Drop `Extrapolated` handling from `tgs_compute_auc()` (`:521-542`) once the source of extrapolated rows is gone.

---

### G.5 Volume units — *"Users can upload data that is either in cubic mm or cm, but mm is most common."*

**Confirms R3.8 as a real problem in the common case.** For mm³, log-volume runs roughly 2–9 (`master_synthetic_data`: 2.26–8.35). The default `prior_strength = "skeptical"` sets the Intercept prior to `normal(0, 0.625)` — about nine prior SDs from where the data actually sit. For cm³ the same prior is roughly reasonable. A fixed prior cannot serve both, and it is mis-specified for the more common of the two.

**Action:**
1. **Scale the Intercept prior to the response**, as brms does by default: `normal(median(y), 2.5 * mad(y))` computed on the modelling scale after transformation. This is unit-agnostic and eliminates the problem for both inputs.
2. **Set `class = "b"` priors per coefficient**, splitting the `Treatment:Day` interaction terms from the main effects (R3.8a). The interaction is the estimand and its plausible scale (~0.01–0.08 per day on the log scale) is unit-*independent* — a growth-rate difference does not change when volumes are rescaled — so a fixed prior is defensible there in a way it is not for the intercept. Re-tune the skeptical/weakly-informative/informative/diffuse ladder against that scale.
3. **Detect and record the unit.** A median volume below ~50 on input is almost certainly cm³; above, mm³. Surface the inference in the result object and the dashboard, and let the user override. This is worth doing independently of priors, because TGI, AUC units (mm³·day vs cm³·day), and the `Net_Weight` tumour-mass adjustment in `analyze_body_weight()` (`wd$Volume / 1000 * tumor_density`, `:69`) **all assume mm³**. A cm³ upload silently under-subtracts tumour mass by a factor of 1000, making the adjustment a no-op and the toxicity result wrong with no warning. That is a net-new finding from this answer — logged as **R3.30** below.

---

### R3.30 The tumour-mass weight adjustment hard-codes mm³ and fails silently on cm³ input — **Major (new, from G.5)**
**Files:** `R/analyze_body_weight.R:69`, `R/weight_loss_threshold.R:50`, `R/therapeutic_window_metric.R:61`

```r
wd$Tumor_Weight <- wd$Volume / 1000 * tumor_density   # mm³ → cm³ → g
```

Correct for mm³. For a cm³ upload the conversion should be `Volume * tumor_density` — the `/ 1000` makes the subtracted mass 1000× too small, so `Net_Weight` is effectively unadjusted body weight while the result is labelled "Net Weight (body − tumor)". Since the whole point of the adjustment is to stop a large tumour masking treatment-induced weight loss, a silent no-op inverts the conclusion for exactly the animals it matters most for.

The same hard-coded `/ 1000` appears in all three toxicity functions. None validates the input scale.

**Fix:** add a `volume_units = c("mm3", "cm3")` argument (or take it from the unit detection in G.5 item 3), and range-check: warn when the implied tumour mass exceeds a plausible fraction of body weight, which catches the error from either direction.

---

### G.6 Frequentist inference vs description — answered separately

The maintainer asked for more information before deciding. See the analysis returned in-session; the choice and its consequences are:

- **Option A (make them inferential):** add a nonparametric bootstrap to `analyze_drug_synergy()` and `therapeutic_window_metric()`. Resample **mice within group** (not observations), including the control arm, and recompute the entire statistic per resample. This produces a percentile CI *and* propagates the TGI denominator uncertainty (R3.6) in the same pass, since the control mean is recomputed on each resample. ~40 lines per function, no new dependency, and `tgs_boot_diff_ci()` (`R/tumor_growth_statistics.R:457`) is already the pattern to follow.
- **Option B (make them descriptive):** keep the point estimates, remove or clearly demote the categorical verdicts, and document that inference comes from `bayesian_synergy()` / `bayesian_therapeutic_window()`.

Recommendation: **A for `analyze_drug_synergy()` and `therapeutic_window_metric()`** — both produce headline results that drive decisions and figures, and both currently emit a categorical verdict from an unqualified point estimate. **B for `total_benefit_area()`, `weight_corrected_tgi()`, and `efficacy_toxicity_bivariate()`**, which are exploratory summaries. Pending confirmation.

---

### G.7 Revised severities and settled items

| ID | Change | Reason |
|---|---|---|
| R3.17 | Major → **Critical** | Default cage handling discards estimable cage variance in the maintainer's stated standard design (G.2 case b); SEs understated ~14–26 % at plausible ICC |
| R3.5 | Major → **Critical** | Confirmed against the stated dropout mechanism; endpoint-day means estimate a survivor-conditional quantity, biased hardest against control (G.3) |
| R3.3 | Fix revised | Model-based endpoint means replace the three truncation/normalisation options originally offered (G.3) |
| R3.9 | **Settled — deprecate** | Statistical extrapolation removed; `plot_tumor_growth()` extrapolation retained (G.4) |
| R3.1 | **Settled — design agreed** | `comparison_family` + working `p_adjust_method`, one implementation shared with the dashboard (G.1) |
| R3.8 | **Settled — data-scale priors** | mm³ confirmed as the common case; fixed intercept prior cannot serve both units (G.5) |
| R3.13 | Scope extended | Cox needs `cluster()`/`frailty()` in case (b), gated on the same structural check as R3.17 (G.2) |
| R3.30 | **New — Major** | mm³ hard-coded in the tumour-mass adjustment; cm³ input silently disables it (G.5) |
| R3.6 | Fix unified | Bootstrap in G.6 Option A resolves denominator uncertainty in the same pass |

---

## R3-H. Resampling strategy — where bootstrap, where permutation (2026-07-29)

**Decision recorded:** G.6 Option A is confirmed. `analyze_drug_synergy()` and `therapeutic_window_metric()` gain a mouse-level nonparametric bootstrap; `total_benefit_area()`, `weight_corrected_tgi()`, and `efficacy_toxicity_bivariate()` are documented as descriptive.

The maintainer then asked whether permutation testing is worth adding elsewhere. It is, in four specific places and not in the others — the two techniques answer different questions and are not interchangeable. The governing distinction:

> **A bootstrap estimates the sampling distribution of an estimator** — use it when you want an interval around a quantity (TGI, TWM, Bliss excess, a mean difference).
> **A permutation test generates a null distribution by re-randomising labels** — use it when the null hypothesis is genuinely "the labels carry no information", and only when the thing you permute is exchangeable under that null.

Applying that test to each candidate site:

### H.1 Dose-response trend — the Jonckheere-Terpstra test is already built and switched off — **Major (new finding, R3.31)**

The canonical permutation-based test for an ordered alternative across dose groups is Jonckheere–Terpstra, and this package already has the whole apparatus for it:

- `clinfun` is declared in `DESCRIPTION` `Suggests` **solely** for this test — nothing else in `R/` references the package.
- `dose_response_statistics()` returns `trend_test$jonckheere_test` (`R/dose_response_statistics.R:90`).
- Two reporting branches consume a non-NULL result (`:629-635`, `:675-679`).
- `tests/testthat/test-dose_response_statistics.R:126-136` calls `clinfun::jonckheere.test()` **directly** and asserts significance on the fixture.

And then:

```r
# Jonckheere-Terpstra test not run due to mentioned issues
jt_result <- NULL                                          # :293-294
```

There are no "mentioned issues" anywhere in the file — the comment refers to something that does not exist in the source. So `jonckheere_test` is permanently `NULL`, both reporting branches are unreachable, `clinfun` is a dependency the package never loads, and the test suite verifies a test the package does not run. That last point is the worst of it: the suite is green on a feature that is disabled, which is exactly the false-assurance pattern flagged in Round 2 §K.2.

This is also the cleanest available fix for Round 2 §G.7 (unequal dose spacing breaking `contr.poly`). Jonckheere–Terpstra uses only the **ordering** of the dose levels, never their spacing, so doses of 0/10/30/100 mg/kg pose it no difficulty at all. It is strictly more robust than the orthogonal-polynomial decomposition for the question "does response change monotonically with dose", and it makes no normality assumption.

**Action:** re-enable it (`clinfun::jonckheere.test(volume, dose, alternative = ...)` guarded by `requireNamespace`), pick the alternative from the observed direction rather than hardcoding, document the field in the `@return` block, and note in `@details` that it is a permutation test insensitive to dose spacing. If it was disabled for a real reason, that reason needs to be written down — as it stands the comment is unactionable.

### H.2 AUC pairwise comparisons — permutation p-value to pair with the existing bootstrap CI — **worth it**

`tgs_path_auc()` uses Welch's t-test on per-mouse AUCs, with an optional percentile bootstrap CI (`auc_bootstrap_n`, `tgs_boot_diff_ci()` at `R/tumor_growth_statistics.R:457`). At n = 8–10 per group, on a derived statistic whose distribution is right-skewed, Welch's t is precisely the approximation a permutation test improves on — and the null ("group label carries no information about a mouse's AUC") is exactly a label-exchangeability null, so permutation is licensed here.

The natural pairing is **bootstrap for the interval, permutation for the p-value**; they are not redundant. Roughly 20 lines reusing the existing helper's structure, no new dependency. Worth gating behind the same `auc_bootstrap_n`-style argument so it stays optional.

### H.3 Global test of a treatment effect on trajectory — **valuable, but blocked on one question**

Permuting treatment labels across whole mouse trajectories and refitting gives an exact test of "treatment has no effect on the growth trajectory", free of two assumptions the current path leans on: the Kenward–Roger / Satterthwaite denominator-df approximation (known to be anti-conservative at these group sizes) and normality of the random effects.

**The blocker is what to permute, and it must mirror how randomisation was actually done.** Given G.2 (one treatment per cage, ≥2 cages per arm), there are two possibilities with different consequences:

- **Mice individually randomised to treatment, then re-housed by group.** Permutation unit is the mouse; 48 mice give ample resolution. But cage then becomes a *post-randomisation* variable, and permuting mice across cages breaks the cage structure that G.2 establishes is real.
- **Cages assigned to treatments.** Permutation unit is the cage. With 12 cages across 6 arms there are 12!/(2!)^6 ≈ 7.5 million distinct assignments, so enumeration is not the constraint — but the *effective* sample size for the test is 12, not 48, and power drops accordingly.

That second number is uncomfortable and it is also honest: if cages were the unit of assignment, 12 is the real replication and the current analysis is overstating precision regardless of whether a permutation test is run. The permutation framing makes that unavoidable rather than optional, which is an argument for doing it.

**Open question:** *In your studies, are individual mice randomised to treatment and then re-housed by group, or are cages assigned to treatments?* This determines the permutation unit here and also refines the G.2 cage recommendation.

### H.4 Survival with separation — exact permutation log-rank — **worth it, and it routes around R3.2**

When a group has zero events, the Cox partial likelihood is not estimable, and the package falls back to Firth (`coxphf`) or an asymptotic log-rank whose chi-square approximation is poor at low event counts. A permutation log-rank (`coin::logrank_test(..., distribution = "exact")`) is valid in exactly that regime: it makes no large-sample appeal and needs no bias correction, because it is not estimating a hazard ratio at all.

Note the division of labour — Firth corrects the **estimate's** small-sample bias, permutation gives a valid **p-value**. They are complementary, and reporting a Firth hazard ratio alongside a permutation p-value is a defensible pairing in the separation case.

This would also make R3.2's broken pairwise fallback moot: rather than repairing a `survdiff()` call that has never executed, replace the branch with pairwise permutation log-rank tests, which are correct by construction in the regime that triggers the branch. Cost: a new `coin` dependency (well-maintained, in Suggests, gated by `requireNamespace`).

### H.5 Synergy — permutation is the **wrong** tool here

This is the case most likely to be reached for by analogy, and it does not work. The Bliss / Loewe null is not a label-exchangeability null; it is a **structural** null — "the combination arm's effect equals what Bliss (or Loewe) predicts from the two monotherapy arms". There is no permutation of treatment labels that generates a null distribution for that hypothesis, because permuting the labels destroys the monotherapy-vs-combination structure in terms of which the null is defined. Shuffling `Drug A` / `Drug B` / `A + B` / `Control` labels produces a world in which "the combination arm" is no longer the combination arm.

The bootstrap in G.6 Option A is the right and sufficient tool: it delivers a credible interval on the Bliss/Loewe excess and, because the control arm is resampled too, propagates the TGI denominator uncertainty (R3.6) in the same pass. If a formal test of the additivity null is wanted later, it needs a model with an explicit interaction term (which is what `bayesian_synergy()` already provides), not resampling.

### H.6 TWM, TGI, total benefit area — bootstrap only

These are **estimands**, not hypotheses. The question is "how large is it, and how sure are we", not "could the labels have produced this by chance". Bootstrap intervals; no permutation.

### H.7 Not a substitute for the mixed model

A permutation test returns a p-value and nothing else — no estimate, no interval, no variance components, no marginal means, and no way to compute the model-based endpoint means that G.3 depends on. It is a complement to the LMM path, never a replacement for it.

### H.8 Caveat — permutation interacts with informative dropout

Permuting treatment labels permutes each mouse's dropout pattern along with its trajectory. Under the **strong global null** (treatment affects nothing, including time on study) that is correct and the test is exact. But bolting a permutation test onto a survivor-conditional endpoint statistic — a mean at `max_day` among survivors — makes the null distribution inherit the same survivor bias as the observed statistic, so the test comes out valid for a null nobody intended. **Sequence matters: fix the estimand (G.3) first, then permute.** Do not use a permutation p-value to paper over Tier 2 metrics.

### H.9 Priority and cost

| # | Site | Technique | Cost | Priority |
|---|---|---|---|---|
| H.1 | Dose-response trend (R3.31) | Permutation (Jonckheere–Terpstra) | ~free — already built, one line disabled | **1st** |
| H.2 | AUC pairwise | Permutation p-value + existing bootstrap CI | ~20 LOC, no new dependency | 2nd |
| — | Synergy / TWM (G.6 Option A) | Bootstrap | ~40 LOC each | 2nd |
| H.4 | Survival under separation | Exact permutation log-rank | new `coin` Suggests; also retires R3.2 | 3rd |
| H.3 | Global trajectory null | Permutation of treatment labels | moderate; **blocked on randomisation unit** | 4th |

---

### R3.31 Jonckheere-Terpstra test is hardwired to NULL; `clinfun` dependency and two reporting branches are dead — **Major (new)**
**File:** `R/dose_response_statistics.R:293-294`, with consumers at `:90`, `:629-635`, `:675-679`

See H.1 above for the full analysis. Summary: the test is disabled by `jt_result <- NULL` under a comment referencing "mentioned issues" that appear nowhere in the file; `clinfun` is a `Suggests` dependency used only by the test suite; the returned `trend_test$jonckheere_test` field is always `NULL` and undocumented in the `.Rd`; two reporting branches are unreachable; and `tests/testthat/test-dose_response_statistics.R:126-136` verifies the test independently, so the suite reports green on functionality the package does not execute.


---

## R3-I. Implementation log — v0.5.0 (2026-07-29)

First batch of Round 3 fixes landed. Selection criterion: findings that are
self-contained, produce a demonstrably wrong number today, and are unblocked by
the §R3-F/§R3-G decisions.

**Closed:** R3.2, R3.10, R3.11, R3.14, R3.16, R3.18, R3.19, R3.20, R3.23,
R3.25, R3.27, R3.29, R3.30, R3.31. **Partial:** R3.1 (survival path done; the
lme4 / gam `comparison_family` work is the larger G.1 change and is still open).

Each fix has a regression test in `tests/testthat/test-code_review_round3.R`
labelled with its finding ID — closing the §K.11 gap for Round 3 at least.
R3.1 and R3.2 assert observable behaviour rather than the presence of code,
because both were cases where a fix existed in the source and never executed.

Test suite: **359 passing / 0 failing**, up from 352 passing / 7 failing. All
seven baseline failures were §K.11 stale tests referencing APIs removed in
v0.3.4 and earlier; removing or repairing them unblocks the `covr` baseline
that §K.10 flagged as blocked.

### Still open, in suggested order

| Priority | Finding | Why it is next | Blocked on |
|---|---|---|---|
| 1 | **R3.1** (lme4/gam) + G.1 | Default model path still reports unadjusted p-values as Bonferroni | — |
| 2 | **R3.17** + R3.13 (cage) | Default discards estimable cage variance; SEs understated 14–26 % | — |
| 3 | **R3.5** + R3.3 (endpoint estimand) | Survivor-biased TGI across 5 functions; needs model-based EMMs | — |
| 4 | **R3.4** (BW GAM) | Marginal-means table silently empty | — |
| 5 | **G.6 Option A** + R3.6 | Bootstrap CIs for synergy / TWM; fixes denominator uncertainty in the same pass | — |
| 6 | **R3.9** / G.4 | Deprecate statistical extrapolation, keep the plotting one | depends on 3 |
| 7 | **R3.8** / G.5 | Data-scaled Bayesian priors; per-coefficient `class = "b"` | — |
| 8 | **H.2, H.4** | Permutation p-value for AUC; permutation log-rank (retires R3.2's branch) | — |
| 9 | **H.3** | Global trajectory permutation test | **randomisation unit question** |
| — | R3.12, R3.15, R3.21, R3.22, R3.24, R3.26, R3.28 | Smaller items | — |

Dashboard items (R3.D1–R3.D5) are mostly blocked on backend G.1; R3.D5 (sort by
day before taking the per-animal last row) is unblocked and small.


---

## R3-J. Implementation log — v0.6.0 (2026-07-29)

Second batch. Two deliberate behaviour changes here, both of which alter numbers
users may already have published — flagged as BREAKING in the changelog.

**Closed:** R3.1 (all paths), R3.4, R3.12, R3.17, R3.21.
**Partial:** R3.13 — `classify_cage_structure()` now exists and gives the Cox
path everything it needs, but the `cluster()` / `frailty()` term is not yet
added to `survival_statistics()`.

### G.1 as built

`comparison_family` (`"vs_reference"` / `"all_pairs"` / `"custom"`) plus a
working `p_adjust_method` (now including `"dunnett"` → exact `mvt`, and
`"tukey"`), resolved once by `resolve_comparison_spec()` and shared by all three
tumour-growth paths and `analyze_body_weight()`. The adjustment covers exactly
the returned family, so the two can no longer diverge. Invalid pairings are
rejected: Dunnett requires `vs_reference`, Tukey requires `all_pairs`, and both
are refused on the AUC path, whose independent Welch t-tests have no joint
covariance to exploit.

For the GAM path the family is every returned (contrast × day) cell rather than
each day in isolation — the five quantile days are reported and read together.
Tukey/Dunnett cannot be re-derived across strata without the full joint
covariance, so those stay within-day and are labelled as such in
`Adjust_Scope`. Verified: Bonferroni over 15 GAM cells, 6 all-pairs, and 3
vs-reference contrasts all differ by exactly the expected factors.

### G.2 as built

`classify_cage_structure()` operates on **mice**, not observations — the old
chi-square counted rows, so its statistic was inflated by the number of
timepoints. The new `handle_cage_effects = "auto"` default maps structure to
placement: crossed → fixed, nested-with-replication → `(1|cage)`,
completely-confounded → omitted with a warning.

Measured impact on a 2-cages-per-arm / 4-mice-per-cage design with a real cage
effect (ICC 0.66): treatment-contrast SEs went from 0.141 under the old default
to 0.305 with cage restored — a factor of 2.2. At the more typical ICC of
0.1–0.2 the factor is ~1.15–1.25. Either way the old default was
anti-conservative, and it was anti-conservative *specifically* in the design the
maintainer identified as standard.

`cage_analysis$icc` now reports the cage ICC so this is visible rather than
inferred; the design effect on effective sample size is 1 + (m − 1) × ICC.

### Still open, in suggested order

| Priority | Finding | Why it is next | Blocked on |
|---|---|---|---|
| 1 | **R3.5** + R3.3 | Survivor-biased endpoint TGI in 5 functions; needs model-based EMMs at day T (G.3) | — |
| 2 | **G.6 Option A** + R3.6 | Bootstrap CIs for synergy / TWM; fixes denominator uncertainty in the same pass | — |
| 3 | **R3.13** (finish) | Add `cluster()` / `frailty()` to the Cox path, gated on the structure check | — |
| 4 | **R3.9** / G.4 | Deprecate statistical extrapolation, keep the plotting one | depends on 1 |
| 5 | **R3.8** / G.5 | Data-scaled Bayesian priors; per-coefficient `class = "b"` | — |
| 6 | **H.2, H.4** | Permutation p-value for AUC; permutation log-rank | — |
| 7 | **R3.15** | Power analysis: multiplicity + attrition | — |
| 8 | **H.3** | Global trajectory permutation test | **randomisation unit question** |
| — | R3.22, R3.24, R3.26, R3.28 | Smaller items | — |

### Dashboard

R3.D1–R3.D4 are now **unblocked** — the backend applies the adjustment and
records `comparison_family` / `p_adjust_method_used`, so the dashboard should
delete its recomputation block (`mod_tumor_growth.R:633-653`) and render
`result$pairwise_comparisons` directly. That also resolves R3.D4 (the GAM
"All Comparisons" tab) as a side effect, since the failing `emmeans()` call on
the raw `gamm4` object goes away with the block. R3.D5 remains unblocked and
small. New dashboard work: surface `cage_analysis$icc`, and relabel the
cage-handling selector for the `"auto"` default.


---

## R3-L. Findings surfaced by installing brms (2026-07-29)

Installing `brms` 2.23 and `coin` 1.4.5 locally made ~150 previously-skipped
assertions execute for the first time. **Six findings fell out immediately, two of
them Critical**, all invisible because the Bayesian test files skip wholesale when
brms is absent — a live and expensive illustration of §K.3 and §K.11.

The headline: **`bayesian_survival()` errored on every call without a cage random
effect, and `bayesian_synergy()` had been entirely non-functional since v0.4.6** —
five releases, with the dashboard advertising both. Neither would have been found
without installing brms.

**The single highest-value process change available to this project is making the
Bayesian tests run in CI.** A suite that skips its heaviest, most complex path by
default is not testing that path, and §K.3/§K.11 predicted exactly this. Adding
brms and coin to a CI lane — even a slow nightly one — is worth more than any
further review round.

### R3.32 `bayesian_survival()` errored on every fit without a cage random effect — **Critical (pre-existing)**
**File:** `R/bayesian_survival.R` (both prior branches)

The function declared a `class = "sd"` prior unconditionally, but the formula only
contains a random effect when `include_cage_effect = TRUE` and a cage column is
present. brms rejects a prior that matches no model parameter:

```
Error in .validate_prior(...): The following priors do not correspond to any
model parameter: sd ~ exponential(1)
```

So **every** `bayesian_survival()` call with `include_cage_effect = FALSE`, or with
no cage column, failed outright. Not a degradation — a hard error on a documented
code path. The `sd` prior is now conditional on `use_cage_re`.

This had been shipped since the function was written. It was invisible because
`test-bayesian_survival.R` skips without brms, and brms was not installed in the
development environment.

### R3.33 `test-bayesian_survival.R` used an argument renamed in Round 1 — **Stale test**
Round 1 §B1.5 renamed `include_frailty` to `include_cage_effect`. The test file
still passed `include_frailty`, producing `unused argument` errors on 38
assertions the moment brms became available. Same class as §K.1 (stale
`"bayes"` / `Lower_CL` assertions) and §K.11 generally: the fix landed, the test
did not follow, and nothing failed because nothing ran.

### R3.35 `bayesian_synergy()` and `bayesian_synergy_over_time()` failed on every call — **Critical (pre-existing)**
**File:** `R/bayesian_synergy.R:710, 1080`

Both exported functions referenced a bare `re_term` while assembling their
summary metadata. `re_term` is built inside `bs_fit_synergy_model()` and returned
as `fit$re_term` — the Round 2 §G.6 refactor (v0.4.6) that extracted the shared
helper moved the construction but left both consumers pointing at the old local
variable. Every call died with `object 're_term' not found` *after* the Stan fit
completed, so the user paid the full 3–12 minute fit cost and then got an error.

**`bayesian_synergy()` has therefore been completely non-functional since
v0.4.6** — through five releases, while the dashboard advertised a Bayesian
synergy mode. It was invisible because `test-bayesian_synergy.R` skips wholesale
without brms, and brms was not installed in the development environment. On
installing brms the file went from 2 passing / 24 errors to 45 passing / 0 errors
with the one-line fix.

### R3.36 The body-weight fixture had zero between-animal variance — **Test defect**
**File:** `tests/testthat/helper-fixtures.R` (`make_bw_simple()`)

Every mouse started at exactly 22 g, so the fixture contained no between-animal
variation at all. That makes the random intercept under
`random_effects_specification = "intercept_only"` unidentifiable — its true SD is
0, the sampler explores an Intercept↔sd funnel, and Intercept ESS collapses
(Bulk 281 / Tail 71 at 3000 draws) while every other parameter sits above 900.

This is the mechanism behind the ESS failures that R3.34 first surfaced, and it
explains why the old mis-located Intercept prior *passed*: a prior tight around
the wrong value pins the intercept and suppresses the funnel. A correctly-located
wider prior lets the intercept move, re-opening it. The old prior sampled better
precisely because it was informative and wrong.

Fixed by giving the fixture a realistic 0.8 g per-animal SD. A body-weight
fixture with no between-animal variation cannot exercise the very random effect
the tests assert on. Body-weight and toxicity suites both clean afterwards
(32/0/0 and 70/0/0).

### R3.37 `bayesian_therapeutic_window()`'s efficacy-vs-safety plot was always NULL — **Major (pre-existing)**
**File:** `R/bayesian_therapeutic_window.R:383-384`

The TWM=1 isoline block computed its x-range from `plot_df`, a data frame that is
never created anywhere in the function — the frame is called `scatter_df`. The
undefined variable threw inside the enclosing `tryCatch`, so `tgi_wl_plot` was
**always** `NULL` and the dashboard's TWM scatter tab was permanently empty.

Note the shape: Round 2 §J.2 flagged this very isoline as "approximate" and it
was reportedly fixed in v0.4.5 by replacing `geom_abline` with an explicit
`geom_path`. That fix introduced the `plot_df` reference. So §J.2 is a third
instance of the pattern in §R3-D5 — a fix that was written, marked resolved, and
never executed.

### R3.34 The first R3.8 prior recalibration was wrong, in the opposite direction — **caught by the tests**

Worth recording because it is the same defect I criticised, mirrored. The initial
fix reinterpreted `b_sd` as a prior SD on the *total log-fold change over the
study* and divided it by the time span to get the slope scale. That reasoning
holds on a log scale but not on a raw one, and reusing the main-effect ladder
multiplier (0.25 for "skeptical") made the slope prior 5–20× too tight. On the
package's own body-weight fixture the true interaction is 0.14 g/day and the
prior SD came out at 0.0089 — **16 prior SDs from the truth**, versus the ~9 SDs
that made the original Intercept prior a finding in the first place.

It surfaced as Bulk_ESS / Tail_ESS falling below 400 on the body-weight fits: the
prior was fighting the likelihood, degrading the geometry. Attribution was
confirmed by stashing the change and re-running (32/0/0 before, failing after) —
not inferred.

The corrected form gives rate coefficients their own scale unit — the slope that
would traverse the entire observed response range over the study, which no real
slope can much exceed — with its own multipliers (skeptical 1.0 → diffuse 5.0).
Verified against known truth on both a gram-scale response and a log-volume
response: "skeptical" now sits 1.2 and 0.2 prior SDs from the true interaction
respectively, rather than 21 and 2.8.

**Lesson for the file:** a data-scaled prior is not automatically a well-calibrated
one. The intercept fix was unambiguous (centre on the data); the slope fix needed
its own scale reasoning and empirical checking against known effect sizes, and
the first attempt failed it. Both directions of mis-scaling are defects.

---

## R3-K. Round 3 closed — final implementation log (v0.9.0, 2026-07-29)

Every Round 3 finding is now closed except **H.3**, which remains blocked on the
randomisation-unit question (see below). Shipped across five releases:

| Release | Content |
|---|---|
| v0.5.0 | R3.2 R3.10 R3.11 R3.14 R3.16 R3.18 R3.19 R3.20 R3.23 R3.25 R3.27 R3.29 R3.30 R3.31; R3.1 (survival) |
| v0.6.0 | R3.1 (all paths) R3.4 R3.12 R3.17 R3.21; G.1 + G.2 |
| v0.7.0 | R3.6 R3.7 R3.13 R3.15 R3.22 R3.24 R3.26 R3.28; G.6 Option A |
| v0.8.0 | R3.3 R3.5; G.3 |
| v0.9.0 | R3.8 R3.9; G.4 + G.5 + H.2 + H.4 |

### What the numbers moved by

Three findings changed results by amounts worth restating, because anything
published from the affected paths needs re-running:

1. **R3.5 (survivor bias).** Against ground truth on a simulation with the real
   dropout mechanism — per-animal growth rates, removal at a 2000 mm³ limit,
   60 % of controls lost by day 28 — mean absolute TGI error was **23.3
   percentage points** under the old survivor-conditional estimand and **5.3**
   under the model-based one. `last_obs` scored 28.9, i.e. worse than the old
   behaviour, which is why it is documented as a fallback rather than an equal
   option.
2. **R3.17 (cage).** On a 2-cages-per-arm, 4-mice-per-cage design with a real
   cage effect, restoring the estimable `(1|cage)` term moved treatment-contrast
   SEs from 0.141 to 0.305 — a factor of **2.2**. At the more typical ICC of
   0.1–0.2 the factor is ~1.15–1.25.
3. **R3.1 (multiplicity).** The default `lme4` path reported unadjusted p-values
   while documenting Bonferroni. With three vs-control contrasts, adjusted
   p-values are 3× the previously-reported ones under Bonferroni.

### Two fixes that had never executed

R3.1 and R3.2 were both written in earlier rounds and never ran — R3.2 because a
formula referenced a full-length `Surv` object against subsetted data and the
error was swallowed into a fallback, R3.1 because `emmeans::contrast()` defaults
to no adjustment for a coefficient list. Both are now covered by tests asserting
observable behaviour rather than the presence of code. That distinction is the
single most transferable lesson from this round and is recorded in
[[code-review-rounds]] terms in §R3-D5: **verify a claimed fix by running it.**

### Test surface

- **Backend:** 352 passing / 7 failing at the start of Round 3 → all green, with
  the suite now genuinely exercising the Bayesian paths (`brms` 2.23) and the
  permutation paths (`coin` 1.4.5) that previously skipped. The seven baseline
  failures were §K.11 stale tests referencing APIs removed in v0.3.4 and earlier;
  clearing them also unblocks the `covr` baseline §K.10 flagged.
- **Dashboard:** 210 passing / 0 failing.
- Every Round 3 fix has a regression test in
  `tests/testthat/test-code_review_round3.R`, labelled with its finding ID —
  closing the §K.11 process gap for this round.

### The one item still open

**H.3 — permutation test of the global "no treatment effect on trajectory"
null.** Still blocked, and deliberately so. A permutation test is a
*randomisation* test and must mirror how randomisation was actually performed:

> Are individual mice randomised to treatment and then re-housed by group, or are
> cages assigned to treatments?

If mice, the permutation unit is the mouse and 48 animals give ample resolution —
but cage then becomes a post-randomisation variable, and permuting mice across
cages breaks the structure §G.2 establishes is real. If cages, the unit is the
cage, and with 12 cages across 6 arms the effective sample size for the test is
12, not 48. That second number is uncomfortable and also honest: if cages were
the unit of assignment, 12 is the real replication and the analysis overstates
precision whether or not a permutation test is ever run.

Implementing this on a guess would produce a test that looks rigorous and is
calibrated for the wrong design, which is worse than not having it.

### Recommended follow-ups (not Round 3 findings)

- **Run the `covr` baseline** now that §K.11 is clear — §K.10's blocker is gone.
- **Add `coin` (and `clinfun`) to the VPS Shiny image.** Both are `Suggests`, so
  the permutation log-rank and the Jonckheere-Terpstra trend test degrade to a
  skip-with-message when absent. They will silently not run in production
  otherwise.
- **Re-check the Bayesian prior-posterior plots on a real mm³ dataset.** The
  data-scaled Intercept prior (R3.8) should remove the prior/data conflict that
  was previously visible on every real study; that plot is now a meaningful
  sensitivity check rather than a guaranteed red flag.
- **B7.1 / B7.2** (Bayesian vs frequentist result-schema harmonisation) remain
  open from Round 1 and are now the largest consistency gap left.
