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
| B1.1 | `model_type_used = "bayes"` should be `"bayes_tg"` | Major | `bayesian_tumor_growth.R:579` | Open |
| B1.2 | `prior_strength` default is `"weakly_informative"` in body_weight, `"skeptical"` elsewhere | Major | `bayesian_body_weight.R:116` | Open |
| B1.3 | `"informative"` preset missing from `bayesian_dose_response` | Major | `bayesian_dose_response.R:123` | Open |
| B1.4 | `Lower_CL`/`Upper_CL` should be `Lower_CrI`/`Upper_CrI` in Bayesian functions | Major | TG, BW, Survival | Open |
| B1.5 | `include_frailty` vs `include_cage_effect` — same concept, different name | Minor | `bayesian_survival.R:149` | Open |
| B1.6 | Intercept prior width: `b_sd * 2.5` in 3 files, `b_sd * 2.0` in synergy | Minor | `bayesian_synergy.R:287` | Open |
| B1.7 | `set_prior()` in synergy vs `prior_string()` everywhere else | Minor | `bayesian_synergy.R:276` | Open |
| B2.1 | `bt()` back-transform helper duplicated in 3 files | Minor | BW, Synergy, TWM | Open |
| B2.2 | MCMC diagnostics DF construction duplicated in 4+ files; schema differs in synergy | Minor | Multiple | Open |
| B2.3 | Prior specification switch duplicated in 5 files | Minor | Multiple | Open |
| B2.4 | Cage placeholder setup duplicated in 4 files | Minor | Multiple | Open |
| B3.1 | Loewe 1e-6 floor arbitrary and undocumented; produces huge CI for near-zero combo effects | Major | `bayesian_synergy.R:413` | Open |
| B3.2 | `bayesian_synergy` MCMC diagnostics missing Bulk_ESS and Tail_ESS | Minor | `bayesian_synergy.R:317` | Open |
| B3.3 | All-zero/negative volume with log transform silently produces `log(Inf)` | Major | TG, BW, Synergy | Open |
| B3.4 | TWM group-name intersection not normalised for whitespace/case | Minor | `bayesian_therapeutic_window.R:225` | Open |
| B3.5 | `bliss_summary` inconsistent rounding (3 vs 4 dp) | Minor | `bayesian_synergy.R:402` | Open |
| B4.1 | Reference-group auto-detection not documented | Minor | `bayesian_tumor_growth.R` | Open |
| B4.2 | Weight-loss endpoints not documented in `bayesian_body_weight` | Minor | `bayesian_body_weight.R` | Open |
| B4.3 | `bayesian_synergy` missing `@section Assumptions and Limitations:` | Minor | `bayesian_synergy.R` | Open |
| B4.4 | `bayes_prior_posterior_plot()` cross-file dependency not indicated | Minor | `bayesian_tumor_growth.R` | Open |
| B5.1 | `bayesian_body_weight` calls `posterior_epred()` 3× instead of 2× | Minor | `bayesian_body_weight.R:416` | Open |
| B6.1 | No `bayesian_power_analysis()` function | Enhancement | — | Open |
| B6.2 | No single-call wrapper for TWM (requires two pre-fitted models) | Enhancement | `bayesian_therapeutic_window.R` | Open |
| B6.3 | No `bayesian_synergy_over_time()` | Enhancement | — | Open |
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
