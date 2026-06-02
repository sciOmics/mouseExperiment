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
