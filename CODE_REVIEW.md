# Code Review: mouseExperiment R Package

**Date:** Updated April 1, 2026  
**Package:** mouseExperiment v0.1.0  
**Scope:** Deep review of all 19 R source files (~7,060 lines)

---

## Executive Summary

mouseExperiment is a specialized R package for statistical analysis of mouse tumor growth experiments. It provides LME4 mixed-effects modeling, AUC analysis, survival analysis, dose-response analysis, drug synergy assessment, and post-hoc power analysis with associated visualizations.

**April 1 update:** All 5 confirmed bugs (B1–B5) are fixed. All deprecated API calls (aes_string, group_by_at, size→linewidth) migrated. NAMESPACE cleaned up with selective importFrom. New S3 class `me_result` with utility functions. Real Monte Carlo power simulation. AUC consolidated into shared utility. 24 new tests passing. Primary remaining technical debt: ~46 `cat()` calls that should be `message()`, and two very large functions (~1,000+ lines each) that could benefit from sub-function extraction.

**Overall assessment:** Production-ready. Remaining items are code hygiene improvements, not correctness issues.

---

## 1. Confirmed Bugs

These were identified through independent manual validation (see VALIDATION.md) and verified by code inspection.

### B1: `log(Volume + 1)` Instead of `log(Volume)` — MEDIUM

**File:** tumor_growth_statistics.R ~L297  
**Impact:** Inflates chi-square values in cage-treatment collinearity test; minor effect on LME4 model coefficients

The log transform uses `log(Volume + 1)` (via `log1p`-like behavior), but the package documentation and model interpretation assume `log(Volume)`. The `+1` constant biases small volumes upward and affects the chi-square test on the cage-treatment contingency table because the transform is applied before that test is run on the transformed data.

```r
# Current (line ~297)
analysis_df[[volume_column]] <- log(analysis_df[[volume_column]] + 1)

# Should be
analysis_df[[volume_column]] <- log(analysis_df[[volume_column]])
```

**Note:** Zero volumes need separate handling (e.g., replace with a small epsilon before log, or use `log(max(V, 0.01))`).

---

### B2: Growth Rates Double Log-Transformed — HIGH

**File:** tumor_growth_statistics.R ~L316  
**Impact:** All growth rate values are wrong (double-logged)

After the volume column is already log-transformed (B1), the growth rate calculation applies `log1p()` again:

```r
# Line ~316: volume_column is ALREADY log-transformed at this point
log_volume <- log1p(subject_data[[volume_column]])
```

This means growth rates are computed as `d/dt[log(1 + log(V+1))]` instead of `d/dt[log(V)]`. The resulting values are mathematically incorrect — they underestimate true exponential growth rates significantly.

**Fix:** Either:
- Compute growth rates from the pre-transformed `auc_df` data (which retains original volumes), or
- Remove the `log1p()` call since `subject_data[[volume_column]]` is already log-transformed.

---

### B3: Pairwise Comparisons Group Ordering Wrong — MEDIUM

**File:** tumor_growth_statistics.R ~L840–L899  
**Impact:** LME4 pairwise comparison labels may have wrong sign/direction

The emmeans contrasts are built using `treatment_groups <- unique(df[[treatment_column]])`, which returns groups in dataset encounter order. However, emmeans internally uses alphabetical factor level order. When these disagree, the contrast coefficient vector `contrast_coef` is assigned to the wrong groups, producing comparisons with swapped signs (e.g., "B − A" labeled as "A − B").

```r
# Current: uses unique() order
treatment_groups <- unique(df[[treatment_column]])
# ...
ref_idx <- which(treatment_groups == reference_group)
group_idx <- which(treatment_groups == group)
contrast_coef[ref_idx] <- -1
contrast_coef[group_idx] <- 1
```

**Fix:** Sort `treatment_groups` alphabetically to match emmeans internal ordering, or use `levels(factor(...))`.

---

### B4: Warning Handler Discards Valid Models — LOW-MEDIUM

**File:** tumor_growth_statistics.R ~L428–L445  
**Impact:** Valid models with boundary warnings are silently discarded

Both the intercept-only and random-slope model fitting blocks use `tryCatch` with a `warning` handler that returns `NULL` on boundary (singular) fit warnings:

```r
}, warning = function(w) {
    if (grepl("boundary", w$message)) {
      warning("Boundary (singular) fit detected in intercept-only model")
    }
    return(NULL)  # discards the model even though it converged
})
```

A boundary warning does not mean the model failed — it means the estimated variance of a random effect is near zero. The model is still valid and usable. Discarding it means both candidate models may be NULL, triggering the fallback path which has its own bug (B5).

**Fix:** Use `withCallingHandlers` to capture the warning while still returning the model, or use `tryCatch` with `warning = function(w) invokeRestart("muffleWarning")` after recording the warning.

---

### B5: `best_model` Undefined in Fallback Path — LOW

**File:** tumor_growth_statistics.R ~L450–L460  
**Impact:** Silent failure in diagnostics when both models are discarded (B4)

When both candidate models fail (are set to NULL by B4), the fallback block fits a new model but never assigns the local variable `best_model`. The diagnostics section then references `best_model` to decide whether to extract random slopes:

```r
diagnostics$random_effects <- list(
  intercepts = lme4::ranef(model)[[id_column]],
  slopes = if (best_model == "random_slope") {
    lme4::ranef(model)[[id_column]][, 2]
  } else NULL
)
```

This will error with "object 'best_model' not found" if the fallback path was taken.

**Fix:** Set `best_model <- "intercept_only"` in the fallback block.

---

## 2. Code Quality Issues

### 2.1 Function Size and Complexity

| Function | File | Lines | Issue |
|----------|------|-------|-------|
| `tumor_growth_statistics` | tumor_growth_statistics.R | 975 | Handles LME4, AUC, model selection, diagnostics, growth rates, cage effects, emmeans contrasts, extrapolation — all in one function |
| `post_power_analysis` | post_power_analysis.R | 1,218 | Monolithic: input parsing, AUC calculation, effect size estimation, power curves, sample size estimation, plotting |
| `survival_statistics` | survival_statistics.R | 710 | Combines Cox, Firth, log-rank, median survival calculation, results formatting |
| `dose_response_statistics` | dose_response_statistics.R | 608 | Orchestrator calling 5 internal helpers — better structured than above |
| `plot_tumor_growth` | plot_tumor_growth.R | 383 | Contains ~150 lines of extrapolation logic embedded in a plotting function |

**Recommendation:** Extract into focused sub-functions. `tumor_growth_statistics` should dispatch to `tgs_fit_lme4()`, `tgs_fit_auc()`, `tgs_compute_growth_rates()`, `tgs_compute_cage_effects()`, etc.

---

### 2.2 Deprecated ggplot2 Usage — ✅ FIXED

All `size` aesthetics for line geoms have been replaced with `linewidth` across all 5 affected plot files. Point geoms correctly retain `size`.

---

### 2.3 `aes_string()` Usage — ✅ FIXED

All `aes_string()` calls in dose_response_statistics.R have been replaced with `aes(.data[[col]])` pattern. Zero occurrences remain in the package.

---

### 2.4 Logic Bug in plot_caterpillar.R — ✅ FIXED

The Main Effect assignment now uses `grepl() & !grepl()` (both logical vectors). Interaction and Intercept assignments use `grep()` for index-based subsetting (single condition, no `&` mixing), which is correct.

---

### 2.5 Hardcoded Group Names in plot_bliss.R — ✅ FIXED

Group labels are now parameterized via `drug_a_label`, `drug_b_label`, `combo_label`, and `expected_label` arguments with sensible defaults.

---

### 2.6 Performance: `rbind()` in Loop — ✅ FIXED

Extrapolation in plot_tumor_growth.R now collects rows in `vector("list")` and uses `do.call(rbind, ...)`. Same pattern applied in tumor_growth_statistics.R and post_power_analysis.R.

---

### 2.7 Side Effects in Plot Functions — ✅ FIXED (plot_growth_rate.R)

`cat()` calls in plot_growth_rate.R replaced with `message()`. However, **`cat()` calls remain** in other non-plot files — see 2.13.

---

### 2.8 CI Calculation Uses `qnorm` Instead of `qt` — ✅ FIXED

plot_caterpillar.R now uses `qt(1 - (1 - ci_level)/2, df = df_resid)` with degrees of freedom extracted from the model.

---

### 2.9 `in_place` Parameter Anti-Pattern — ✅ FIXED

Both `calculate_volume.R` and `calculate_dates.R` now emit a deprecation warning when `in_place = TRUE` is used, and default to `in_place = FALSE`. The `assign()` into parent frame code is retained for backward compatibility but warned against.

---

### 2.10 Identical Non-Linear Models in dose_response_statistics.R — ✅ FIXED

`dr_model_decr` now uses `LL.4()` (4-parameter) and `dr_model_incr` uses `LL.5()` (5-parameter with asymmetry), providing genuinely different model fits.

---

### 2.11 Loewe Additivity Oversimplified — ✅ PARTIALLY FIXED

User-facing labels and column names have been renamed from "Loewe" to "Additive (Mean)" — variables `additive_mean_tgi`, `additive_mean_difference`, columns `Additive_Mean_Expected_TGI`, `Additive_Mean_Difference`. A clarifying comment explains this is not true Loewe Additivity.

**Remaining:** The return list key is still `loewe_additivity` (API compatibility) and some roxygen comments still reference "Loewe additivity". Renaming the key would be a breaking API change. ⏳

---

### 2.12 Deprecated dplyr Usage — ✅ FIXED

All `group_by_at(vars(...))` and `arrange_at()` calls in dose_response_statistics.R replaced with `group_by(.data[[col]])` and `arrange()`. Zero deprecated dplyr calls remain.

---

### 2.13 `cat()` vs `message()` for Output — ⏳ PARTIALLY FIXED

Fixed in plot_growth_rate.R (converted to `message()`). **Remaining `cat()` calls (~46):**

| File | Count | Context |
|------|-------|---------|
| tumor_growth_statistics.R | ~18 | Mostly `if (verbose)` gated — acceptable but should be `message()` |
| analyze_drug_synergy.R | ~18 | Result summary printing |
| analyze_drug_synergy_over_time.R | ~8 | Result summary printing |
| dose_response_statistics.R | ~15 | Summary output block (L485–L554) |

Additionally, dose_response_statistics.R has ~8 bare `print()` calls (L228–L443) for model summaries that should use `message(paste(capture.output(...)))` pattern.

**Recommendation:** Convert all `cat()`→`message()` and wrap `print()` in `capture.output()` → `message()`. Gate behind `verbose` parameter where appropriate.

---

### 2.14 "Simulation" Power Method Is Not a Simulation — ✅ FIXED

The `method = "simulation"` path now runs real Monte Carlo simulations: generates `rnorm()` data per group, runs `t.test()` per simulation, and computes empirical power as the fraction of rejections. The `method = "analytical"` path correctly uses `power.t.test()` for closed-form computation.

---

## 3. Architecture Issues

### 3.1 NAMESPACE: Full Package Imports — ✅ FIXED

NAMESPACE converted from wholesale `import()` to selective `importFrom()` directives: dplyr (8), ggplot2 (44), stats (26), survival (4), drc (3), rlang (1), utils (1). Dead `import(survminer)` removed.

### 3.2 Heavy Dependency Tree — ✅ FIXED

`anytime`, `clinfun`, and `ggpubr` moved from Imports to Suggests. `calculate_dates.R` now uses `requireNamespace("anytime")` guard.

### 3.3 Duplicated AUC Calculation — ✅ FIXED

Consolidated into a single exported `calculate_auc()` in `R/utils_auc.R`. All three former inline implementations (tumor_growth_statistics.R, tumor_auc_analysis.R, post_power_analysis.R) now call this shared function.

### 3.4 Duplicated Composite ID Pattern — ✅ FIXED

All composite ID creation now uses `sep = "|||"` with `strsplit(..., fixed = TRUE)` for safe parsing. Applied consistently in tumor_growth_statistics.R (4 locations), tumor_auc_analysis.R, and post_power_analysis.R.

### 3.5 No Formal Test Suite — ✅ PARTIALLY FIXED

`tests/testthat/test-utils_and_me_result.R` added with 24 passing tests covering `calculate_auc()`, `new_me_result()`, `print/summary/plot.me_result`, `export_diagnostics()`, `tumor_doubling_time()`, and `repeated_measures_anova()`.

**Remaining:** No tests for existing functions (tumor_growth_statistics, survival_statistics, dose_response_statistics, plot functions). Coverage is limited to new utility functions.

### 3.6 Composite ID Parsing Fragility — ✅ FIXED

All composite IDs now use `"|||"` separator with `strsplit(..., fixed = TRUE)`, eliminating collisions from underscored names.

---

## 4. Enhancement Opportunities

### 4.1 Additional Functionality

| Feature | Description | Priority | Status |
|---------|-------------|----------|--------|
| **Repeated-measures ANOVA** | Alternative to LME4 for simpler designs | Medium | ✅ `repeated_measures_anova()` in me_result.R |
| **Body weight analysis** | Toxicity assessment from weight data | High | ⏳ Deferred |
| **Tumor doubling time** | Common metric not currently computed | Low | ✅ `tumor_doubling_time()` in me_result.R |
| **Multiple comparison methods** | Add Holm, FDR in addition to Bonferroni/Tukey | Medium | ✅ `p_adjust_method` parameter (bonferroni/holm/fdr/none) |
| **Model diagnostics export** | Return formal diagnostic test results (Shapiro-Wilk, Levene's) | Medium | ✅ `export_diagnostics()` in me_result.R |
| **Bayesian analysis option** | brms-based mixed models as alternative to lme4 | Low | ⏳ Deferred |
| **Formal simulation power** | Real Monte Carlo power analysis to replace fake "simulation" method | Medium | ✅ Real Monte Carlo in post_power_analysis.R |

### 4.2 API Improvements

1. ✅ **Consistent return structures:** S3 `me_result` class provides standardized `print`/`summary`/`plot` methods for analysis results. Existing functions can progressively adopt `new_me_result()`.

2. ✅ **S3 class for results:** `me_result` class defined in `me_result.R` with `print.me_result`, `summary.me_result`, `plot.me_result` methods plus `export_diagnostics()`.

3. **Formula interface:** Allow users to specify the model formula directly (e.g., `formula = Volume ~ Day * Treatment + (Day | ID)`) as an alternative to the parameter-based interface. *(Deferred — would require significant refactoring of tumor_growth_statistics.R.)*

4. **Progress reporting:** Long-running functions (power analysis with simulations) should support progress callbacks, especially when called from Shiny. *(Mitigated by `withSpinner()` in the dashboard.)*

---

## 5. File-by-File Summary

| File | Lines | Status | Key Issues |
|------|-------|--------|------------|
| tumor_growth_statistics.R | 981 | Good | B1–B5 ✅; composite ID ✅; p_adjust_method ✅; `cat()` remains (verbose-gated) |
| survival_statistics.R | 709 | OK | Well-structured adaptive method; could use refactoring |
| dose_response_statistics.R | 609 | Warn | aes_string ✅; group_by_at ✅; LL.4/LL.5 ✅; **`print()`/`cat()` calls remain** |
| post_power_analysis.R | 1,219 | Good | Real Monte Carlo ✅; AUC consolidated ✅; still very long |
| analyze_drug_synergy.R | 367 | Good | Additive rename ✅; `loewe_additivity` key retained for API compat; `cat()` remains |
| analyze_drug_synergy_over_time.R | 513 | Good | Column names updated ✅; `cat()` remains |
| tumor_auc_analysis.R | 404 | Good | Uses shared `calculate_auc()` ✅ |
| calculate_dates.R | 189 | Good | `in_place` deprecated ✅; `requireNamespace("anytime")` guard ✅ |
| calculate_volume.R | 135 | Good | `in_place` deprecated ✅ |
| data.R | 110 | Good | Clean dataset documentation |
| plot_tumor_growth.R | 385 | Good | rbind ✅; linewidth ✅; still long (extrapolation logic embedded) |
| plot_growth_rate.R | 259 | Good | cat→message ✅; linewidth ✅ |
| plot_auc.R | 274 | Good | linewidth ✅ |
| plot_bliss.R | 125 | Good | Parameterized labels ✅; linewidth ✅ |
| plot_caterpillar.R | 140 | Good | grep/grepl ✅; qt() ✅ |
| plot_combination_index.R | 142 | Good | linewidth ✅ |
| plot_treatments.R | 121 | Good | Cleanest plot file; modern `.data[[]]` usage |
| utils_auc.R | 46 | Good | **New** — consolidated trapezoidal AUC |
| me_result.R | 337 | Good | **New** — S3 class + utility functions |

---

## 6. Priority Action Items

### Critical (Fix Before Release) — ✅ ALL RESOLVED
1. ✅ **B2:** Fixed double-log growth rate calculation (commit `e3a0bb6`)
2. ✅ **B3:** Fixed emmeans contrast ordering (commit `074e52e`)
3. ✅ **B1:** Uses `log(Volume)` with proper zero handling (commit `e3a0bb6`)

### High Priority — ✅ ALL RESOLVED
4. ✅ **B4/B5:** Fixed warning handler and `best_model` fallback (commit `e3a0bb6`)
5. ✅ Consolidated 3 AUC implementations into `utils_auc.R::calculate_auc()`
6. ✅ Replaced deprecated `aes_string()` and `group_by_at()` in dose_response_statistics.R
7. ✅ Fixed `grep()` vs `grepl()` logic bug in plot_caterpillar.R
8. ✅ Added testthat tests for new functions (`test-utils_and_me_result.R`, 24 tests passing)

### Medium Priority
9.  ⏳ Refactor `tumor_growth_statistics()` into sub-functions (deferred — function is large but stable)
10. ⏳ Refactor `post_power_analysis()` into sub-functions (deferred — function is large but stable)
11. ✅ Replaced deprecated `size` with `linewidth` across all 5 plot files
12. ✅ Standardized return structures via S3 `me_result` class (`me_result.R`)
13. ✅ Converted full `import()` to specific `importFrom()` in NAMESPACE (dplyr:8, ggplot2:44, stats:26, drc:3)
14. ✅ Moved `anytime`, `clinfun`, `ggpubr` to Suggests in DESCRIPTION
15. ⏳ Convert remaining `cat()` → `message()` (~46 calls across 4 files) — see 2.13
16. ⏳ Convert `print()` → `message(capture.output())` in dose_response_statistics.R (~8 calls)

### Low Priority
17. ✅ `in_place` parameter deprecated with warning in `calculate_dates()` and `calculate_volume()`
18. ✅ S3 class `me_result` with `print`, `summary`, `plot` methods (`me_result.R`)
19. ✅ Renamed Loewe → "Additive (Mean)" in user-facing labels; `loewe_additivity` key retained for API compat
20. ✅ Implemented real Monte Carlo simulation for power analysis (`post_power_analysis.R`)
21. ⏳ `.Rbuildignore` improvements (deferred — low impact)
22. ⏳ Rename `loewe_additivity` return key to `additive_model` (breaking API change — defer to next major version)
23. ⏳ Add tests for existing functions (tumor_growth_statistics, survival, plot functions)

### New Functions Added
- `calculate_auc()` — Consolidated vectorized trapezoidal AUC (`utils_auc.R`)
- `new_me_result()` — S3 constructor for standardized analysis results (`me_result.R`)
- `export_diagnostics()` — Export model diagnostics to CSV/data frame (`me_result.R`)
- `tumor_doubling_time()` — Per-subject exponential growth doubling time (`me_result.R`)
- `repeated_measures_anova()` — Treatment × Time interaction via lmerTest (`me_result.R`)
- `palette_colors()` — Named colorblind-friendly palette presets (dashboard `helpers.R`)

---

*Review based on complete source code examination and independent validation against all 4 demo datasets.*
