# Code Review: mouseExperiment R Package

**Date:** Updated March 25, 2026  
**Package:** mouseExperiment v0.1.0  
**Scope:** Deep review of all 17 R source files (~6,800 lines)

---

## Executive Summary

mouseExperiment is a specialized R package for statistical analysis of mouse tumor growth experiments. It provides LME4 mixed-effects modeling, AUC analysis, survival analysis, dose-response analysis, drug synergy assessment, and post-hoc power analysis with associated visualizations. The package is well-documented and follows R package conventions, but contains **5 confirmed bugs** (validated against independent manual calculations), significant **function complexity issues**, and several areas where code quality and architecture can be improved.

**Overall assessment:** Functional with known bugs. Prioritize bug fixes, then refactoring.

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

### 2.2 Deprecated ggplot2 Usage

**`size` aesthetic deprecated in ggplot2 >= 3.4** — affects 5 of 7 plot files:

| File | Lines | Issue |
|------|-------|-------|
| plot_auc.R | L175, L199 | `size = 1` in `stat_summary`/`geom_segment` — use `linewidth` |
| plot_bliss.R | L78, L79, L98 | `aes(size = Group)` + `scale_size_manual` — use `linewidth`/`scale_linewidth_manual` |
| plot_combination_index.R | L78 | `geom_line(size = 1.5)` — use `linewidth` |
| plot_growth_rate.R | L201 | `geom_hline(size = 0.5)` — use `linewidth` |
| plot_tumor_growth.R | L296, L330, L341, L346 | Multiple `size` → `linewidth` in `geom_line`/`stat_summary` |

---

### 2.3 `aes_string()` Usage

`aes_string()` is deprecated since ggplot2 3.0.0 (2018). Found in:
- dose_response_statistics.R: `create_dose_plots()` uses `aes_string()` throughout
- Should migrate to `aes(.data[[col]])` pattern (already used in plot_auc.R, plot_treatments.R, plot_tumor_growth.R)

---

### 2.4 Logic Bug in plot_caterpillar.R

Line ~87: `grep()` returns integer indices, `grepl()` returns logical vector. Using `&` between them does not work as intended:

```r
effect_type[grep(main_effect_pattern, coef_names) & !grepl("^\\(Intercept\\)$", coef_names)] <- "Main Effect"
```

**Fix:** Use `grepl()` for both to get same-length logical vectors.

---

### 2.5 Hardcoded Group Names in plot_bliss.R

Group labels "Drug A", "Drug B", "Combination", "Bliss Value" are hardcoded. If `analyze_drug_synergy()` is called with custom drug names, the plot labels won't match.

---

### 2.6 Performance: `rbind()` in Loop

plot_tumor_growth.R L260–L275: Grows a data frame row-by-row inside a nested loop — O(n^2) performance. Should collect in a list and call `do.call(rbind, ...)` or use `dplyr::bind_rows()`.

---

### 2.7 Side Effects in Plot Functions

plot_growth_rate.R L161–L165: Uses `cat()` to print mouse counts. Plot functions should not have console side effects. Use `message()` or remove.

---

### 2.8 CI Calculation Uses `qnorm` Instead of `qt`

plot_caterpillar.R L72: Uses `qnorm(0.975)` (1.96) for confidence intervals on all model types, including finite-sample models where `qt()` with appropriate degrees of freedom would be more accurate.

---

### 2.9 `in_place` Parameter Anti-Pattern

Both `calculate_volume.R` and `calculate_dates.R` implement `in_place = TRUE` by using `assign()` into the parent frame:

```r
if (in_place) {
    parent_frame <- parent.frame()
    df_name <- deparse(substitute(df))
    assign(df_name, result, envir = parent_frame)
}
```

This is fragile — it fails when the input is a complex expression (e.g., `calculate_volume(my_list$data)`), inside pipes, or when called from a different environment depth. The function also always returns the result, making the `in_place` flag redundant in practice since callers can just do `df <- calculate_volume(df)`.

**Recommendation:** Deprecate `in_place` parameter. The standard R pattern is to return a modified copy.

---

### 2.10 Identical Non-Linear Models in dose_response_statistics.R

`try_nonlinear_models()` fits two models (`dr_model_decr` and `dr_model_incr`) using identical formulas and `LL.4()` function — they will always produce the same result. The intent was likely to try a decreasing vs increasing curve, but the implementation doesn't differentiate them.

---

### 2.11 Loewe Additivity Oversimplified

`analyze_drug_synergy.R` implements Loewe additivity as `(Effect_A + Effect_B) / 2`. This is a simple average, not the Loewe additivity model, which requires dose-response curve parameters. The code comment acknowledges this ("simplified approach without dose information"), but the result is labeled "Loewe Additivity" in output, which could mislead users.

**Recommendation:** Either implement proper Loewe additivity (requires dose-response curves) or rename to "Additive Model (mean)" to avoid confusion.

---

### 2.12 Deprecated dplyr Usage

dose_response_statistics.R uses `dplyr::group_by_at()` and `dplyr::vars()`, deprecated since dplyr 1.0.0. Should migrate to `group_by(across(...))` or `group_by(.data[[col]])`.

---

### 2.13 `cat()` vs `message()` for Output

Multiple functions use `cat()` for console output (dose_response_statistics.R, survival_statistics.R, generate_user_report, etc.). `cat()` output cannot be suppressed with `suppressMessages()` and is not capturable in testing. Use `message()` for informational output.

---

### 2.14 "Simulation" Power Method Is Not a Simulation

post_power_analysis.R `method = "simulation"` (~L870) contains:
```r
# In a real implementation, this would run actual simulations
```
It actually calls `power.t.test()` — identical to the parametric method. Either implement actual Monte Carlo simulation (resample data, fit models, count rejections) or remove the misleading method option.

---

## 3. Architecture Issues

### 3.1 NAMESPACE: Full Package Imports

The NAMESPACE uses `import(dplyr)`, `import(ggplot2)`, `import(stats)`, `import(survival)`, `import(survminer)`, `import(drc)` — importing entire namespaces. This creates potential function name conflicts:

- `dplyr::filter()` masks `stats::filter()`
- `dplyr::lag()` masks `stats::lag()`
- `drc::gaussian()` masks `stats::gaussian()` (visible in runtime warnings)

**Recommendation:** Convert `import()` to specific `importFrom()` calls, especially for dplyr and stats.

### 3.2 Heavy Dependency Tree

The package has 22+ direct dependencies in Imports. Several are only used in one function:
- `anytime`: only used in `calculate_dates()` fallback
- `clinfun`: only used for Jonckheere-Terpstra test (currently disabled with comment "not run due to mentioned issues")
- `ggpubr`: only used in `plot_synergy_combined()` for `ggarrange()`
- `rlang`: only used for `sym()` and `!!` in 2 files

**Recommendation:** Move rarely-used dependencies to Suggests and use `requireNamespace()` checks.

### 3.3 Duplicated AUC Calculation

AUC (trapezoidal rule) is implemented in three places:
1. `tumor_growth_statistics.R` ~L540–L560 (inline `calculate_auc` function)
2. `tumor_auc_analysis.R` ~L170–L200 (`calculate_subject_auc` method)
3. `post_power_analysis.R` ~L1100–L1150 (`calculate_auc_values` function)

All three implement the same trapezoidal rule but with slightly different handling of extrapolation, composite IDs, and edge cases. Any bug fix must be applied in three places.

**Recommendation:** Consolidate into a single exported `calculate_auc()` utility.

### 3.4 Duplicated Composite ID Pattern

Creating composite IDs from `paste(ID, Treatment, Cage, sep = "_")` appears in:
- tumor_growth_statistics.R (AUC section)
- tumor_auc_analysis.R
- post_power_analysis.R (`calculate_auc_values`)
- survival_statistics.R (event counting)

**Recommendation:** Extract to a helper: `make_composite_id(df, id_col, treatment_col, cage_col)`.

### 3.5 No Formal Test Suite

Tests exist as ad-hoc scripts in `/temp/` but are not integrated into R's testing infrastructure. No `tests/testthat/` directory, no `tests/testthat.R` runner.

**Recommendation:** 
1. Create `tests/testthat/` with `usethis::use_testthat()`
2. Migrate validation script results into test assertions
3. Add to `.Rbuildignore`: `^temp$`

### 3.6 Composite ID Parsing Fragility

tumor_growth_statistics.R (~L570): Parses composite IDs using `strsplit(unique_id, "_")`. This breaks when treatment names contain underscores (e.g., "Drug_A"). There is a partial workaround that checks for known column values, but it's fragile and doesn't handle all cases.

**Recommendation:** Use a non-ambiguous separator (e.g., `"|||"`) or store components as separate columns rather than concatenating into a string.

---

## 4. Enhancement Opportunities

### 4.1 Additional Functionality

| Feature | Description | Priority |
|---------|-------------|----------|
| **Repeated-measures ANOVA** | Alternative to LME4 for simpler designs | Medium |
| **Body weight analysis** | Toxicity assessment from weight data | High (dashboard has placeholder tab) |
| **Tumor doubling time** | Common metric not currently computed | Low |
| **Multiple comparison methods** | Add Holm, FDR in addition to Bonferroni/Tukey | Medium |
| **Model diagnostics export** | Return formal diagnostic test results (Shapiro-Wilk, Levene's) | Medium |
| **Bayesian analysis option** | brms-based mixed models as alternative to lme4 | Low |
| **Formal simulation power** | Real Monte Carlo power analysis to replace fake "simulation" method | Medium |

### 4.2 API Improvements

1. **Consistent return structures:** LME4 path returns `pairwise_comparisons` as an emmeans object; AUC path returns `posthoc$pairwise` as a data frame. Downstream code must handle both cases. Standardize the return slot name and format.

2. **S3 class for results:** Define `class(result) <- "tgs_result"` with `print.tgs_result`, `summary.tgs_result`, `plot.tgs_result` methods.

3. **Formula interface:** Allow users to specify the model formula directly (e.g., `formula = Volume ~ Day * Treatment + (Day | ID)`) as an alternative to the parameter-based interface.

4. **Progress reporting:** Long-running functions (power analysis with simulations) should support progress callbacks, especially when called from Shiny.

---

## 5. File-by-File Summary

| File | Lines | Status | Key Issues |
|------|-------|--------|------------|
| tumor_growth_statistics.R | 975 | Bug | B1–B5, function too large, growth rate double-log |
| survival_statistics.R | 710 | OK | Well-structured adaptive method; could use refactoring |
| dose_response_statistics.R | 608 | OK | Deprecated `aes_string()`/`group_by_at()`; identical non-linear models |
| post_power_analysis.R | 1,218 | Warn | Extremely long; duplicated AUC calc; simulations are fake |
| analyze_drug_synergy.R | 363 | Warn | Loewe oversimplified; hardcoded combo name check |
| analyze_drug_synergy_over_time.R | 508 | OK | Good validation checks; deprecated `size` in plots |
| tumor_auc_analysis.R | 409 | OK | Clean; duplicated AUC logic |
| calculate_dates.R | 176 | OK | Fragile `in_place` pattern |
| calculate_volume.R | 130 | OK | Fragile `in_place` pattern; could use `pmax`/`pmin` |
| data.R | 111 | Good | Clean dataset documentation |
| plot_tumor_growth.R | 383 | Warn | Too long; contains extrapolation logic; deprecated `size` |
| plot_growth_rate.R | 226 | OK | Side-effect `cat()`; deprecated `size` |
| plot_auc.R | 216 | OK | Deprecated `size`; otherwise clean |
| plot_bliss.R | 107 | Warn | Hardcoded group names; deprecated `size`; legend issues |
| plot_caterpillar.R | 118 | Bug | `grep()` vs `grepl()` logic error; uses `@import ggplot2` |
| plot_combination_index.R | 125 | Warn | Constants in `aes()`; deprecated `size` |
| plot_treatments.R | 115 | Good | Cleanest plot file; modern `.data[[]]` usage |

---

## 6. Priority Action Items

### Critical (Fix Before Release)
1. **B2:** Fix double-log growth rate calculation
2. **B3:** Fix emmeans contrast ordering (use sorted group names)
3. **B1:** Use `log(Volume)` instead of `log(Volume + 1)`, with proper zero handling

### High Priority
4. **B4/B5:** Fix warning handler to not discard valid boundary-fit models; assign `best_model` in fallback
5. Consolidate 3 separate AUC implementations into one
6. Replace deprecated `aes_string()` and `group_by_at()` in dose_response_statistics.R
7. Fix `grep()` vs `grepl()` logic bug in plot_caterpillar.R
8. Add formal test suite (testthat)

### Medium Priority
9. Refactor `tumor_growth_statistics()` into sub-functions
10. Refactor `post_power_analysis()` into sub-functions
11. Replace deprecated `size` with `linewidth` across all plot files
12. Standardize return structure between LME4 and AUC paths
13. Convert full `import()` to specific `importFrom()` (especially dplyr/stats)
14. Move rarely-used dependencies to Suggests

### Low Priority
15. Deprecate `in_place` parameter
16. Add S3 class for result objects
17. Rename Loewe label to "Additive (mean)" or implement proper Loewe
18. Implement real Monte Carlo simulation for power analysis
19. Improve `.Rbuildignore` to exclude development files

---

*Review based on complete source code examination and independent validation against all 4 demo datasets.*
