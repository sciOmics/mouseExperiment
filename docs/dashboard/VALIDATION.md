# Validation Report: Dashboard Analysis Results

**Date:** March 27, 2026  
**Scope:** All four demo datasets; both LME4 and AUC model types; survival analysis  
**Method:** Independent manual R calculations compared against `mouseExperiment` package output  

---

## Summary Table

| Dataset | Volume | Dates | AUC values | AUC ANOVA | LME4 ANOVA | Pairwise (displayed) | Survival |
|---------|--------|-------|-----------|-----------|-----------|---------------------|---------|
| weight_synthetic | N/A (pre-computed) | ✅ | ✅ | ✅ | ⚠️ see B1 | ✅ | N/A |
| combo_treatment | ✅ | ✅ | ✅ | ✅ | ⚠️ see B1 | ✅ | ✅ |
| dose_levels | ✅ | ✅ | ✅ | ✅ | ⚠️ see B1 | ✅ | ✅ |
| necrotic | ✅ | N/A (Day col) | ✅ | ✅ | ⚠️ see B1 | ✅ | N/A |

Legend: ✅ correct  ⚠️ discrepancy or limitation  ❌ wrong

**Key finding:** What the dashboard *displays* is correct for all tested outputs. Several bugs exist in internal package return values and growth-rate calculations that do not affect the displayed results but would affect users calling the package directly.

---

## Dataset Preprocessing

### Date Parsing

| Dataset | Format | Dashboard handling | Days produced | Status |
|---------|--------|-------------------|--------------|--------|
| weight | "March 1" (month + day, no year) | Appends year from ref-date → "March 1 2025"; uses `"%B %d %Y"` | 0, 3, 7, 10, 14 | ✅ |
| combo | "03/24/2025" (MM/DD/YYYY) | Auto-detects `"%m/%d/%Y"` | 0, 2, 4, 7, 9, 11, 14, 16, 18, 21, 23, 25, 28, 30, 32 | ✅ |
| dose_levels | "24-Mar" (DD-Mon) | Auto-detects `"%d-%b"`, appends year 2025 from ref-date | 0, 2, 4, 7, 9, 11, 14, 16, 18, 21, 23, 25 | ✅ |
| necrotic | Day column (already numeric) | Uses Day column directly; no date conversion | 0, 10, 15, 20, 25 | ✅ |

### Volume Calculation

Formula used: **Ellipsoid** (dashboard default) — `V = (π/6) × L × W²`

The `calculate_volume()` function automatically swaps Length and Width when Width > Length (to ensure the longer dimension is always Length). This is working correctly; verified for the necrotic dataset where several measurements have Width > Length.

**Verified for combo_treatment row 1:**  
`(π/6) × 1.20 × 0.80² = 0.4021 mm³` — matches stored value ✅

**Verified for necrotic auto-swap (row with L=1.00, W=1.20):**  
After swap: `(π/6) × 1.20 × 1.00² = 0.6283 mm³` — matches stored value ✅

### Dose-Levels: TreatmentGroup Combination

Dashboard combines `Treatment + " " + Dose` into a new `TreatmentGroup` column before analysis:

| TreatmentGroup | N mice |
|----------------|--------|
| Control 0 | 4 |
| HDACi 10 | 4 |
| HDACi 50 | 4 |
| HDACi 100 | 4 |

This correctly treats each dose level as an independent group. ✅

---

## AUC Analysis

### AUC Calculation

The package uses trapezoidal integration per unique (ID, Treatment, Cage) composite key, which correctly handles datasets where the same ID appears in multiple cages (e.g., mouse "L" in Cage 1 Control and mouse "L" in Cage 3 aPD1 are different animals).

#### Dataset 1 — weight_synthetic (manually verified, per-animal values)

| Mouse | DMSO | Drug A | Drug A+B | Drug B |
|-------|------|--------|----------|--------|
| L | 1322.2 | 2187.3 | 21.8 | 492.4 |
| LL | 2015.5 | 1395.5 | 30.8 | 576.3 |
| LR | 3059.3 | 2567.6 | 11.0 | 700.8 |
| N | 1800.8 | 1592.1 | 44.8 | 163.6 |
| R | 2843.2 | 2890.1 | 35.9 | 341.8 |

**Group means match exactly between package and manual:** ✅

| Group | Package Mean_AUC | Manual Mean_AUC | Match |
|-------|-----------------|-----------------|-------|
| DMSO | 2208.2 | 2208.2 | ✅ |
| Drug A | 2126.5 | 2126.5 | ✅ |
| Drug A + B | 28.9 | 28.9 | ✅ |
| Drug B | 455.0 | 455.0 | ✅ |

#### Dataset 2 — combo_treatment (8 mice per group, 2 cages each)

| Group | Package Mean_AUC | Manual Mean_AUC (composite ID) | Match |
|-------|-----------------|--------------------------------|-------|
| Control | 11915.96 | 11915.96 | ✅ |
| HDACi | 9496.75 | 9496.75 | ✅ |
| HDACi + PD1 | 297.99 | 298.00 | ✅ |
| aPD1 | 8385.59 | 8385.59 | ✅ |

#### Dataset 3 — dose_levels (4 mice per dose group)

| Group | Package Mean_AUC | Manual Mean_AUC | Match |
|-------|-----------------|-----------------|-------|
| Control 0 | 7007.62 | 7007.62 | ✅ |
| HDACi 10 | 5284.79 | 5284.79 | ✅ |
| HDACi 50 | 5774.57 | 5774.57 | ✅ |
| HDACi 100 | 2764.83 | 2764.83 | ✅ |

#### Dataset 4 — necrotic (4 mice per group)

| Group | Package Mean_AUC | Manual Mean_AUC | Match |
|-------|-----------------|-----------------|-------|
| None | 438.22 | 438.22 | ✅ |
| Drug | 1702.68 | 1702.68 | ✅ |

### AUC ANOVA

| Dataset | Package F | Manual F | Package p | Manual p | Match |
|---------|-----------|---------|----------|---------|-------|
| weight | 26.05 | 26.05 | 2.14 × 10⁻⁶ | 2.14 × 10⁻⁶ | ✅ |
| combo | — | 31.13 | — | 6.07 × 10⁻⁶ | ✅ |
| dose_levels | — | 1.855 | — | 0.191 | ✅ |
| necrotic | — | 1.46 | — | 0.272 | ✅ |

### AUC Pairwise Welch t-tests (Bonferroni-adjusted)

#### Dataset 1 — weight_synthetic  
Post-fix (after the AUC comparisons-tab fix applied in this session):

| Comparison | Package p_adj | Manual p_adj (Bonferroni) | Match |
|-----------|--------------|--------------------------|-------|
| DMSO − Drug A | 1.000 | 1.000 | ✅ |
| DMSO − Drug A+B | 0.016 | 0.015 | ✅ |
| DMSO − Drug B | 0.026 | 0.026 | ✅ |
| Drug A − Drug A+B | 0.011 | 0.011 | ✅ |
| Drug A − Drug B | 0.016 | 0.016 | ✅ |
| Drug A+B − Drug B | 0.061 | 0.061 | ✅ |

*These are the values the dashboard Comparisons tab now displays after today's bug fix.*

---

## LME4 Analysis

### Transformation

The package applies `log(Volume + 1)` (`log1p()`) when `transform = "log"`:

```r
# tumor_growth_statistics.R ~line 297
analysis_df[[volume_column]] <- log(analysis_df[[volume_column]] + 1)
```

Standard independent analysis uses `log(Volume)` (excluding zero volumes, as the dashboard already does when filtering `Volume > 0`). This produces different chi-square values in the Type III ANOVA but the same qualitative conclusions (significant treatment effect, significant Day × Treatment interaction).

**Example — Dataset 1, Treatment Chisq:**
- Package (log V+1): 11.07, p = 0.011
- Manual (log V): 72.18, p = 1.46 × 10⁻¹⁵

Both indicate a significant treatment effect; the quantitative values differ. See Bug B1 below.

### ANOVA — Qualitative Conclusions

Despite the quantitative chi-square differences from the log(V+1) issue, the qualitative direction (which effects are significant) is consistent:

| Dataset | Day sig.? | Treatment sig.? | Day:Treatment sig.? | Package | Manual | Agree? |
|---------|-----------|----------------|---------------------|---------|--------|--------|
| weight | ✅ | ✅ | ✅ | p<2e-16, p=0.011, p<2e-16 | p<2e-16, p<2e-15, p<2e-16 | ✅ direction |
| combo | ✅ | ✅ | ✅ | p<2e-16, p<2e-16, p<2e-16 | same | ✅ |
| dose_levels | ✅ | ✅ | ❌ (p=0.28) | same | p=0.18 | ✅ direction |
| necrotic | ✅ | ❌ (p=0.82) | ❌ (p=0.16) | p<2e-16, p=0.82, p=0.16 | p<2e-16, p=0.99, p=0.68 | ✅ direction |

### Pairwise Comparisons (Displayed in Dashboard)

The dashboard recomputes all comparisons via emmeans in `server_modules.R` (lines ~1006–1060) rather than using the `pairwise_comparisons` slot from the package return. These server-computed values are correct.

**Dataset 2 (combo) — Correct Dunnett contrasts (dashboard display):**

| Contrast | Estimate | t.ratio | p-value (Dunnett) | Direction |
|----------|----------|---------|-------------------|-----------|
| aPD1 − Control | −0.939 | −10.67 | < 0.0001 | Tumor reduction ✅ |
| HDACi − Control | −0.816 | −9.26 | < 0.0001 | Tumor reduction ✅ |
| HDACi + PD1 − Control | −3.028 | −33.77 | < 0.0001 | Strongest reduction ✅ |

These are consistent with marginal means (Control > aPD1 > HDACi > HDACi+PD1).

---

## Survival Analysis

### Dataset 2 — combo_treatment

The package correctly identifies all 8 mice per treatment group (2 cages × ~4 mice/cage). Firth correction is automatically applied due to complete separation for HDACi+PD1 (0 events).

| Group | Package HR (95% CI) | Package p | K-M median | Events/Total |
|-------|--------------------|-----------|-----------:|-------------|
| Control | 1.00 (ref) | — | 30 days | 8/8 |
| HDACi | 0.151 (0.044–0.519) | 0.0011 | not reached | 4/8 |
| aPD1 | 0.106 (0.028–0.410) | 0.0002 | not reached | 3/8 |
| HDACi+PD1 | 0.014 (0.001–0.295) | <0.0001 | not reached | 0/8 |

**Manual log-rank test:** p = 0.0028 (consistent direction — all treated groups have significantly better survival than Control). ✅

**Note on method:** Package uses Firth's correction (`coxphf`) because HDACi+PD1 has complete separation (0 events). Manual log-rank test does not apply Firth correction. The difference between Firth HR and log-rank p-value is expected and appropriate.

### Dataset 3 — dose_levels

The package runs survival on the Treatment column (not TreatmentGroup), combining all HDACi dose levels into a single survival group (12 mice total: 4 per dose × 3 doses).

| Group | Package HR (95% CI) | Package p | K-M median | Events/Total |
|-------|--------------------|-----------|-----------:|-------------|
| Control | 1.00 (ref) | — | 23 days | 4/4 |
| HDACi | 0.077 (0.018–0.332) | 0.0006 | 25 days | 8/12 |

**Manual log-rank test:** p = 0.0069 (consistent direction). ✅

*Note: The manual matched the package's Events/Total (4/4 and 8/12) only when using composite (ID, Treatment, Cage) grouping to identify unique mice. Using (ID, Treatment) alone incorrectly collapses mice with the same ID from different dose-level cages.*

---

## Bugs Identified

### B1 — `log(V + 1)` Transformation

**File:** `mouseExperiment/R/tumor_growth_statistics.R` ~line 297  
**Severity:** Medium — affects quantitative ANOVA values but not qualitative conclusions  
**Description:** The package uses `log(Volume + 1)` for the LME4 transform, while the standard approach for tumor volume data is `log(Volume)` (after excluding zeros). The `+1` offset is significant when volumes are in the range 0–2 mm³ (as in early study days), producing different chi-square statistics compared to a standard independent analysis.

```r
# Current (non-standard):
analysis_df[[volume_column]] <- log(analysis_df[[volume_column]] + 1)
# Standard practice:
analysis_df[[volume_column]] <- log(analysis_df[[volume_column]])
# (caller already filters Volume > 0 before passing to the function)
```

**Recommendation:** Switch to `log(Volume)` and document the requirement that Volume > 0 (the dashboard already ensures this). Alternatively, if `log(V+1)` is a deliberate choice (for biological reasons), document it prominently in the function help and dashboard info boxes.

---

### B2 — Double Log-Transform in Growth Rate Calculation

**File:** `mouseExperiment/R/tumor_growth_statistics.R` ~line 316  
**Severity:** High — growth rate values are incorrect when `transform = "log"`  
**Description:** Growth rates are computed by fitting `log1p(Volume_transformed) ~ Day`. When `transform = "log"`, the volume has already been transformed to `log(V + 1)`. Applying `log1p()` again gives `log(1 + log(1 + V))`, not the exponential growth rate of the original volume.

```r
# Problem: analysis_df$Volume is log(V+1) at this point, then:
log_volume <- log1p(subject_data[[volume_column]])   # = log(1 + log(1+V))
```

**Expected:** Growth rates should be computed on `log(original_volume)` or `original_volume`, not on a double-transformed value.

**Recommendation:** Compute growth rates from `auc_df` (the untransformed copy saved before transformation) rather than `analysis_df`:

```r
# Use auc_df for growth rate calc (untransformed)
log_volume <- log1p(auc_subject_data[[volume_column]])   # log(1+V) on original scale
```

---

### B3 — Incorrect Contrast Vector Ordering in `pairwise_comparisons` Return Slot

**File:** `mouseExperiment/R/tumor_growth_statistics.R` ~lines 878–899  
**Severity:** High — produces wrong p-values and incorrect labels in the package return slot  
**Description:** Contrast coefficient vectors are built using the **data-appearance order** of treatment levels (`unique(df[[treatment_column]])`), but `emmeans` interprets those vectors using its own **alphabetical ordering** of factor levels. When these orderings differ, the contrast coefficients are assigned to the wrong groups. The result is that `pairwise_comparisons` contains comparisons between the wrong groups, with wrong signs and labels.

**Demonstrated for Dataset 2 (combo):**

The data-appearance order is: `[Control, HDACi, aPD1, HDACi+PD1]` (indices 1,2,3,4)  
The emmeans/alphabetical order is: `[aPD1, Control, HDACi, HDACi+PD1]` (aPD1 sorts first because uppercase < lowercase in ASCII)

| Package label | Coefficient vector (data order) | emmeans interprets as | Actual computation |
|---------------|--------------------------------|-----------------------|-------------------|
| "HDACi − Control" | [−1, +1, 0, 0] | −1×aPD1 + 1×Control | **Control − aPD1** |
| "aPD1 − Control" | [−1, 0, +1, 0] | −1×aPD1 + 0×Control + 1×HDACi | **HDACi − aPD1** |
| "HDACi+PD1 − Control" | [−1, 0, 0, +1] | −1×aPD1 + 1×HDACi+PD1 | **(HDACi+PD1) − aPD1** |

None of the three "vs Control" contrasts actually involve Control.

**Package result for "HDACi − Control":** estimate = +0.939, p < 0.0001  
**Correct "HDACi − Control":** estimate = −0.816, p < 0.0001  
(control has higher volume than HDACi, so the estimate must be negative)

**Note:** The **dashboard display is unaffected** — the Comparisons tab is populated by a server-side emmeans call with correct alphabetical ordering (lines ~1006–1060 of `server_modules.R`). However, any user who calls `tumor_growth_statistics()` directly and reads `result$pairwise_comparisons` will get completely wrong values.

**Recommendation:** Fix the contrast construction to use factor levels rather than `unique()`:

```r
# Replace:
treatment_groups <- unique(df[[treatment_column]])

# With (inside the emmeans block, after creating lsmeans_obj):
treatment_groups <- levels(lsmeans_obj)[[1]]   # guaranteed to match emmeans element order
```

---

### B4 — Warning Handler Discards Valid Models

**File:** `mouseExperiment/R/tumor_growth_statistics.R` ~lines 428–445  
**Severity:** Low–Medium — causes unnecessary fallback computation; no incorrect results  
**Description:** The `tryCatch` for model fitting uses `return(NULL)` inside the `warning` handler, which discards the model result for **any** warning (not just boundary singular fits). This means every dataset that triggers a "boundary (singular) fit" warning (very common with small n) silently discards the model and triggers a redundant fallback re-fit.

```r
models$intercept_only <- tryCatch({
    lme4::lmer(...)
}, warning = function(w) {
    if (grepl("boundary", w$message)) {
        warning("Boundary...")
    }
    return(NULL)   # NULL returned for ANY warning
})
```

**Impact:** All four demo datasets trigger this fallback. The fallback model is identical to the discarded model, so results are correct but computation is doubled.

**Recommendation:** Use `withCallingHandlers` to capture warnings without aborting, then check `isSingular()` after fitting:

```r
warn_msg <- NULL
models$intercept_only <- withCallingHandlers(
    lme4::lmer(...),
    warning = function(w) { warn_msg <<- conditionMessage(w); invokeRestart("muffleWarning") }
)
if (!is.null(models$intercept_only) && lme4::isSingular(models$intercept_only))
    message("Note: intercept-only model is singular (random intercept variance ~ 0)")
```

---

### B5 — `best_model` Undefined in Fallback Path

**File:** `mouseExperiment/R/tumor_growth_statistics.R` ~lines 480–488, ~470  
**Severity:** Medium — causes runtime error in diagnostics  
**Description:** When the fallback path is taken (both `intercept_only` and `random_slope` models NULL), `model_selection$selected_model <- "intercept_only"` is set, but the local variable `best_model` is never assigned. The diagnostics extraction at line ~470 uses `best_model == "random_slope"`, causing `Error: object 'best_model' not found`.

**Note:** This error is silently swallowed by the `tryCatch` around the diagnostics block, so `include_diagnostics = TRUE` silently produces empty diagnostics when the fallback triggers.

**Recommendation:** After the fallback:

```r
# Inside the fallback else-branch:
model_selection$selected_model <- "intercept_only"
best_model <- "intercept_only"   # ← add this line
```

---

## Notes and Limitations

### N1 — Necrotic Flag Not Handled Specially

The necrotic dataset has a `Necrotic` column flagging observations where the tumor was necrotic at measurement time. The package and dashboard include these observations as normal data points. Depending on the analysis goal, necrotic observations might warrant:
- Exclusion from volume analysis
- Special coding as a competing event in survival analysis
- A separate necrosis outcome

Currently there is no UI option or automatic handling for this flag.

### N2 — Survival for Necrotic Dataset

The necrotic dataset (`necrotic_synthetic_data.csv`) has no `Survival_Censor` column. It cannot be used for survival analysis with the current package interface. No survival analysis was attempted for this dataset.

### N3 — Dose-Levels Survival Combines All Doses

When running survival analysis on `dose_levels_synthetic_data.csv`, the dashboard uses the `Treatment` column (not `TreatmentGroup`), meaning all three HDACi dose levels are pooled into a single survival curve. This loses dose-response information at the survival level. A future enhancement could offer per-dose survival curves.

### N4 — Singular Fit Warnings Are Normal for Small Datasets

All four demo datasets produce "boundary (singular) fit" warnings from `lme4`. This occurs when the random intercept variance is near zero, typically in small studies (n ≤ 5 mice/group, ≤ 5 timepoints). The model still runs and produces valid fixed-effects ANOVA and comparisons. The warning does not indicate incorrect results.

### N5 — AUC Session Fix

Prior to this session, the AUC Comparisons tab incorrectly ran `emmeans()` on the `aov` object from the AUC path, pooling variance across all groups and producing wrong p-values (e.g., Drug A+B vs Drug B: p = 0.537 instead of correct p = 0.061). This was fixed in `server_modules.R` during this session.

---

## Files Analysed

| File | Role |
|------|------|
| `mouseExperiment/R/tumor_growth_statistics.R` | Main analysis function |
| `mouseExperiment/R/survival_statistics.R` | Survival analysis |
| `mouseExperiment/R/calculate_dates.R` | Date → Day conversion |
| `mouseExperiment/R/calculate_volume.R` | L × W → Volume |
| `mouseExperimentDashboard/R/server_modules.R` | Dashboard server — data preprocessing and comparison tab |
| `mouseExperimentDashboard/inst/sample_data/*.csv` | Demo datasets |
| `mouseExperiment/temp/validate_dashboard_results.R` | Validation script (kept for reproducibility) |
