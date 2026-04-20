# Toxicity & Efficacy-Toxicity Analysis — Implementation Plan

## Phase 1: Infrastructure / Prerequisites

### 1.1 Fix weight column sync bug
`col_weight` is configured in the Data Upload UI but **never synced to `rv$settings`** in `mod_data_upload.R`. Must add `col_weight = input$col_weight` to the settings observer, and add weight pattern detection to `detect_columns()` in `helpers.R`.

### 1.2 Add Sex column as optional input
No sex/gender column currently exists. Add `col_sex` to the Data Upload UI (`ui_modules.R`), sync to `rv$settings`, and add auto-detection for common names (`Sex`, `Gender`, `M/F`).

### 1.3 Resolve `ensure_preprocessed()` for new columns
Extend to resolve `rv$weight_col` and `rv$sex_col` from `rv$settings`, same pattern used for dose/cage/id.

---

## Phase 2: Package Functions (`mouseExperiment`)

### 2.1 `analyze_body_weight()` — Mixed Modeling
- Compute estimated tumor weight: w_tumor = V_mm³ × ρ / 1000 where ρ defaults to 1.0 g/cm³
- Net mass: m_net = m_body − w_tumor (toggleable)
- Fit via `lme4::lmer()` with REML (default) or ML:
  ```
  net_mass ~ Treatment * Day + [optional: Sex + Initial_Mass + Volume] + (1 + Day | MouseID)
  ```
- Return: model summary, fixed effects table, random effects variance, marginal means via `emmeans`, diagnostic plots (residuals, Q-Q)
- Handle convergence failures gracefully (singular fits, boundary warnings)

### 2.2 `body_weight_auc()` — AUC Analysis
- Trapezoidal AUC of net mass (or raw mass) per mouse, similar to existing `tumor_auc_analysis()`
- Normalize as % change from baseline: (m_t − m_0) / m_0 × 100
- Return: per-mouse AUC, group summaries (mean, SD, SEM), pairwise comparisons (t-tests or Wilcoxon)

### 2.3 `weight_loss_threshold()` — Time-to-Threshold
- Per mouse: find first day where m_t ≤ m_0 × (1 − threshold)
- Default threshold: 20%, user-adjustable
- Baseline day: default = first measurement, user-adjustable (with note that first treatment day is recommended)
- Right-censored if threshold never reached
- Fit Kaplan-Meier + log-rank test via `survival::survfit()` / `survival::survdiff()`
- Optionally fit Cox PH model if ≥2 treatment groups

### 2.4 `therapeutic_window_metric()` — TWM
- TWM = TGI / MaxWeightLoss% per treatment group
- Higher TWM = better therapeutic window
- Return: per-group TWM, ranking, confidence intervals via bootstrap

### 2.5 `efficacy_toxicity_bivariate()` — Safety-Efficacy Data
- Compute per-group: toxicity metric (max % weight loss or AUC of weight change) + efficacy metric (selectable: Final TGI, Tumor AUC, or Log-Cell Kill)
- Return: data frame for plotting + correlation statistics

### 2.6 `total_benefit_area()` — Integrated Efficacy-Toxicity
- B = AUC_efficacy − λ × AUC_toxicity
- λ = user-adjustable trade-off parameter (default 1.0)
- Return: per-group benefit score, ranking

### 2.7 `weight_corrected_tgi()` — Weight-Corrected TGI
- Filter to mice whose weight never dropped below safety threshold
- Recompute TGI on this subset only
- Return: corrected TGI, n excluded per group, comparison to uncorrected TGI

---

## Phase 3: Dashboard Module (`mouseExperimentDashboard`)

### 3.1 `toxicityUI()` / `toxicity_server()` — new `R/toxicity_modules.R`

**Sidebar parameters:**
- Reference/control group selector
- Tumor weight adjustment toggle (default: on) + density input (default: 1.0 g/cm³)
- Covariates checkboxes: Tumor Volume (default: on), Sex, Initial Mass
- Estimation method: REML / ML radio
- Weight loss threshold % slider (default: 20%)
- Baseline day selector

**Result tabs:**
1. **Body Weight Mixed Model** — summary text, fixed effects table, forest plot of treatment effects
2. **Body Weight AUC** — per-group AUC bar chart, downloadable table
3. **Weight Loss Survival** — Kaplan-Meier curve (time to threshold), risk table
4. **Therapeutic Window** — TWM bar chart by treatment group
5. **Safety-Efficacy Plot** — bivariate scatter (X = toxicity, Y = selectable efficacy metric), with quadrant shading
6. **Total Benefit** — integrated benefit score bar chart with λ slider
7. **Weight-Corrected TGI** — side-by-side comparison of standard vs corrected TGI

### 3.2 Enable the toxicity tab
In `app.R` — remove the disabled styling and redirect guard, wire up `toxicity_server()`.

---

## Open Questions

**Q1 — Log-Cell Kill metric:** This requires knowledge of tumor doubling time. The package has `tumor_doubling_time()` — should Log-Cell Kill be computed as LCK = log₂(V_control / V_treated) / T_doubling, or is there a preferred formula?

Preferred Formula:$$LCK = \frac{T - C}{3.32 \times T_{d}}$$$T$: Time (days) for the treated group to reach a target volume.$C$: Time (days) for the control group to reach the same target volume.$T_{d}$: Control group doubling time.3.32: Constant to convert $\log_{10}$ to doubling (since $\log_{10}(2) \approx 0.301$ and $1/0.301 \approx 3.32$).Recommendation: Implement the "Growth Delay" version above, as it is the most common in the literature for longitudinal mouse data.

**Q2 — TWM denominator:** If a treatment group has zero or negligible weight loss, TWM → ∞. Should we cap it, use a small floor (e.g., max(weight_loss, 1%)), or flag those groups specially?

Present it as a "Safety Score" where if weight loss is $\leq$ baseline noise, the score is simply set to the max inhibition value ($TGI \times 1$).

**Q3 — Mixed model complexity vs convergence:** With small group sizes (typical mouse studies: 5–10/group), the full model `Treatment * Day + Sex + Initial_Mass + Volume + (1 + Day | ID)` will likely hit convergence issues. Suggestion:
- Default to random intercept only `(1 | ID)` — add random slope as a checkbox that warns about small sample requirements
- Auto-fall back: try full model → if singular, reduce to intercept-only and warn

Recommendation: Auto-fallback is the most user-friendly path.
Attempt: (1 + Day | ID) (Random Slope + Intercept).
If failed/singular: Fall back to (1 | ID) (Random Intercept only).
UI Feedback: Display a small "Model simplified due to sample size" badge next to the results.

**Q4 — Bivariate plot granularity:** Each point represents a treatment group mean. With only 3–5 groups in a typical study, the plot will be sparse. Should we also offer a per-mouse version (each point = one mouse) with group coloring?

Recommendation: Default to the Per-Mouse version.Mean-only plots with $N=4$ groups look like a mistake to most users.Plotting individual mice (points) with a translucent "ellipse" or centroid for the group mean provides much better insight into the drug's consistency. It also helps identify if one "super-responder" is skewing a group's average.

**Q5 — Tab organization:** Two options:
- **A)** Single "Toxicity" tab with sub-tabs for all 7 result panels
- **B)** Split into two top-level tabs: "Toxicity" (mixed model, AUC, weight loss survival) and "Efficacy-Toxicity" (TWM, bivariate, total benefit, corrected TGI)

Option B keeps each tab focused but adds more top-level navigation.

Recommendation: Option B

**Q6 — Weight data without tumor volume:** Should the toxicity module work even if no tumor volume data is available (skipping tumor weight adjustment and tumor-dependent covariates), or require both?

Recommendation: Make it flexible.
The module should work without tumor data. Many researchers run "Tox-only" pilot studies before they even transplant tumors to find the MTD (Maximum Tolerated Dose). If Volume columns are missing, gray out the "Tumor Adjustment" and "Efficacy-Toxicity" features, but keep the core weight loss Mixed Model and Time-to-Threshold active.

---

## Suggested Improvements

1. **Body weight trajectory plot** — a spaghetti / mean±SEM line plot of body weight over time by treatment group (similar to tumor growth curves). This is the most commonly used toxicity visualization and would be expected by users before any statistical output.

Recommendation: Yes, implement this. 

2. **Nadir weight analysis** — report each mouse's minimum weight and the day it occurred. Useful clinically and computationally cheap.

Recommendation: Yes, implement this. 

3. **Correlation heatmap** — for the bivariate analysis, if multiple efficacy metrics are available, show a correlation matrix of efficacy vs toxicity metrics across groups.

Recommendation: Yes, implement this. 

4. **Export combined report** — a single downloadable summary combining key toxicity + efficacy metrics per group in one table.

Recommendation: Yes, implement this. 
