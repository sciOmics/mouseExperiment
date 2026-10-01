# Statistical methods reference

For each analysis module: what the method does, when to choose it, what the result fields mean, and the assumptions you accept by using it.

This document complements the per-function roxygen — read it when you're trying to decide which model to use, not when you're looking up a specific argument. The roxygen is canonical for argument-level detail; this is canonical for the "why" and the "what's in the box".

For the Bayesian-specific diagnostics surface (Rhat, ESS, NUTS, LOO, Bayes R², PPC coverage, etc.), see [`BAYESIAN.md`](BAYESIAN.md).

---

## Endpoint TGI and evaluable days

Synergy, the therapeutic window, dose-response and the AUC all read each arm's volume at an endpoint from one model, the endpoint model (v0.26.0, `CODE_REVIEW.md` R20.1, R20.29):

- **Scale and curve:** log volume, with a natural spline in time (3 df) for each arm, or a straight line when an arm has fewer than 4 measured days.
- **Random effects:** per-animal random slopes, falling back to uncorrelated slopes and then to a random intercept when a fit fails.
- **Data:** fitted to every measurement of every animal, so animals removed at the volume limit still inform their arm. This is valid when removal depends on a measured volume (missing at random).
- **Zero volumes:** those measured before a tumour was palpable are left out; a zero after a positive volume is set to the smallest positive volume.
- **Reported mean:** each arm's geometric mean at the endpoint day, exp of its log-scale mean.
- **Intervals:** any function of the arm means (TGI, the Bliss excess, an AUC) gets its interval from draws of the model's fixed effects.

**Evaluable days** (`evaluable_days()`; maintainer decision R20-K, 2026-09-30). An arm is evaluable on a day when at least 50 % of its enrolled animals, and at least 3, are still on study (measured on that day or later). Every endpoint uses only days on which every arm it compares is evaluable. The default is the last such day, and a requested day that is not evaluable is an error. Past that day an arm's mean would be extrapolated beyond most of its own animals; the straight-line model that preceded this one put the Combo demo's control at 28,544 mm³ against about 3,000 observed. Each result reports `evaluability`: the rule, the day used, and the days and arms it excluded.

## Tumor Growth

Function family: `tumor_growth_statistics()`, `bayesian_tumor_growth()`, `tumor_doubling_time()`.

The tumor growth pipeline supports three models:

| Model | Function | Use when |
|---|---|---|
| `"lme4"` | `tumor_growth_statistics(..., model_type = "lme4")` | Default. Linear mixed-effects on log-volume with per-animal random slopes `(Day \| animal)` since v0.26.0, falling back to uncorrelated slopes and then a random intercept (with a warning) when a fit fails. Type III F-tests with Satterthwaite degrees of freedom. |
| `"auc"` | `tumor_growth_statistics(..., model_type = "auc")` | When total tumour burden over a window is the question. The area under each arm's fitted curve over the window in which every arm is evaluable, compared by AUC ratios with intervals from model draws (v0.26.0). |
| Bayesian LMM | `bayesian_tumor_growth(...)` | When you want posterior probability statements, credible intervals with direct probability interpretation, or you have small N where the prior matters. The first fit of a model structure in an R session compiles Stan code (about 20 s); later fits of the same structure reuse it. |

### Required columns

| Column | Notes |
|---|---|
| `id_column` | Identifies the animal within its arm and cage. IDs may restart in each arm or cage: every entry point groups animals by treatment + ID + cage (v0.25.0, `CODE_REVIEW.md` T1) |
| `time_column` | Numeric (day) or `Date` — date detection auto-converts |
| `volume_column` | Tumor volume (mm³). If you have `Length`/`Width`, run `calculate_volume()` first |
| `treatment_column` | Treatment group label |

Optional: `cage_column` (random intercept), `dose_column` (when crossed with treatment), `necrotic_column` (binary indicator → handled per `necrotic_handling`).

### Transforms (`transform` argument)

- `"log"` (default) — exponential growth is linear on log scale; standard.
- `"sqrt"` — variance-stabilizing for Poisson-like growth; rarely the right choice for tumors but available.
- `"none"` — for already-transformed data or when explicitly testing additive effects.

### LME4 path: what the result list contains

`tumor_growth_statistics(model_type = "lme4")` returns:

| Field | Description |
|---|---|
| `model` | The fitted `lmerMod` object (or `NULL` when `return_model = FALSE`) |
| `model_type_used` | `"lme4"` |
| `transform_used` | The transform that was applied |
| `treatment_effects` | One row per treatment group: `Adjusted_Mean`, `SE`, `Lower_CL`, `Upper_CL` (95% CI). Marginal means from `emmeans` at the mean study day |
| `pairwise_comparisons` | Treatment-pair contrasts. Columns: `contrast`, `estimate`, `SE`, `df`, `p.value`, `lower.CL`, `upper.CL` |
| `all_pairwise_comparisons` | All pairs, not Dunnett-filtered (used by the dashboard's forest plot) |
| `reference_comparisons_dunnett` | Only present when `reference_group` is specified; Dunnett-style multiplicity-adjusted contrasts vs reference |
| `anova`, `anova_method` | Type III F-tests with Satterthwaite degrees of freedom (`lmerTest`); `car`'s Wald χ² only if that fails |
| `random_effects` | The random-effects structure requested and used, and the reason for any fallback |
| `growth_rates` | Per-animal exponential growth rate (= slope on log scale); one row per animal |
| `data_summary` | N per group, mean baseline volume, mean final volume |
| `tumor_growth_plot` | ggplot of trajectories with group means overlaid (only when `plots = TRUE`) |
| `diag_qq_plot`, `diag_resid_fitted_plot` | Diagnostic plots (only when `include_diagnostics = TRUE`) |

### AUC path

Since v0.26.0 (`CODE_REVIEW.md` R20.4) the AUC is **model-based**. It is the area under each arm's fitted geometric-mean curve, from the endpoint model described under "Endpoint TGI and evaluable days" below. The window runs from the first study day to the last day on which every arm is evaluable (`auc_window`).
- `treatment_effects`: each arm's AUC with a 95 % interval.
- Pairwise comparisons: the AUC ratio and difference, with intervals and Wald p-values from 4,000 draws of the model's fixed effects, adjusted over the requested family.
- `anova`: a Wald test that every arm has the same log AUC.

`auc_analysis$individual` keeps each animal's trapezoidal AUC for description only. It covers the animal's own follow-up, so controls removed early at the volume limit get small AUCs; compared directly, that reversed the sign of the effect in 53.5 % of simulated studies.

### Assumptions you're accepting

- **LME4:** linearity in time on the transformed scale within each animal; animals differ in intercept and slope (random slopes); homoscedastic residuals. With a random intercept only, differing growth rates land in the residual and the Treatment × Day test rejected a true null in 73 % of simulated studies (R20.5).
- **AUC:** the endpoint model's assumptions (log-normal volumes, missing at random given the observed volumes, which holds when removal is triggered by a measured volume); independence between animals.
- **Bayesian LMM:** prior choice is honest; see [`BAYESIAN.md`](BAYESIAN.md)

### When to pick what

- **Default:** LME4. Familiar; testable; fast.
- **Small N (< 4 per group)**: Bayesian LMM with a `weakly_informative` or `informative` prior. Frequentist CIs at small N are notoriously narrow / over-confident; the Bayesian path is more honest.
- **Direct comparison of total tumor burden:** AUC. Use it as a confirmatory metric alongside the rate-based model, not as a replacement.

---

## Survival

Function family: `survival_statistics()`, `bayesian_survival()`.

### `survival_statistics()` — frequentist

Returns:

| Field | Description |
|---|---|
| `km_fit` | `survival::survfit` object (Kaplan-Meier) |
| `cox_model` | `survival::coxph` model |
| `log_rank_test` | `survival::survdiff` chi-squared test |
| `ph_test_table` | `cox.zph()` proportional-hazards check (one row per term + global). Reject the null (p < 0.05) means PH is violated for that term |
| `concordance` | C-index (Harrell's). 0.5 = random; 1.0 = perfect; report alongside HRs |
| `results` | Treatment-effect summary table |
| `survival_data` | The prepared data (Time, Event, Treatment) |

Note: when PH is violated, the Cox HR is no longer a constant time-ratio — it's an average over the time window. Switch to the Bayesian AFT path or stratify the Cox model in that case.

**Unit of randomisation (`randomisation_unit`, v0.25.0).** Say what treatment was assigned to:

- `"mouse"` (default): Cox, or Firth when an arm has no events, with ordinary model-based standard errors; cage is not modelled. When each cage holds one arm, the result carries a `cage_caveat`. The p-values treat cage-mates as independent; in simulation with a moderate cage effect that gave 16–19 % false positives at α = 0.05, and 2–4 % with none.
- `"cage"`: each arm is compared with the reference by a cage-level permutation log-rank that moves whole cages between the two arms, exactly when there are at most 5,000 assignments. `cage_permutation` records the assignments and the smallest attainable p-value: 1/3 with 2 cages per arm, 0.1 with 3 and 0.029 with 4. Hazard ratios are point estimates only, without an interval. Cages holding more than one treatment are refused, because treatment cannot then have been assigned to whole cages.

Before v0.25.0 the function added `cluster(cage)` whenever cages were replicated. A sandwich variance from 4–10 clusters gave 28 % false positives under the null (`CODE_REVIEW.md` R20.6).

### `bayesian_survival()` — Bayesian AFT

Two parametric families (exponential and gamma were removed in v0.23.0; exponential is Weibull with shape fixed at 1):

| Family | Use when |
|---|---|
| `"weibull"` (default) | Flexible — handles increasing, decreasing, and constant hazards. Reportable as both AFT (time ratio) and PH (hazard ratio). |
| `"lognormal"` | Hazard rises then falls. Biologically plausible for treated tumors where late-time animals tend to survive longer once they've cleared an initial dose |

Returns the standard `treatment_effects` shape (Group, Time_Ratio, Lower_CrI, Upper_CrI, HR, Median_Survival, Events, Total, Event_Rate, Note) plus the Bayesian diagnostics block (Rhat, ESS, NUTS, LOO).

It takes one row per animal, like `survival_statistics()`; a longitudinal frame is an error (v0.28.0; it used to be fitted with every measurement as a censored animal, `CODE_REVIEW.md` R20.39). The Weibull shape and log-normal sigma have weakly informative priors that do not follow `prior_strength` (R20.38; see [`BAYESIAN.md`](BAYESIAN.md)).

`include_cage_effect = TRUE` adds a cage-level frailty `(1 | cage)`. Worth using when you have ≥ 5 cages per group and reason to believe between-cage variation is non-trivial.

### Assumptions

- **Cox PH:** the hazard ratio is constant in time. Violated → switch to Bayesian AFT
- **Weibull AFT:** time on the log scale is linear in the covariates; the shape parameter is constant
- **All survival models:** non-informative censoring (animals censored for reasons unrelated to their hazard)

---

## Body Weight / Toxicity

Function family: `analyze_body_weight()`, `weight_loss_threshold()`, `therapeutic_window_metric()`.

Body-weight AUC, weight-corrected TGI, the efficacy–toxicity bivariate metric, total benefit area, and the Bayesian body-weight and therapeutic-window models were removed in v0.23.0: each had an estimand flaw that biased it toward calling a toxic arm safe or efficacious (`CODE_REVIEW.md` R20.9–R20.13, R20-K).

### `analyze_body_weight()` — frequentist LME4

Linear mixed-effects model with random intercept per animal (and optional cage random effect). Computes:

- Treatment-time interaction (does weight diverge over time?)
- Per-group adjusted means via `emmeans`
- Optionally `adjust_tumor_weight = TRUE`: subtracts estimated tumor weight (volume × `tumor_density`) before modelling. **`volume_units` (`"mm3"` or `"cm3"`) must then be declared** (v0.25.0, also for `weight_loss_threshold()` and `therapeutic_window_metric()`). The data are checked against the declared unit using the 90th percentile of volume, and a tumour mass above half the body weight stops the run. Units used to be inferred from the median volume, which read small-tumour mm³ studies as cm³ (`CODE_REVIEW.md` R20.83).
- **Weighings between calliper days** (v0.27.0, all three functions): animals are usually weighed more often than they are callipered. For the mass correction, each animal's volume is interpolated linearly between its calliper days, and held at the nearest measured value before the first and after the last. No weight row is dropped for lack of a same-day volume. Before, such rows were dropped, and weight-loss nadirs between calliper days were missed: in one scenario the mean worst loss read 12.4 % against a true 22.6 % (`CODE_REVIEW.md` R20.22).

Returns the same `treatment_effects` shape as `tumor_growth_statistics()` plus a `weight_trajectory_plot`.

### `weight_loss_threshold()` — time to a weight-loss threshold

Time to ≥ X % weight loss from each animal's own baseline, with Kaplan-Meier curves, a log-rank test and a Cox model per arm. When an arm has no events, standard Cox does not converge, and Firth's penalised Cox (`coxphf`) is used (`cox_method = "coxphf"`). Before v0.27.0 that fallback never ran, and a non-converged fit was reported as "cox" (HR 6 × 10⁹ against a log-rank p of 10⁻⁵; `CODE_REVIEW.md` R20.21).

**Why an animal's record ended** decides how it is counted. It can be read from an optional removal-reason column (`removal_reason_column`):

| Reason | Treated as |
|---|---|
| In `weight_loss_reasons` (e.g. "Body condition") | A weight-loss event at the animal's last day |
| In `planned_end_reasons` (end of study, scheduled sacrifice, data cut) | Censoring |
| Any other non-empty reason (tumour burden, death) | A competing removal |
| No reason column | Censoring at the last day |

With competing removals, `cuminc` gives the Aalen-Johansen cumulative incidence of weight loss, which `1 − KM` overstates. Without a reason column nothing marks a removal, so an early end is censoring, and `assumption` says so. Before v0.27.0 every record ending before the last study day was a competing removal, so staggered enrolment or a data cut halved the Aalen-Johansen incidence (`CODE_REVIEW.md` R20.23).

### Therapeutic window (two axes)

`therapeutic_window_metric()` reports two quantities per arm, each with a 95 % interval:

- **Efficacy:** the endpoint TGI described under "Endpoint TGI and evaluable days", at the last day on which every arm is evaluable. Its interval comes from draws of the endpoint model, or from a bootstrap of animals under the per-animal estimands.
- **Tolerability:** the worst weight loss, in percent of each animal's own first weighing, averaged over the arm's animals (`Worst_Loss`), with an interval from a bootstrap of animals. It is taken over each animal's whole record, not only up to the efficacy day. `N_Over_Threshold` counts the animals that reached the threshold on their own.

The tolerability flag compares the weight-loss interval with a declared threshold (`tolerability_threshold`, in percent, default 20):

| Flag | Rule |
|---|---|
| Tolerated | the upper bound is below the threshold |
| Not tolerated | the lower bound is at or above it |
| Unclear | the interval contains it |

The flag is `NA` without an interval. It describes the arm's average animal, so a "Tolerated" arm can still have an animal over the threshold.

Before v0.27.0 the function reported a ratio, TWM = TGI / weight loss, and ranked arms by it. It was replaced (`CODE_REVIEW.md` R20-K) because it needed an arbitrary floor for arms that lose no weight, set below the scales' own noise (R20.70), and because its two parts came from different estimators, so it had no matching interval (R20.2).

### Assumptions

- **BW LME4:** linearity on natural weight scale; per-mouse random effect captures individual differences
- **Adjusted weight:** the tumor density assumption (default 1.0 g/cm³) — most tumors are within 10% of this; large discrepancies matter
- **TGI:** the comparison is on the chosen endpoint day; trajectories with different shapes can have the same TGI at one day and differ on another
- **Therapeutic window:** the two axes come from separate analyses and are read side by side; no joint interval is implied
- **Weight-loss threshold:** without a removal-reason column, censoring at an early end is assumed non-informative

---

## Drug Synergy

Function family: `analyze_drug_synergy()`, `analyze_drug_synergy_over_time()`. (The Bayesian synergy models were removed in v0.23.0.)

### Bliss Independence

Tests whether the combo effect exceeds what's expected from independent action of each monotherapy:

```
expected_combo = effect_A + effect_B - effect_A × effect_B   (on fractional scale)
synergy        = observed_combo - expected_combo
synergy > 0    = supra-additive (Bliss synergy)
synergy < 0    = sub-additive (Bliss antagonism)
synergy ≈ 0    = independent
```

The fractional effects are the endpoint TGIs described under "Endpoint TGI and evaluable days", at the last day on which all four arms are evaluable. `synergy_ci` gives 95 % intervals for the reported estimates: from draws of the endpoint model under the default estimand, or from a bootstrap of animals under the per-animal estimands. The combination-versus-agent tests use the same estimand.

**The verdict** (`overall_assessment`, v0.27.0) comes from the 95 % interval for the Bliss excess:

| Verdict | Rule |
|---|---|
| Synergy | the interval lies above 0 |
| Antagonism | the interval lies below 0 |
| Additive (no departure from Bliss detected) | the interval contains 0 |

Without an interval (`n_boot = 0`) the point estimate is compared with symmetric bands, ± `additivity_margin` (default 0.1), and the verdict says it has no interval. When a single agent did not inhibit growth, Bliss does not apply and no verdict is given. The bands used to be asymmetric: any positive excess was "Synergy" but negative ones down to −0.1 were "Additivity", so a Bliss-additive combination was called synergistic in about half of simulated studies (`CODE_REVIEW.md` R20.19).

### Why there is no Combination Index

A Loewe-style CI was removed in v0.21.0. The formula above is the real Loewe
definition — dose-equivalence — and it needs a dose-response curve per agent so
the IC50s are known. Single-dose designs do not provide that, so what the package
actually computed was a stand-in:

```
CI = min(FE_A + FE_B, 1) / FE_combo
```

That is *response additivity*, a different null, and it fails the
sham-combination test: combine a drug with itself and it predicts twice the
fractional effect, so an agent at 50 % inhibition should reach 100 %. It will
not, so the method calls a drug antagonistic with itself.

Measured across the (FE_A, FE_B) grid with the combination set to **exactly
Bliss-additive**, 42 % of cells were labelled antagonistic and none synergistic —
a one-directional bias covering most of the range where active single agents sit.

If a Loewe analysis is required, collect per-agent dose-response curves and use
`drc::isobole()` or a full Loewe surface. It is a study-design change, not a code
change.

### `_over_time` variants

Same metric on every evaluable day, from one fit of the endpoint model, returning a `synergy_summary` table with one row per evaluable day and interval columns. Days on which an arm has thinned out are not analysed; `evaluability$excluded` lists them with the arms that failed the rule.

### Assumptions

- **Bliss:** the two drugs act independently (no shared targets / pathways). When the drugs hit the same pathway, "Bliss synergy" can look high but is mechanistically expected
- **Single-dose designs:** Bliss is the only synergy null this package computes, because it is the only one the data support. See "Why there is no Combination Index" above.

---

## Dose-Response

Function family: `dose_response_statistics()`. (The Bayesian dose-response model was removed in v0.23.0.)

### Frequentist (`drc` backend)

Fits a 4-parameter logistic (Hill / Emax) via `drc::drm`:

```
volume(dose) = lower + (upper - lower) / (1 + (dose / EC50)^slope)
```

Returns: `linear_model`, `anova_model`, `statistics$ec50`, `hill_slope`, `lower_limit`, `upper_limit`, `growth_dose_p_value`, plus the analysis data.

**One agent's dose series** (v0.27.0). The analysis takes one agent at several doses plus its control. The control is `control_group_name`, which must exist and have dose 0 (a missing dose is read as 0), or by default the rows with dose 0. The agent's arms are `treatments`. Without it, the data may hold only one other arm name, as when one treatment name is used at every dose. A dose held by more than one arm is always an error. Before, every arm went onto one dose axis: on the Master demo, Drug_B and the Drug_A + Drug_B combination were fitted as doses of Drug_A (`CODE_REVIEW.md` R20.16).

The analysis uses one day: by default the last day on which every dose group is evaluable (`endpoint_day`), with the animals measured that day. A requested day that is not evaluable, for example because the control has thinned out, is an error. Before v0.26.0 the default was each animal's own last observation, which gave removed controls capped, earlier volumes (EC50 7.5 against a true 2.2; R20.17). `tgi_table` reports TGI per dose group at `endpoint_day` from the endpoint model, with intervals, against the dose-0 control.

**The EC50** (v0.27.0) is `ED(model, 50)`, the dose giving half the fitted response range. Its 95 % interval is computed on log dose, so it stays positive. The lower asymptote is constrained to be at least 0. The 5-parameter (asymmetric) curve is considered only with at least 6 dose levels, and only if its lower asymptote is not negative. `ec50_in_range` is FALSE when the EC50 lies outside the tested doses, and `ec50_note` also flags a curve with as many parameters as dose levels, which nothing checks. The EC50 used to be the curve's `e` parameter, which is not the ED50 under the 5-parameter curve, and its interval was symmetric on the dose scale: [−394.7, 592.5] on the dashboard demo (`CODE_REVIEW.md` R20.18).

**Trend test.** The Jonckheere–Terpstra test is two-sided, and `direction` reports whether the dose-group means rise or fall. Choosing the one-sided alternative from the data had doubled its false-positive rate, to 0.096 (`CODE_REVIEW.md` R20.20).

**Growth rates.** Each animal's growth rate is the slope of log volume on day, from days with a positive measured volume; an animal needs three such days. Missing and zero volumes used to be set to half the animal's smallest volume, which turned the missing rows after a death into apparent shrinkage, and an animal with no positive volume stopped the whole analysis (`CODE_REVIEW.md` R20.14). `growth_rate_animals_left_out` counts the animals without enough days.

### Assumptions

- Smooth, monotonic dose-response (no biphasic / hormetic effects)
- Independence of measurements across animals
- The Hill equation is the right functional form (it usually is for cytotoxic response; check the diagnostic fit plot)

---

## Power Analysis

Function family: `apriori_power_analysis()`, `apriori_power_simulation()`. (`bayesian_power_analysis()` was removed in v0.23.0.)

### Analytic (`apriori_power_analysis()`)

Two-sample t-test power via `stats::power.t.test`, for one treated-vs-control comparison. With k ≥ 3 groups the study is analysed as k − 1 such comparisons, so k enters only through the per-comparison alpha: `alpha / n_comparisons` under the default Bonferroni correction (`p_adjust_method = "bonferroni"`, `n_comparisons = k − 1`). `dropout_rate` converts the analysable N into the number to enrol (`Enroll_N`). Fast (< 1 second). Use when:
- You have an effect-size estimate (Cohen's d, treated arm vs control) from a prior study
- You're sketching the study design — first-pass numbers

Before v0.24.0 the k ≥ 3 path powered a one-way ANOVA with f = d/√2 (CODE_REVIEW.md R20.7). The conversion for its own configuration is d/√(2k), and sample sizes came out 2.6–4× too small. Dunnett's exact correction (the Tumor Growth default) is slightly less conservative than Bonferroni, so the Bonferroni N slightly overestimates what a Dunnett analysis needs.

### LMM simulation (`apriori_power_simulation()`)

Simulates data under an LME4 model with specified variance structure, fits, and counts the proportion of simulated runs where the treatment effect is significant. The analysis powered is the default tumour-growth model, with random slopes (since v0.26.0; `CODE_REVIEW.md` R20.34). Slower (1-10 min depending on `n_sim`). Use when:
- Downstream analysis is the LMM path (you almost always want this for repeated-measures designs)
- You have an estimate of within-mouse variance from prior data
- You need to size for the right test (not the t-test approximation)

### Assumptions

- Effect size is on the same scale you'll analyze on (log volume for `transform = "log"`; raw volume for `"none"`)
- The simulated variance structure matches reality (the biggest assumption — sensitivity-analyze across plausible variance values)

---

## Utilities

| Function | Notes |
|---|---|
| `calculate_volume()` | Tumor volume from length × width. Default is the standard ellipsoid `V = (L × W²) × π / 6`. Other formulas: modified ellipsoid, cylinder, sphere. Dimensions are auto-corrected so the longer measurement is always L (handles "I measured them in the wrong order" data) |
| `calculate_dates()` | Date-to-day conversion using a reference date. Two methods: `"direct"` (study day already in data) and `"computed"` (compute from a date column) |
| `calculate_auc()` | Trapezoidal AUC over a time-volume vector |
| `tumor_doubling_time()` | Doubling time from exponential growth rate (returns `log(2) / rate`) |
| `tg_priors()`, `tg_mcmc()` | Config helpers (v0.4.7+) — bundle prior and MCMC arguments. See [`BAYESIAN.md`](BAYESIAN.md) |

---

## Choosing between frequentist and Bayesian

A short decision guide:

| Situation | Pick |
|---|---|
| N ≥ 6 per group, normal-ish residuals, simple comparison | Frequentist |
| N < 5 per group | Bayesian with `weakly_informative` or `informative` prior |
| You want "P(effect > 0)" or direct probability statements | Bayesian |
| You want to formally incorporate prior data | Bayesian |
| You need to publish in a venue that doesn't yet accept Bayesian methods | Frequentist (or both — Bayesian as a supplementary analysis) |
| You're short on compute (CI runs, exploratory analysis) | Frequentist |
| Frequentist CIs feel too narrow for the data you have | Bayesian (the discrepancy usually means you're under-N; the Bayesian prior makes uncertainty honest) |
| The frequentist model failed to converge | Bayesian (better-behaved when N is small) |

When in doubt, run both. The frequentist result tells you what a colleague reviewing the study will compute; the Bayesian result tells you what you actually believe given the data.

---

## See also

- [`BAYESIAN.md`](BAYESIAN.md) — diagnostics, priors, MCMC interpretation
- Roxygen `?<function_name>` — argument-level reference
- Vignettes in `vignettes/` — worked end-to-end examples
- Dashboard `docs/ARCHITECTURE.md` — how the UI exposes these methods
- `CODE_REVIEW.md` (repo root) — Round 2 audit; explains many design choices
