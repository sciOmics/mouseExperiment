# Bayesian guide

How to use the package's Bayesian models, interpret the diagnostics, and choose priors. This document is focused — it assumes you've already decided to use a Bayesian model (see [`METHODS.md`](METHODS.md) for that decision) and now need to do it well.

## Quick reference

```r
fit <- bayesian_tumor_growth(
  df     = df,
  priors = tg_priors(strength = "weakly_informative"),
  mcmc   = tg_mcmc(chains = 4, warmup = 1000, iter = 500, seed = 42))

# Before interpreting results, check diagnostics:
print(fit$mcmc_diagnostics)        # Rhat, ESS_Bulk, ESS_Tail, Converged
print(fit$nuts_diagnostics)        # divergences, max_treedepth, E-BFMI
print(fit$bayes_R2)                # Bayes R² with 95% CrI
print(fit$loo_diagnostics)         # elpd_loo + n_high_k (Pareto-k above k_threshold)
fit$pp_check_plot                  # posterior predictive density overlay
fit$prior_posterior_plot           # how much the data updated the prior

# Then look at effects:
print(fit$treatment_effects)       # adjusted means with 95% CrI
print(fit$pairwise_comparisons)    # contrasts with posterior P direction
```

If any diagnostic is bad, the effect estimates are unreliable. Read the next section before reporting anything.

---

## The Bayesian diagnostics checklist

Run through this list every time. The package surfaces each piece via a result field; the dashboard's MCMC Diagnostics tab renders them via `helpers_bayes.R::bayes_diagnostics_panel()`.

### 1. Rhat (potential scale reduction factor)

**Threshold:** all parameters Rhat ≤ 1.01. Older guidance said 1.1 — that's too lax.

**What it measures:** convergence across chains. Rhat = 1.0 means chains have converged to the same distribution; > 1.01 means at least one chain hasn't.

**Found in:** `result$mcmc_diagnostics$Rhat` (one row per fixed-effect parameter) + the `Converged` flag.

**If it's bad:**
- Increase `n_iter` (start: 1000 post-warmup)
- Increase `n_warmup` (start: 2000)
- Check `nuts_diagnostics` — divergences usually drive high Rhat
- Reparameterize if the issue is identifiability (e.g., centered → non-centered random effects)

### 2. ESS (effective sample size)

**Threshold:** Bulk_ESS ≥ 400 and Tail_ESS ≥ 400 for every parameter you care about. ≥ 1000 is comfortable.

**What it measures:** how many independent draws your MCMC effectively produced. Lower = more autocorrelation = noisier estimates.

**Found in:** `result$mcmc_diagnostics$Bulk_ESS` and `Tail_ESS`.

**If it's bad:** more iterations. The fix is usually mechanical.

### 3. NUTS sampler diagnostics

**Threshold:** zero divergent transitions; zero max_treedepth hits; E-BFMI ≥ 0.3.

**What they measure:**
- **Divergences** — the sampler couldn't follow the posterior's geometry. Even one divergence in a small / moderate fit usually means the posterior estimates are biased.
- **Max treedepth** — the sampler's adaptive step ran out of room. Less catastrophic than divergences but indicates the posterior is hard to sample.
- **E-BFMI** — energy stability across iterations. Below 0.3 indicates the sampler isn't exploring the posterior efficiently.

**Found in:** `result$nuts_diagnostics` (one row per chain: `chain_id`, `n_divergent`, `n_max_treedepth`, `ebfmi`, etc.).

**If divergences are present:**
- Raise `adapt_delta` (pass to `brms::brm()` directly via `control = list(adapt_delta = 0.95)`; default is 0.8)
- For random-effects models, switch to a non-centered parameterization
- Tighten the prior (a too-diffuse prior on variance components is the most common cause)
- As a last resort, reduce model complexity (drop a random slope, etc.)

**Dashboard surfaces this:** the "Sampler diagnostic warnings" alert at the top of the MCMC Diagnostics tab fires when any of these flags are non-zero.

### 4. Bayes R²

**Threshold:** there isn't one, but report it always. A Bayes R² of 0.7 with a tight CrI means the model explains 70% of variance reliably; 0.7 with a wide CrI (e.g., 0.4 to 0.9) means you don't know how well it fits.

**What it measures:** posterior variance-explained ratio, model-side. Direct analog of frequentist R² but with credible-interval uncertainty.

**Found in:** `result$bayes_R2` (single-row data frame: `Estimate`, `Lower_95_CrI`, `Upper_95_CrI`).

### 5. Posterior predictive coverage

**Threshold:** empirical coverage should match nominal. 50% intervals should contain ~50% of observed data; 80% → ~80%; 95% → ~95%. Off by 5 percentage points is acceptable; more than that suggests model misspecification.

**What it measures:** does the posterior predictive distribution actually cover the observed data at the nominal rate?

**Found in:** `result$ppc_coverage` (single-row data frame: `cov_50`, `cov_80`, `cov_95`, `n_obs`). For survival models only uncensored times count, because a censored time is a lower bound rather than an observation; the predictive draws come from the fit's seed, not the session's random-number stream (v0.28.0).

**If it's off:** the model's predictive distribution is too narrow (under-coverage) or too wide (over-coverage). Check the residuals; consider a richer likelihood (e.g., Student-t instead of Normal if outliers are an issue).

### 6. PSIS-LOO + Pareto-k

**Threshold:** all Pareto-k below `k_threshold`, which depends on the number of posterior draws S: min(1 − 1/log10 S, 0.7), so 0.667 at S = 1,000 and 0.64 at S = 600 (the loo 2.6 rule). Before v0.28.0 the threshold was a fixed 0.7, and any loo warning, which is exactly what an influential observation raises, discarded the whole LOO result, so `n_high_k` could only ever read 0 (`CODE_REVIEW.md` R20.35).

**What it measures:** leave-one-out cross-validation via Pareto-smoothed importance sampling. Pareto-k flags observations the LOO approximation can't reliably estimate — usually highly influential observations.

**Found in:** `result$loo_diagnostics` (`elpd_loo`, `se_elpd`, `k_threshold`, `n_high_k`, `warnings`, plus the list-column `pareto_k` with the per-observation Pareto-k vector).

`bayes_influential_obs(result)` (in `helpers_bayes.R`) filters for k > threshold and returns a small table for display. Used by the dashboard's MCMC Diagnostics tab.

**If you have high-k observations:**
- Inspect them — are they obvious outliers / data-entry errors?
- If they're genuine but extreme, the model may be too rigid; consider a more flexible likelihood

### 7. Posterior P direction

**Threshold:** there isn't one, but values near 1.0 indicate a clear direction (most posterior mass is on one side of zero).

**What it measures:** P(effect > 0) under the posterior. A frequentist analog would be a one-sided p-value but with direct probability interpretation.

**Found in:** the `P_Direction` column on `result$pairwise_comparisons` for the Bayesian path.

**Use it for:** reporting effect-direction confidence without p-value semantics.

---

## Prior strength presets

`prior_strength` accepts five values. The four presets scale their priors to the data, so the same preset means the same thing in mm³ or cm³, in days or weeks. The priors actually used, with their values, are in `fit$summary$methods` (`prior_b`, `prior_intercept`, `prior_sd`, `prior_sigma`; for survival also `prior_aux`).

**Tumour growth** (`bayesian_tumor_growth()`), with y the modelled volume:

| Parameter | Prior | skeptical | informative | weakly_informative | diffuse |
|---|---|---|---|---|---|
| Intercept | normal(median y, 2.5 × MAD y) | | | | |
| Treatment main effects | normal(0, w) | w = 0.25 | 0.5 | 1 | 2.5 |
| Day slope and Treatment × Day | normal(0, m × range y / study span) | m = 1 | 1.5 | 2 | 5 |
| Animal and cage SDs, sigma | exponential(r / MAD y) | r = 2 | 2 | 1 | 0.5 |

w is a log-fold change; with `transform = "sqrt"` or `"none"` it is multiplied by MAD y (v0.28.0; before, a raw-mm³ treatment effect was pinned near 0, `CODE_REVIEW.md` R20.41).

**Survival** (`bayesian_survival()`), with t the event or censoring time:

| Parameter | Prior | skeptical | informative | weakly_informative | diffuse |
|---|---|---|---|---|---|
| Intercept (log time) | normal(median log t, 2.5 × MAD log t) | | | | |
| Treatment (log time ratio) | normal(0, w) | w = 0.25 | 0.5 | 1 | 2.5 |
| Cage frailty SD | exponential(r) | r = 2 | 2 | 1 | 0.5 |
| Weibull shape | lognormal(1, 1) under every preset | | | | |
| Log-normal sigma | exponential(1) under every preset | | | | |

The shape and sigma describe how event times are spread, not the treatment effect, so the ladder does not apply to them (v0.28.0). Under the ladder's exponential(2) the shape was shrunk toward 0.5, which pulls the Weibull hazard ratio, exp(−shape × b), toward 1: with a true shape of 6 and HR 0.088 the skeptical fit reported shape 3.86 and HR 0.246 (`CODE_REVIEW.md` R20.38).

`"manual"` uses `prior_b`, `prior_intercept`, `prior_sd` and `prior_sigma` (survival: `prior_aux`) as given.

### When to use which preset

| Situation | Pick |
|---|---|
| Default, you're not sure | `"skeptical"`. Requires strong data to support large estimated effects. Robust to small N. |
| You have prior data suggesting non-trivial effects | `"weakly_informative"` |
| You have explicit prior estimates from a prior trial | `"informative"` or `"manual"` (specify your prior centers) |
| You want to express genuine ignorance and have N ≥ 20 per group | `"diffuse"`. Don't use at small N — the result will look frequentist with wide CrIs |
| You're studying a known effect's magnitude (not whether it's non-zero) | `"manual"` with a prior centered on the expected effect |

### How the priors flow through

The `priors` argument to `bayesian_tumor_growth()` and `bayesian_survival()` takes a `tg_priors()` object. Inside the function, `.resolve_priors()` unpacks it into the individual prior arguments. This is purely a signature-cleanup mechanism.

The individual priors (`b`, `intercept`, `sd`, `sigma`) apply only with `strength = "manual"`, and then all of them must be given; the error names any that are missing. Under a preset they are ignored. (An earlier version of this guide suggested overriding one prior of a preset, which never did anything.)

```r
tg_priors(strength = "manual",
          b = "normal(0, 0.5)", intercept = "normal(6, 2)",
          sd = "exponential(1)", sigma = "exponential(1)")
```

---

## MCMC configuration

`tg_mcmc()` bundles the four common knobs.

```r
tg_mcmc(
  chains  = 4L,        # ≥ 4 for diagnostic reliability; 2 is enough for testing
  warmup  = 1000L,     # raise to 2000 if Rhat doesn't converge
  iter    = 500L,      # post-warmup. Raise to 1000-2000 if ESS is low
  seed    = 42L,       # set for reproducibility
  backend = "rstan"    # or "cmdstanr" for 3-10x faster Stan compilation
)
```

### Backend choice

| Backend | When to pick |
|---|---|
| `"rstan"` (default) | Works out of the box. |
| `"cmdstanr"` | Faster compilation. Needs the cmdstanr package and a CmdStan toolchain (`cmdstanr::install_cmdstan()`); without them the fit stops with installation instructions rather than falling back to rstan. |

**Compiled models are reused** (v0.28.0, `CODE_REVIEW.md` R20.78). Compiling the Stan model takes about 20 s, while sampling a study of this size takes a few seconds. The prior values are passed to Stan as data, so models with the same structure (formula, family, number of arms, random effects) share their Stan code, and a later fit in the same R session reuses the compiled model. Its draws are identical to a fresh fit's with the same seed. `fit$model_reused` says whether it happened. `clear_compiled_model_cache()` empties the cache, and `options(mouseExperiment.cache_compiled_models = FALSE)` turns it off.

### Chain count

- **4 chains** (default): the minimum for reliable Rhat / ESS computation. Use this unless you're sure why you're not.
- **2 chains**: faster; fine for development / testing. Don't report from a 2-chain fit.
- **8+ chains**: more parallelism if you have cores. Doesn't make the model fit better; just gives you more total samples in less wall-time.

### Iteration count

- **500 post-warmup × 4 chains = 2,000 draws**: the default. Adequate for most reportable runs.
- **1,000-2,000 post-warmup × 4 chains**: tighten CrIs when ESS is below 400 on any parameter you care about.
- **5,000+ post-warmup**: rarely needed; if you think you need this, the model is more likely ill-specified than under-sampled.

### Warmup count

- **1,000** (default): adequate for well-specified models. Stan auto-adapts during warmup so this is usually enough.
- **2,000**: bump if Rhat > 1.01 after the default. Often resolves "almost-converged" runs.
- Less than 1,000 isn't recommended.

---

## Prior–posterior overlays

`result$prior_posterior_plot` overlays the prior density (grey) and posterior density (blue) for each treatment-effect coefficient, with a vertical dashed line at zero. It answers the question "how much did the data update my prior?"

**Reading the overlay:**

- **Posterior centered far from prior center, narrow posterior:** strong evidence; data dominated the prior. ✓
- **Posterior close to prior center, narrow posterior:** prior and data agree, or the data was uninformative and the prior is what you're seeing. Check the data summary
- **Posterior similar to prior, wide posterior:** data was uninformative; the prior is doing all the work. Consider getting more data or pre-registering the prior choice
- **Posterior far from prior, but wide:** the data has pulled away from the prior but is itself uncertain. Honest uncertainty

The dashboard shows this plot in every Bayesian module's "Prior vs Posterior" tab.

---

## Worked examples

### Same data, three priors

```r
df <- combo_treatment_synthetic_data

skeptical <- bayesian_tumor_growth(
  df = df, priors = tg_priors(strength = "skeptical"),
  mcmc = tg_mcmc(seed = 1))

weakly <- bayesian_tumor_growth(
  df = df, priors = tg_priors(strength = "weakly_informative"),
  mcmc = tg_mcmc(seed = 1))

diffuse <- bayesian_tumor_growth(
  df = df, priors = tg_priors(strength = "diffuse"),
  mcmc = tg_mcmc(seed = 1))

# Compare treatment effects across the three:
skeptical$treatment_effects   # tight CrIs, smaller estimates (pulled toward 0)
weakly$treatment_effects      # moderate
diffuse$treatment_effects     # wide CrIs, close to maximum-likelihood estimates
```

What you should see: the skeptical prior shrinks estimates toward zero, narrowing CrIs but moving the central estimate. The diffuse prior approximates MLE, with wide CrIs. The weakly informative is in between. Pick the one that matches your prior beliefs honestly — don't pick the one that gives the answer you want.

### Diagnosing a non-converging fit

```r
fit <- bayesian_tumor_growth(df = small_df, priors = tg_priors(strength = "diffuse"))

# Check Rhat
print(fit$mcmc_diagnostics[fit$mcmc_diagnostics$Converged == FALSE, ])
# If anything appears here, the fit hasn't converged.

# Check NUTS
print(fit$nuts_diagnostics)
# If n_divergent > 0, the diffuse prior is likely the culprit at small N.

# Re-fit with a tighter prior + more warmup
fit2 <- bayesian_tumor_growth(
  df = small_df,
  priors = tg_priors(strength = "skeptical"),
  mcmc = tg_mcmc(warmup = 2000, iter = 1000))
```

### Reporting checklist (for a publication / report)

Always report:
- **Prior:** the priors in `fit$summary$methods`, with their values, e.g. "Treatment effects N(0, 0.25) on log volume; growth-rate terms N(0, 0.147); SDs Exponential(2.4)"
- **MCMC:** "4 chains × 1500 iterations (1000 warmup, 500 post-warmup) via `brms` 2.21 / Stan 2.32 with `cmdstanr` backend"
- **Convergence:** "All parameters had Rhat ≤ 1.01 and ESS ≥ 400. Zero divergent transitions."
- **Predictive coverage:** "Posterior predictive 95% intervals covered 94% of held-out data (n = X)."
- **LOO:** "PSIS-LOO elpd = X (SE Y); 2 observations had Pareto-k above 0.667 (inspected; not data-entry errors)."
- **Effect:** "Posterior median treatment effect = Z (95% CrI: A to B), P(effect > 0) = 0.98."

---

## Common pitfalls

| Pitfall | What to do |
|---|---|
| Reading effect estimates before checking diagnostics | Don't. If Rhat > 1.01 or there are divergences, the numbers are wrong. |
| Choosing a diffuse prior because it "feels more objective" | Diffuse priors aren't more objective — they encode a specific prior belief that effects can be arbitrarily large. At small N, this leads to wide CrIs that look like high uncertainty but are actually prior-driven |
| Reporting "P > 0.95" as a p-value substitute | It's not. It's a posterior probability of direction. Report it as `P(effect > 0) = 0.97`, not "p = 0.03" |
| Comparing models by AIC/BIC across the Bayesian and frequentist paths | Use LOO/WAIC for Bayesian model comparison. AIC/BIC apply to MLE estimates |
| Treating credible intervals as confidence intervals | They have different interpretations. A 95% CrI is "the true value is in this range with 95% posterior probability" — a 95% CI is "if I ran this study many times, 95% of intervals would contain the true value." The CrI is what you actually want to report |
| Re-running until you get the "right" answer | Don't. Set the seed, pick the prior up front, and report whatever the model says. If you don't like the result, run a sensitivity analysis across priors transparently |

---

## See also

- [`METHODS.md`](METHODS.md) — what each Bayesian fit actually fits (model formula, family choices)
- Roxygen `?bayesian_tumor_growth` and friends — argument-level reference
- `bayes_diagnostics_panel()` in `R/helpers_bayes.R` — the dashboard's diagnostic rendering
- `tg_priors()`, `tg_mcmc()` in `R/bayesian_config.R` — config helpers (CODE_REVIEW.md D.2)
- brms documentation: https://paul-buerkner.github.io/brms/
- Stan reference manual: https://mc-stan.org/users/documentation/
