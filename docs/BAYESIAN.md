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
print(fit$loo_diagnostics)         # elpd_loo + n_high_k (Pareto-k > 0.7)
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

**Found in:** `result$ppc_coverage` (single-row data frame: `cov_50`, `cov_80`, `cov_95`, `n_obs`).

**If it's off:** the model's predictive distribution is too narrow (under-coverage) or too wide (over-coverage). Check the residuals; consider a richer likelihood (e.g., Student-t instead of Normal if outliers are an issue).

### 6. PSIS-LOO + Pareto-k

**Threshold:** all Pareto-k < 0.7. A handful of observations with 0.5 < k < 0.7 are acceptable but worth listing.

**What it measures:** leave-one-out cross-validation via Pareto-smoothed importance sampling. Pareto-k flags observations the LOO approximation can't reliably estimate — usually highly influential observations.

**Found in:** `result$loo_diagnostics` (`elpd_loo`, `se_elpd`, `n_high_k`, plus the list-column `pareto_k` with the per-observation Pareto-k vector).

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

`prior_strength` accepts five values. Each one sets a specific brms prior across fixed effects (`b`), intercept, and variance components (`sd`, `sigma`). The presets are calibrated for tumor-growth-scale data (log mm³), so they make sense as defaults; for other scales, set `prior_strength = "manual"` and supply your own.

```r
# Skeptical (default) — strong shrinkage toward zero effects
b ~ N(0, 0.25)
sd, sigma ~ Exponential(2)

# Weakly informative — moderate shrinkage
b ~ N(0, 1)
sd, sigma ~ Exponential(1)

# Informative — light shrinkage; assumes effects are moderate
b ~ N(0, 0.5)
sd, sigma ~ Exponential(2)

# Diffuse — minimal shrinkage; close to a flat prior on the log scale
b ~ N(0, 2.5)
sd, sigma ~ Exponential(0.5)

# Manual — supply prior_b, prior_intercept, prior_sd, prior_sigma directly
```

### When to use which preset

| Situation | Pick |
|---|---|
| Default, you're not sure | `"skeptical"`. Requires strong data to support large estimated effects. Robust to small N. |
| You have prior data suggesting non-trivial effects | `"weakly_informative"` |
| You have explicit prior estimates from a prior trial | `"informative"` or `"manual"` (specify your prior centers) |
| You want to express genuine ignorance and have N ≥ 20 per group | `"diffuse"`. Don't use at small N — the result will look frequentist with wide CrIs |
| You're studying a known effect's magnitude (not whether it's non-zero) | `"manual"` with a prior centered on the expected effect |

### How the priors flow through

The `priors` argument to `bayesian_tumor_growth()` (and friends) takes a `tg_priors()` object. Inside the function, `.resolve_priors()` unpacks it into the individual prior arguments. This is purely a signature-cleanup mechanism — the priors that hit brms are the same as if you'd passed `prior_strength = "weakly_informative", prior_b = NULL, ...`.

To override only some priors:

```r
# Use the "weakly_informative" preset but tighten the intercept prior
tg_priors(strength = "weakly_informative",
          intercept = "normal(0, 0.5)")
```

Any field you leave NULL takes its value from the strength preset.

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

| Backend | Compile time | Run time | When to pick |
|---|---|---|---|
| `"rstan"` (default) | Slow (60-120 s first time, cached after) | Same | Default. Works out of the box. |
| `"cmdstanr"` | Fast (10-30 s first time) | Same | Production. VPS. Anywhere you'll rerun the same formula often. Install via `cmdstanr::install_cmdstan()` |

`resolve_brms_backend()` in `R/utils_bayes.R` falls back to rstan if cmdstanr is requested but unavailable — you don't need to feature-detect at the call site.

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
- **Prior:** "We used a weakly informative prior (b ~ N(0, 1); sd, sigma ~ Exponential(1))"
- **MCMC:** "4 chains × 1500 iterations (1000 warmup, 500 post-warmup) via `brms` 2.21 / Stan 2.32 with `cmdstanr` backend"
- **Convergence:** "All parameters had Rhat ≤ 1.01 and ESS ≥ 400. Zero divergent transitions."
- **Predictive coverage:** "Posterior predictive 95% intervals covered 94% of held-out data (n = X)."
- **LOO:** "PSIS-LOO elpd = X (SE Y); 2 observations had Pareto-k > 0.7 (inspected; not data-entry errors)."
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
