# Copyright (c) 2026 mouseExperiment Contributors
# Licensed under the MIT License — see LICENSE file in the repo root.

#' Bayesian prior configuration object
#'
#' Bundles the five prior-related arguments accepted by every
#' `bayesian_*` entry point into a single named list so callers don't
#' have to pass `prior_strength = ..., prior_b = ..., prior_intercept =
#' ..., prior_sd = ..., prior_sigma = ...` five at a time. CODE_REVIEW.md
#' Round 2 D.2.
#'
#' The individual `prior_*` arguments on every Bayesian function continue
#' to work; when a `priors = tg_priors(...)` argument is supplied, its
#' fields take precedence over the legacy arguments. New callers should
#' prefer the helper:
#'
#' \preformatted{
#' fit <- bayesian_tumor_growth(df,
#'   priors = tg_priors(strength = "weakly_informative", sd = 1),
#'   mcmc   = tg_mcmc(chains = 4, iter = 800))
#' }
#'
#' @param strength One of "skeptical" (default), "weakly_informative",
#'   "informative", "diffuse", "manual". See each Bayesian function's
#'   `prior_strength` docs for what each mode does.
#' @param b,intercept,sd,sigma Optional manual `brmsprior` objects or
#'   string specifications, matching the `prior_b`, `prior_intercept`,
#'   `prior_sd`, `prior_sigma` arguments on `bayesian_tumor_growth()` et
#'   al. `NULL` (default) means "use the strength-derived default".
#' @return An object of class `tg_priors` (a tagged list).
#' @seealso [tg_mcmc()]
#' @examples
#' tg_priors()
#' tg_priors(strength = "weakly_informative", sd = 1)
#' @export
tg_priors <- function(strength  = c("skeptical", "weakly_informative",
                                    "informative", "diffuse", "manual"),
                      b         = NULL,
                      intercept = NULL,
                      sd        = NULL,
                      sigma     = NULL) {
  strength <- match.arg(strength)
  structure(
    list(strength = strength, b = b, intercept = intercept,
         sd = sd, sigma = sigma),
    class = c("tg_priors", "list")
  )
}

#' MCMC sampler configuration object
#'
#' Bundles the four-or-five MCMC-related arguments accepted by every
#' `bayesian_*` entry point into a single named list. CODE_REVIEW.md
#' Round 2 D.2.
#'
#' The individual `n_chains`, `n_warmup`, `n_iter`, `seed`, and `backend`
#' arguments on every Bayesian function continue to work; when an `mcmc =
#' tg_mcmc(...)` argument is supplied, its fields take precedence.
#'
#' @param chains Number of MCMC chains (default 4).
#' @param warmup Warm-up iterations per chain (default 1000).
#' @param iter Post-warm-up iterations per chain (default 500).
#' @param seed PRNG seed (default 42).
#' @param backend "rstan" (default) or "cmdstanr". Resolved by
#'   `resolve_brms_backend()` before being passed to `brms::brm()`, so
#'   if `cmdstanr` is requested but unavailable the call falls back to
#'   `rstan` rather than erroring.
#' @return An object of class `tg_mcmc` (a tagged list).
#' @seealso [tg_priors()]
#' @examples
#' tg_mcmc()
#' tg_mcmc(chains = 2, iter = 1000, backend = "cmdstanr")
#' @export
tg_mcmc <- function(chains  = 4L,
                    warmup  = 1000L,
                    iter    = 500L,
                    seed    = 42L,
                    backend = c("rstan", "cmdstanr")) {
  backend <- match.arg(backend)
  structure(
    list(chains  = as.integer(chains),
         warmup  = as.integer(warmup),
         iter    = as.integer(iter),
         seed    = as.integer(seed),
         backend = backend),
    class = c("tg_mcmc", "list")
  )
}

# ── Internal: resolve `priors`/`mcmc` config args against legacy individual
#    parameters. When `cfg` (a tg_priors / tg_mcmc list) is non-NULL its
#    fields take precedence; when `cfg` is NULL, the legacy parameters are
#    used. Returns a list keyed by the field names of the helper. -----------
.resolve_priors <- function(priors, prior_strength, prior_b, prior_intercept,
                            prior_sd, prior_sigma) {
  if (!is.null(priors)) {
    if (!inherits(priors, "tg_priors")) {
      stop("`priors` must be a `tg_priors()` object or NULL.", call. = FALSE)
    }
    return(list(strength  = priors$strength,
                b         = priors$b,
                intercept = priors$intercept,
                sd        = priors$sd,
                sigma     = priors$sigma))
  }
  list(strength  = prior_strength,
       b         = prior_b,
       intercept = prior_intercept,
       sd        = prior_sd,
       sigma     = prior_sigma)
}

.resolve_mcmc <- function(mcmc, n_chains, n_warmup, n_iter, seed, backend) {
  if (!is.null(mcmc)) {
    if (!inherits(mcmc, "tg_mcmc")) {
      stop("`mcmc` must be a `tg_mcmc()` object or NULL.", call. = FALSE)
    }
    return(list(chains  = mcmc$chains,
                warmup  = mcmc$warmup,
                iter    = mcmc$iter,
                seed    = mcmc$seed,
                backend = mcmc$backend))
  }
  list(chains  = n_chains,
       warmup  = n_warmup,
       iter    = n_iter,
       seed    = seed,
       backend = backend)
}

#' @export
print.tg_priors <- function(x, ...) {
  cat("<tg_priors>\n")
  cat("  strength :", x$strength, "\n")
  for (k in c("b", "intercept", "sd", "sigma")) {
    v <- x[[k]]
    if (!is.null(v)) cat(sprintf("  %s: %s\n", format(k, width = 9), format(v)))
  }
  invisible(x)
}

#' @export
print.tg_mcmc <- function(x, ...) {
  cat("<tg_mcmc>\n")
  cat(sprintf("  chains : %d, warmup : %d, iter : %d, seed : %d, backend : %s\n",
              x$chains, x$warmup, x$iter, x$seed, x$backend))
  invisible(x)
}
