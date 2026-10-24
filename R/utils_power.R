# Shared pieces of the analytic power calculations (CODE_REVIEW.md R22-F, F6).

# Dunnett's critical values, by (alpha, m, df). qmvt() integrates numerically,
# and the power curves ask for the same few values many times.
.dunnett_cache <- new.env(parent = emptyenv())

#' Dunnett's per-comparison alpha
#'
#' The two-sided alpha at which each of `m` comparisons with a shared control
#' is tested under Dunnett's correction: with equal group sizes the test
#' statistics have correlation 0.5, and the familywise critical value \eqn{c}
#' solves \eqn{P(\max_j |T_j| \le c) = 1 - \alpha}. The per-comparison alpha is
#' then \eqn{2 P(T > c)}. With `df = Inf` the statistics are normal, as for
#' the log-rank test. One comparison returns `alpha`.
#'
#' \eqn{c} is computed exactly (seeded, so reproducibly) at the degrees of
#' freedom of the classic tables and interpolated linearly in \eqn{1/df}
#' between them, as those tables are read: an exact value takes about 40 ms,
#' and a power curve needs one per point.
#' @param alpha Familywise two-sided alpha.
#' @param m Number of comparisons with the control.
#' @param df Error degrees of freedom; `Inf` for normal statistics.
#' @return The per-comparison alpha, between `alpha / m` (Bonferroni) and
#'   `alpha`.
#' @noRd
me_dunnett_alpha <- function(alpha, m, df = Inf) {
  m <- as.integer(m)
  if (m <= 1L) return(alpha)
  grid <- c(2:10, 12, 15, 20, 24, 30, 40, 60, 120, Inf)
  crit <- vapply(grid, function(g) me_dunnett_crit(alpha, m, g), numeric(1))
  q <- if (is.infinite(df)) crit[length(grid)] else
    stats::approx(1 / grid, crit, xout = 1 / max(df, 2))$y
  if (is.finite(df)) 2 * stats::pt(-q, df) else 2 * stats::pnorm(-q)
}

#' Dunnett's two-sided critical value, exactly, cached
#' @noRd
me_dunnett_crit <- function(alpha, m, df) {
  key <- paste(alpha, m, df)
  if (!is.null(hit <- .dunnett_cache[[key]])) return(hit)
  R <- matrix(0.5, m, m)
  diag(R) <- 1
  q <- me_with_seed(20260930L, if (is.finite(df)) {
    mvtnorm::qmvt(1 - alpha, tail = "both.tails", df = df, corr = R)$quantile
  } else {
    mvtnorm::qmvnorm(1 - alpha, tail = "both.tails", corr = R)$quantile
  })
  .dunnett_cache[[key]] <- q
  q
}
