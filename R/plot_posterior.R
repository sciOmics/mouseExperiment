#' Posterior distributions of a Bayesian fit's treatment effects
#'
#' Draws each treatment coefficient's posterior density, shading its
#' highest-density interval (HDI), with the prior's density optionally
#' overlaid (CODE_REVIEW.md R23.12). The dashed vertical line marks no effect
#' (0, on the model's scale: log volume for [bayesian_tumor_growth()], log
#' time for [bayesian_survival()]).
#'
#' * **The HDI** is the narrowest interval holding `hdi` of the posterior
#'   draws: every value inside it is more probable than every value outside.
#'   For a skewed posterior it differs from the equal-tailed credible interval
#'   the tables report. Each panel's title gives its limits.
#' * **The prior** is drawn over the posterior's range only, so a weak prior,
#'   much wider than the posterior, does not flatten the posterior to a spike;
#'   it then shows as a nearly flat line, which is what it is there.
#'
#' @param x A result of [bayesian_tumor_growth()] or [bayesian_survival()]
#'   (from 0.55.0, which keeps `prior_posterior_data`), or that data frame
#'   (`Parameter`, `Value`, `Source`).
#' @param show_prior Overlay each prior's density.
#' @param hdi The HDI's probability, between 0 and 1.
#' @param parameters The coefficients to draw (names in the `Parameter`
#'   column); all by default.
#' @return A ggplot.
#' @seealso [posterior_hdi()] for the limits as a table.
#' @export
#' @examples
#' \donttest{
#' if (requireNamespace("brms", quietly = TRUE)) {
#'   fit <- bayesian_tumor_growth(master_synthetic_data, reference_group = "Vehicle",
#'                                cage_column = "Cage", n_chains = 2, n_iter = 500,
#'                                plots = FALSE, verbose = FALSE)
#'   plot_posterior(fit, show_prior = TRUE, hdi = 0.9)
#' }
#' }
plot_posterior <- function(x, show_prior = FALSE, hdi = 0.95, parameters = NULL) {
  d <- me_posterior_data(x, parameters)
  lim <- posterior_hdi(d, hdi)
  label <- stats::setNames(sprintf("%s\n%s%% HDI [%s, %s]", lim$Parameter, format(100 * hdi),
                                   signif(lim$Lower, 3), signif(lim$Upper, 3)), lim$Parameter)
  curves <- do.call(rbind, lapply(unique(d$Parameter), function(p) {
    post <- d$Value[d$Parameter == p & d$Source == "Posterior"]
    span <- diff(range(post))
    from <- min(post) - 0.25 * span
    to <- max(post) + 0.25 * span
    dens <- stats::density(post, from = from, to = to, n = 256)
    out <- data.frame(Parameter = p, Source = "Posterior", x = dens$x, y = dens$y,
                      stringsAsFactors = FALSE)
    prior <- d$Value[d$Parameter == p & d$Source == "Prior"]
    if (isTRUE(show_prior) && length(prior) > 1L) {
      pd <- stats::density(prior, from = from, to = to, n = 256)
      out <- rbind(out, data.frame(Parameter = p, Source = "Prior", x = pd$x, y = pd$y,
                                   stringsAsFactors = FALSE))
    }
    out
  }))
  curves$Label <- factor(label[curves$Parameter], levels = label)
  post <- curves[curves$Source == "Posterior", , drop = FALSE]
  inside <- post[post$x >= lim$Lower[match(post$Parameter, lim$Parameter)] &
                   post$x <= lim$Upper[match(post$Parameter, lim$Parameter)], , drop = FALSE]
  p <- ggplot2::ggplot(mapping = ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_area(data = inside, fill = "steelblue", alpha = 0.35) +
    ggplot2::geom_line(data = post, colour = "steelblue4", linewidth = 0.7) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey30", linewidth = 0.5) +
    ggplot2::facet_wrap(~ Label, scales = "free") +
    ggplot2::labs(title = "Posterior distributions of the treatment effects",
                  subtitle = paste0("Shaded: the ", format(100 * hdi), "% highest-density interval",
                                    if (isTRUE(show_prior)) "; grey dashed: the prior" else "",
                                    "; dashed line: no effect (0)"),
                  x = "Coefficient value (model scale)", y = "Density") +
    ggplot2::theme_classic(base_size = 14)
  if (isTRUE(show_prior) && any(curves$Source == "Prior")) {
    p <- p + ggplot2::geom_line(data = curves[curves$Source == "Prior", , drop = FALSE],
                                colour = "grey45", linetype = "dashed", linewidth = 0.7)
  }
  p
}

#' Highest-density intervals of a Bayesian fit's treatment effects
#'
#' The narrowest interval holding `hdi` of each coefficient's posterior draws
#' (CODE_REVIEW.md R23.12).
#' @inheritParams plot_posterior
#' @return A data frame: `Parameter`, `Lower`, `Upper`, `Level`.
#' @export
#' @examples
#' d <- data.frame(Parameter = "x", Source = "Posterior", Value = stats::qnorm(ppoints(2000)))
#' posterior_hdi(d, 0.95)   # about -1.96 to 1.96
posterior_hdi <- function(x, hdi = 0.95, parameters = NULL) {
  if (!is.numeric(hdi) || length(hdi) != 1L || !(hdi > 0 && hdi < 1)) {
    stop("`hdi` must be one number between 0 and 1, e.g. 0.95.", call. = FALSE)
  }
  d <- me_posterior_data(x, parameters)
  out <- do.call(rbind, lapply(unique(d$Parameter), function(p) {
    v <- d$Value[d$Parameter == p & d$Source == "Posterior"]
    h <- me_hdi(v, hdi)
    data.frame(Parameter = p, Lower = h[1], Upper = h[2], Level = hdi, stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  out
}

#' The narrowest interval holding `prob` of the draws
#' @noRd
me_hdi <- function(x, prob) {
  x <- sort(x[is.finite(x)])
  n <- length(x)
  if (n < 2L) return(c(NA_real_, NA_real_))
  k <- max(1L, min(n - 1L, ceiling(prob * n) - 1L))
  width <- x[(k + 1L):n] - x[seq_len(n - k)]
  i <- which.min(width)
  c(x[i], x[i + k])
}

#' A result's prior and posterior draws, checked
#' @noRd
me_posterior_data <- function(x, parameters = NULL) {
  d <- if (is.data.frame(x)) x else x$prior_posterior_data
  if (!is.data.frame(d) || !nrow(d) || !all(c("Parameter", "Value", "Source") %in% names(d))) {
    stop("No posterior draws to draw: `x` must be a Bayesian result from mouseExperiment 0.55.0 ",
         "or later, or its prior_posterior_data.", call. = FALSE)
  }
  if (!is.null(parameters)) {
    unknown <- setdiff(parameters, d$Parameter)
    if (length(unknown)) stop("Unknown parameter(s): ", paste(unknown, collapse = ", "), call. = FALSE)
    d <- d[d$Parameter %in% parameters, , drop = FALSE]
  }
  d
}
