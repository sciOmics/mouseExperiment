# Copyright (c) 2026 mouseExperiment Contributors
# Licensed under the MIT License — see LICENSE file in the repo root.

#' Build a standard set of residual-diagnostic ggplots for an LMM / LM fit.
#'
#' Returns a list of three ggplot objects (Q-Q of residuals,
#' residuals vs fitted, scale-location). Used by every frequentist
#' fitting path so the dashboard can render them uniformly. Mirrors the
#' pattern established in `analyze_body_weight.R` (CODE_REVIEW.md J.12);
#' extracted to a helper in v0.4.8 so `tumor_growth_statistics` and
#' `tgs_path_auc` could adopt it.
#'
#' @param model An `lmerMod`, `lm`, or any object supporting
#'   `stats::residuals()` and `stats::fitted()`.
#' @param title_prefix Character prefix prepended to each plot title
#'   (e.g. "Body-weight LMM", "Tumour growth LMM").
#' @return Named list with `diag_qq_plot`, `diag_resid_fitted_plot`,
#'   `diag_scale_location_plot`. Any element is `NULL` if construction
#'   fails or `ggplot2` is unavailable.
#' @keywords internal
build_residual_diagnostic_plots <- function(model, title_prefix = "Model") {
  out <- list(diag_qq_plot             = NULL,
              diag_resid_fitted_plot   = NULL,
              diag_scale_location_plot = NULL)
  if (is.null(model)) return(out)
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(out)

  tryCatch({
    resid_vec <- stats::residuals(model)
    fit_vec   <- stats::fitted(model)
    qq_data   <- stats::qqnorm(resid_vec, plot.it = FALSE)
    qq_df     <- data.frame(theoretical = qq_data$x, sample = qq_data$y)

    # Use a compact theme so the plot fits a constrained plotOutput
    # without the title eating most of the canvas. The dashboard's h4
    # heading above each plot already names it; the ggplot title is
    # subordinate identification (useful when downloaded standalone).
    compact_theme <- ggplot2::theme_classic() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 11, face = "plain"),
        plot.subtitle = ggplot2::element_text(size = 9),
        plot.margin = ggplot2::margin(4, 6, 4, 6))

    out$diag_qq_plot <- ggplot2::ggplot(
      qq_df, ggplot2::aes(x = .data[["theoretical"]],
                          y = .data[["sample"]])
    ) +
      ggplot2::geom_point(alpha = 0.6) +
      ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                           colour = "red") +
      compact_theme +
      ggplot2::labs(title = paste0(title_prefix, ": Q-Q of residuals"),
                    x = "Theoretical Quantiles",
                    y = "Sample Quantiles")

    rf_df <- data.frame(Fitted = fit_vec, Residuals = resid_vec)
    out$diag_resid_fitted_plot <- ggplot2::ggplot(
      rf_df, ggplot2::aes(x = .data[["Fitted"]],
                          y = .data[["Residuals"]])
    ) +
      ggplot2::geom_point(alpha = 0.5) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                          colour = "red") +
      ggplot2::geom_smooth(method = "loess", se = FALSE,
                           colour = "steelblue", linewidth = 0.8,
                           formula = y ~ x) +
      compact_theme +
      ggplot2::labs(title = paste0(title_prefix, ": Residuals vs Fitted"),
                    x = "Fitted Values", y = "Residuals")

    # Scale-Location plot — sqrt(|residuals|) vs fitted. Heteroscedasticity
    # check: a clear trend in the loess line indicates non-constant variance.
    sl_df <- data.frame(Fitted = fit_vec,
                        SqrtAbsRes = sqrt(abs(resid_vec)))
    out$diag_scale_location_plot <- ggplot2::ggplot(
      sl_df, ggplot2::aes(x = .data[["Fitted"]],
                          y = .data[["SqrtAbsRes"]])
    ) +
      ggplot2::geom_point(alpha = 0.5) +
      ggplot2::geom_smooth(method = "loess", se = FALSE,
                           colour = "steelblue", linewidth = 0.8,
                           formula = y ~ x) +
      compact_theme +
      ggplot2::labs(title = paste0(title_prefix, ": Scale-Location"),
                    x = "Fitted Values",
                    y = expression(sqrt(abs("Residuals"))))
  }, error = function(e) NULL)

  out
}

#' Build a Q-Q plot for the random-effects BLUPs of an lme4 model.
#'
#' LMM assumes the random intercepts (and slopes) are normally
#' distributed. The fitted BLUPs are a proxy for this check. Returns
#' `NULL` when the model has no `id_column` random effect or when
#' `lme4` is unavailable.
#'
#' @param model A fitted `lmerMod` object.
#' @param id_column Character name of the grouping factor for the
#'   random intercept. Defaults to the first random-effect grouping
#'   in the model.
#' @param title_prefix Character prefix for the plot title.
#' @return A ggplot, or `NULL` on failure.
#' @keywords internal
build_random_effects_qq_plot <- function(model,
                                         id_column = NULL,
                                         title_prefix = "Model") {
  if (is.null(model)) return(NULL)
  if (!inherits(model, "lmerMod")) return(NULL)
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(NULL)
  if (!requireNamespace("lme4", quietly = TRUE))    return(NULL)

  tryCatch({
    re_list <- lme4::ranef(model)
    if (length(re_list) == 0L) return(NULL)
    grp <- if (!is.null(id_column) && id_column %in% names(re_list)) {
      id_column
    } else {
      names(re_list)[1L]
    }
    re_df <- as.data.frame(re_list[[grp]])
    # Stack random-effects columns (intercept + slope when present) into
    # long form so we can facet by component.
    if (ncol(re_df) == 0L) return(NULL)
    long <- do.call(rbind, lapply(colnames(re_df), function(cn) {
      vals <- re_df[[cn]]
      qq   <- stats::qqnorm(vals, plot.it = FALSE)
      data.frame(component  = cn,
                 theoretical = qq$x,
                 sample      = qq$y,
                 stringsAsFactors = FALSE)
    }))
    p <- ggplot2::ggplot(long,
                        ggplot2::aes(x = .data[["theoretical"]],
                                     y = .data[["sample"]])) +
      ggplot2::geom_point(alpha = 0.7) +
      ggplot2::geom_abline(slope = 1, intercept = 0,
                           linetype = "dashed", colour = "red") +
      ggplot2::theme_classic() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 11, face = "plain"),
        plot.subtitle = ggplot2::element_text(size = 9),
        plot.margin = ggplot2::margin(4, 6, 4, 6)) +
      ggplot2::labs(title = paste0(title_prefix, ": Q-Q of random effects"),
                    subtitle = paste0("Grouping: ", grp),
                    x = "Theoretical Quantiles",
                    y = "BLUP Sample Quantiles")
    if (length(unique(long$component)) > 1L) {
      p <- p + ggplot2::facet_wrap(~ component, scales = "free")
    }
    p
  }, error = function(e) NULL)
}

#' Compute LMM influence diagnostics (Cook's distance + DFBETAS).
#'
#' Uses `lme4::influence.merMod()` to compute case-deletion diagnostics.
#' Returns a list with two components: `cooks_distance` (data frame of
#' per-observation Cook's distance with a threshold flag) and
#' `dfbetas` (data frame of per-observation DFBETAS for each fixed
#' effect). Returns `NULL` for non-lmerMod input or when `lme4`'s
#' influence machinery fails.
#'
#' Cook's distance threshold: `4/n` is the conventional cut-off
#' (Fox 1991); a per-observation `cooks > 4/n` triggers the
#' `is_influential` flag.
#'
#' Caveat: `influence.merMod` refits the model leaving each group out
#' in turn. Cost is O(n_subjects); for large designs this can be slow.
#' We compute on the `obs = TRUE` (per-observation) basis for the most
#' useful Cook's distance display.
#'
#' @param model A fitted `lmerMod` object.
#' @return A list with `cooks_distance` and `dfbetas`, or `NULL`.
#' @keywords internal
build_lmm_influence <- function(model) {
  if (is.null(model)) return(NULL)
  if (!inherits(model, "lmerMod")) return(NULL)
  if (!requireNamespace("lme4", quietly = TRUE)) return(NULL)

  tryCatch({
    # Per-observation influence — most informative case-deletion view.
    infl <- lme4::influence(model, obs = TRUE)
    cd   <- stats::cooks.distance(infl)
    df   <- stats::dfbetas(infl)

    n <- length(cd)
    threshold <- 4 / max(n, 1L)

    cd_df <- data.frame(
      Obs           = seq_along(cd),
      Cooks_D       = round(as.numeric(cd), 4L),
      Threshold     = round(threshold, 4L),
      Is_Influential = as.numeric(cd) > threshold,
      stringsAsFactors = FALSE
    )

    # DFBETAS is one column per fixed effect.
    if (is.matrix(df)) {
      df_df <- as.data.frame(df)
      df_df$Obs <- seq_len(nrow(df_df))
      df_df <- df_df[, c("Obs", setdiff(names(df_df), "Obs")), drop = FALSE]
      # Round for display
      num_cols <- vapply(df_df, is.numeric, logical(1L))
      df_df[, num_cols] <- lapply(df_df[, num_cols, drop = FALSE], round, 4L)
    } else {
      df_df <- NULL
    }

    list(cooks_distance = cd_df, dfbetas = df_df)
  }, error = function(e) NULL)
}
