#' Accelerated failure time (AFT) survival models
#'
#' Compares each arm's time to the endpoint with the reference arm's by a
#' parametric accelerated failure time model ([survival::survreg()]), and
#' reports **time ratios**: how many times longer (or shorter) an arm's
#' animals take to reach the endpoint than the reference's. A time ratio of
#' 1.8 means the treatment stretches the time to the endpoint by 80 %.
#' (CODE_REVIEW.md R23.11.)
#'
#' Unlike the Cox model, an AFT model does not assume proportional hazards; it
#' assumes a shape for the distribution of times instead.
#' * **The distributions:** Weibull, log-normal and log-logistic. With
#'   `distribution = "auto"` all three are fitted and the one with the lowest
#'   AIC is used; `distribution_comparison` reports each one's AIC.
#' * **The time ratio** is `exp(coefficient)`, with a Wald interval and
#'   p-value on the log scale, the p-values corrected over the arms compared
#'   with the reference (`p_adjust_method`).
#' * **An arm with no events** has no finite time ratio (its times are all
#'   censored): it is left out of the fit and reported with a note.
#' * **Each arm's median,** from the fitted model, is in `Median_Predicted`;
#'   `fitted_curves` holds each arm's fitted survival curve, to compare with
#'   the Kaplan-Meier curves.
#' * **Animals are treated as independent.** When whole cages were
#'   randomised (`randomisation_unit = "cage"`) the p-values overstate the
#'   evidence; the function warns.
#'
#' @param df One row per animal: its time and whether the endpoint was
#'   reached.
#' @param time_column,event_column,treatment_column,id_column Column names.
#'   The event is 1 (reached) or 0 (censored).
#' @param cage_column Cage column, or NULL; used only to key animals.
#' @param dose_column Dose column, or NULL. When a treatment is given at more
#'   than one dose, each dose is its own arm, as in [survival_statistics()].
#' @param reference_group The reference arm. NULL picks the control by name,
#'   as the other analyses do.
#' @param distribution `"auto"` (lowest AIC), `"weibull"`, `"lognormal"` or
#'   `"loglogistic"`.
#' @param p_adjust_method Correction of the p-values over the arms compared
#'   with the reference: `"holm"`, `"bonferroni"`, `"fdr"` or `"none"`.
#' @param randomisation_unit `"mouse"` (default) or `"cage"`; with `"cage"`
#'   the function warns that its p-values treat animals as independent.
#' @param conf_level Level of the intervals.
#' @param return_model Keep the `survreg` fit in `$model`.
#' @return A list of class `me_aft`:
#'   * `time_ratios`: per arm, `N`, `Events`, `Time_Ratio`, `Lower_CL`,
#'     `Upper_CL`, `P_Value`, `P_Adjusted`, `Median_Predicted` and `Note`;
#'   * `distribution`, the one used, and `distribution_comparison` (AIC and
#'     log-likelihood of each distribution fitted);
#'   * `scale` (the survreg scale; for the Weibull, the shape is `1 / scale`);
#'   * `fitted_curves`: `Treatment`, `Time`, `Survival`;
#'   * `survival_data`: the data as fitted (`Time`, `Event`, `Treatment`);
#'   * `reference_group`, `p_adjust_method`, `conf_level`, `model`, `meta`.
#' @seealso [survival_statistics()] for Cox and log-rank, [bayesian_survival()]
#'   for the Bayesian Weibull or log-normal AFT model.
#' @export
#' @examples
#' ev <- survival_events(master_synthetic_data, volume_column = "Volume",
#'                       volume_limit = 1000, cage_column = "Cage")
#' aft <- aft_survival(ev, reference_group = "Vehicle", cage_column = "Cage")
#' aft$time_ratios
aft_survival <- function(df, time_column = "Time", event_column = "Event",
                         treatment_column = "Treatment", id_column = "ID",
                         cage_column = NULL, dose_column = NULL, reference_group = NULL,
                         distribution = c("auto", "weibull", "lognormal", "loglogistic"),
                         p_adjust_method = c("holm", "bonferroni", "fdr", "none"),
                         randomisation_unit = c("mouse", "cage"),
                         conf_level = 0.95, return_model = FALSE) {
  distribution <- match.arg(distribution)
  p_adjust_method <- match.arg(p_adjust_method)
  randomisation_unit <- match.arg(randomisation_unit)
  for (col in c(time_column, event_column, treatment_column)) {
    if (!has_col(df, col)) stop("Column '", col, "' is not in the data.", call. = FALSE)
  }
  if (has_col(df, dose_column)) {
    df[[treatment_column]] <- me_arms_with_dose(df[[treatment_column]], df, dose_column, "aft_survival")
  }
  d <- data.frame(Time = suppressWarnings(as.numeric(df[[time_column]])),
                  Event = suppressWarnings(as.integer(as.character(df[[event_column]]))),
                  Treatment = as.character(df[[treatment_column]]),
                  stringsAsFactors = FALSE)
  if (!all(d$Event %in% c(0L, 1L, NA))) {
    stop("'", event_column, "' must contain only 0 (censored) and 1 (event).", call. = FALSE)
  }
  validate_one_row_per_subject(df, id_column, treatment_column,
                               if (has_col(df, cage_column)) cage_column,
                               caller = "aft_survival()")
  d <- d[!is.na(d$Treatment) & is.finite(d$Time) & !is.na(d$Event), , drop = FALSE]
  if (any(d$Time <= 0)) {
    stop("AFT models need positive times; ", sum(d$Time <= 0), " animal(s) have a time of 0 or less.",
         call. = FALSE)
  }
  reference_group <- me_resolve_reference(d$Treatment, reference_group, "aft_survival")
  if (!any(d$Event[d$Treatment == reference_group] == 1L)) {
    stop("The reference arm ('", reference_group, "') has no events, so no time ratio can be estimated.",
         call. = FALSE)
  }
  if (randomisation_unit == "cage") {
    warning("Whole cages were randomised, but the AFT model treats animals as independent: ",
            "its p-values and intervals overstate the evidence.", call. = FALSE)
  }

  arms <- c(reference_group, setdiff(sort(unique(d$Treatment)), reference_group))
  n_by <- table(factor(d$Treatment, levels = arms))
  ev_by <- tapply(d$Event, factor(d$Treatment, levels = arms), sum)
  no_events <- arms[ev_by == 0L]
  fit_d <- d[!d$Treatment %in% no_events, , drop = FALSE]
  fit_d$Treatment <- factor(fit_d$Treatment, levels = setdiff(arms, no_events))

  dists <- if (distribution == "auto") c("weibull", "lognormal", "loglogistic") else distribution
  fits <- lapply(dists, function(dd) tryCatch(
    survival::survreg(survival::Surv(Time, Event) ~ Treatment, data = fit_d, dist = dd),
    error = function(e) e))
  names(fits) <- dists
  ok <- !vapply(fits, inherits, logical(1), "error")
  if (!any(ok)) stop("The AFT model could not be fitted: ", conditionMessage(fits[[1]]), call. = FALSE)
  comparison <- data.frame(
    Distribution = dists,
    LogLik = vapply(seq_along(fits), function(i) if (ok[i]) as.numeric(stats::logLik(fits[[i]])) else NA_real_, 0),
    AIC = vapply(seq_along(fits), function(i) if (ok[i]) stats::AIC(fits[[i]]) else NA_real_, 0),
    stringsAsFactors = FALSE)
  chosen <- comparison$Distribution[which.min(comparison$AIC)]
  comparison$Used <- comparison$Distribution == chosen
  fit <- fits[[chosen]]

  # Time ratios: exp(coefficient), Wald on the log scale.
  co <- summary(fit)$table
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  terms <- paste0("Treatment", setdiff(levels(fit_d$Treatment), reference_group))
  est <- co[terms, "Value"]
  se <- co[terms, "Std. Error"]
  p <- co[terms, "p"]
  p_adj <- if (length(p)) stats::p.adjust(p, method = p_adjust_method) else numeric(0)
  med <- stats::predict(fit, newdata = data.frame(Treatment = factor(levels(fit_d$Treatment),
                                                                     levels = levels(fit_d$Treatment))),
                        type = "quantile", p = 0.5)
  names(med) <- levels(fit_d$Treatment)
  fitted_arm <- setdiff(levels(fit_d$Treatment), reference_group)
  time_ratios <- data.frame(
    Treatment = arms,
    N = as.integer(n_by[arms]),
    Events = as.integer(ev_by[arms]),
    Time_Ratio = NA_real_, Lower_CL = NA_real_, Upper_CL = NA_real_,
    P_Value = NA_real_, P_Adjusted = NA_real_,
    Median_Predicted = unname(med[arms]),
    Note = "", stringsAsFactors = FALSE)
  i <- match(fitted_arm, arms)
  time_ratios$Time_Ratio[i] <- exp(est)
  time_ratios$Lower_CL[i] <- exp(est - z * se)
  time_ratios$Upper_CL[i] <- exp(est + z * se)
  time_ratios$P_Value[i] <- p
  time_ratios$P_Adjusted[i] <- p_adj
  time_ratios$Time_Ratio[1] <- 1
  time_ratios$Note[1] <- "reference"
  time_ratios$Note[time_ratios$Treatment %in% no_events] <-
    "no events: every time is censored, so no time ratio can be estimated"
  rownames(time_ratios) <- NULL

  # Each arm's fitted survival curve, over the observed times.
  grid <- seq(0, max(d$Time), length.out = 101)[-1]
  lp <- stats::predict(fit, newdata = data.frame(Treatment = factor(levels(fit_d$Treatment),
                                                                    levels = levels(fit_d$Treatment))),
                       type = "lp")
  fitted_curves <- do.call(rbind, lapply(seq_along(lp), function(k) data.frame(
    Treatment = levels(fit_d$Treatment)[k], Time = grid,
    Survival = 1 - survival::psurvreg(grid, mean = lp[k], scale = fit$scale, distribution = chosen),
    stringsAsFactors = FALSE)))
  rownames(fitted_curves) <- NULL

  out <- list(
    time_ratios = time_ratios,
    distribution = chosen,
    distribution_comparison = comparison,
    scale = fit$scale,
    fitted_curves = fitted_curves,
    survival_data = d,
    reference_group = reference_group,
    p_adjust_method = p_adjust_method,
    conf_level = conf_level,
    randomisation_unit = randomisation_unit,
    model = if (isTRUE(return_model)) fit,
    meta = me_result_meta(
      analysis_type = paste0("Accelerated failure time model (", chosen, ")"),
      model_type_used = "aft", inference = "frequentist", interval_type = "confidence",
      estimate_scale = "time ratio", comparison_family = "vs reference",
      p_adjust_method = p_adjust_method))
  me_result(out, "me_aft")
}

#' Kaplan-Meier curves with an AFT model's fitted curves
#'
#' The check on an AFT model's distribution (R23.11): each arm's Kaplan-Meier
#' curve (steps) with the model's fitted curve (dashed). Where the dashed
#' curves follow the steps, the distribution fits; where they part, compare
#' the other distributions' AIC.
#' @param res An [aft_survival()] result.
#' @param colors Named colours by arm, or NULL.
#' @return A ggplot.
#' @keywords internal
aft_fit_plot <- function(res, colors = NULL) {
  d <- res$survival_data
  km <- survival::survfit(survival::Surv(Time, Event) ~ Treatment, data = d)
  arm <- sub("^Treatment=", "", rep(names(km$strata), km$strata))
  steps <- do.call(rbind, lapply(unique(arm), function(a) {
    i <- arm == a
    data.frame(Treatment = a, Time = c(0, km$time[i]), Survival = c(1, km$surv[i]),
               stringsAsFactors = FALSE)
  }))
  fitted <- res$fitted_curves
  p <- ggplot2::ggplot(mapping = ggplot2::aes(x = Time, y = Survival, colour = Treatment)) +
    ggplot2::geom_step(data = steps, linewidth = 0.8) +
    ggplot2::geom_line(data = fitted, linetype = "dashed", linewidth = 0.8) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(title = sprintf("Kaplan-Meier (steps) and the %s AFT fit (dashed)", res$distribution),
                  x = "Time (days)", y = "Proportion without the event", colour = NULL) +
    ggplot2::theme_classic(base_size = 14)
  if (length(colors)) p <- p + ggplot2::scale_colour_manual(values = colors)
  p
}
