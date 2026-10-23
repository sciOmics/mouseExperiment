#' A priori power for a survival endpoint (log-rank)
#'
#' Animals per arm, or power, for comparing each treated arm's time to an event
#' (usually the tumour volume limit) with the control by the log-rank test, by
#' Schoenfeld's approximation (Schoenfeld 1983).
#'
#' * **Events needed** for one comparison with equal arms:
#'   \eqn{D = 4 (z_{1 - \alpha / 2} + z_{1 - \beta})^2 / (\log HR)^2}.
#' * **Animals per arm:** \eqn{D} divided by the expected number of events per
#'   pair of animals, \eqn{p_C + p_T}. The control's proportion with an event
#'   by the end of the study, \eqn{p_C}, is an input; under proportional
#'   hazards the treated arm's is \eqn{p_T = 1 - (1 - p_C)^{HR}}.
#' * **Power** for \eqn{n} per arm:
#'   \eqn{\Phi(\sqrt{n (p_C + p_T) / 4} \, |\log HR| - z_{1 - \alpha / 2})}.
#'
#' Each of the \eqn{k - 1} treated arms is compared with the control, at a
#' per-comparison alpha from the multiplicity correction.
#'
#' Schoenfeld's formula is asymptotic. With few events, as in studies of 8-12
#' animals per arm, the log-rank test's real power is usually somewhat lower,
#' so treat the N as a minimum.
#'
#' @param hazard_ratio The hazard ratio to detect, treated vs control (a single
#'   positive number other than 1; below 1 means treatment delays the event).
#' @param control_event_prob The proportion of control animals expected to have
#'   the event by the end of the study, in (0, 1].
#' @param n_groups Number of groups, the control included.
#' @param alpha Familywise two-sided alpha levels.
#' @param target_power Target power levels (`mode = "find_n"`).
#' @param n_per_group Animals per arm to evaluate (`mode = "find_power"`).
#' @param mode `"find_n"` or `"find_power"`.
#' @param n_comparisons Comparisons with the control; default `n_groups - 1`.
#' @param p_adjust_method `"bonferroni"` (default), `"dunnett"` (the
#'   comparisons' correlation is taken as 0.5, as with equal arms and a shared
#'   control) or `"none"`.
#' @param dropout_rate Expected proportion of animals lost without an event
#'   (censored early), in \[0, 1). `Enroll_N` allows for them.
#' @return A list:
#'   * `scenario_table`: for `"find_n"`, the events and animals per arm needed
#'     for each alpha and target power (`Required_Events`, `Required_N`,
#'     `Enroll_N`); for `"find_power"`, the expected events and power at each
#'     `n_per_group`;
#'   * `power_curve_data`: power against N per arm, 2 to 80;
#'   * `method_note`, `hazard_ratio`, `event_probs`.
#' @references Schoenfeld DA (1983). Sample-size formula for the
#'   proportional-hazards regression model. *Biometrics* 39, 499-503.
#' @seealso [apriori_power_analysis()] for a volume endpoint.
#' @export
#' @examples
#' # Animals per arm to detect a halving of the hazard, most controls reaching
#' # the volume limit, three arms with Dunnett's correction
#' apriori_survival_power(hazard_ratio = 0.3, control_event_prob = 0.9, n_groups = 3,
#'                        p_adjust_method = "dunnett")$scenario_table
apriori_survival_power <- function(hazard_ratio, control_event_prob, n_groups = 2L,
                                   alpha = c(0.05, 0.01), target_power = c(0.80, 0.90),
                                   n_per_group = NULL, mode = c("find_n", "find_power"),
                                   n_comparisons = NULL,
                                   p_adjust_method = c("bonferroni", "dunnett", "none"),
                                   dropout_rate = 0) {
  mode <- match.arg(mode)
  p_adjust_method <- match.arg(p_adjust_method)
  if (!is.numeric(hazard_ratio) || length(hazard_ratio) != 1L || !is.finite(hazard_ratio) ||
      hazard_ratio <= 0 || hazard_ratio == 1) {
    stop("'hazard_ratio' must be one positive number other than 1.", call. = FALSE)
  }
  if (!is.numeric(control_event_prob) || length(control_event_prob) != 1L ||
      !(control_event_prob > 0 && control_event_prob <= 1)) {
    stop("'control_event_prob' must be one number in (0, 1].", call. = FALSE)
  }
  if (!is.numeric(dropout_rate) || length(dropout_rate) != 1L ||
      dropout_rate < 0 || dropout_rate >= 1) {
    stop("'dropout_rate' must be a single number in [0, 1).", call. = FALSE)
  }
  n_groups <- as.integer(n_groups)
  if (n_groups < 2L) stop("'n_groups' must be at least 2.", call. = FALSE)
  m <- if (!is.null(n_comparisons)) as.integer(n_comparisons) else n_groups - 1L
  if (is.na(m) || m < 1L) stop("'n_comparisons' must be a positive integer.", call. = FALSE)

  p_c <- control_event_prob
  p_t <- 1 - (1 - p_c)^hazard_ratio
  per_pair <- p_c + p_t          # expected events from one animal in each arm
  log_hr <- abs(log(hazard_ratio))
  a_eff <- switch(p_adjust_method,
    bonferroni = alpha / m,
    dunnett    = vapply(alpha, me_dunnett_alpha, numeric(1), m = m),   # log-rank z: df = Inf
    none       = alpha)
  power_at <- function(n, a) stats::pnorm(sqrt(n * per_pair / 4) * log_hr - stats::qnorm(1 - a / 2))
  enroll <- function(n) ceiling(round(n / (1 - dropout_rate), 8))
  base <- function(i) data.frame(Hazard_Ratio = hazard_ratio, Control_Event_Prob = p_c,
                                 Treated_Event_Prob = p_t, N_Groups = n_groups,
                                 Alpha = alpha[i], Alpha_Per_Comparison = a_eff[i],
                                 N_Comparisons = m, stringsAsFactors = FALSE)

  scenario_table <- if (mode == "find_n") {
    do.call(rbind, unlist(lapply(seq_along(alpha), function(i) lapply(target_power, function(tp) {
      events <- 4 * (stats::qnorm(1 - a_eff[i] / 2) + stats::qnorm(tp))^2 / log_hr^2
      n <- max(2L, as.integer(ceiling(round(events / per_pair, 8))))
      cbind(base(i), Target_Power = tp, Required_Events = ceiling(events), Required_N = n,
            Enroll_N = enroll(n), Dropout_Rate = dropout_rate)
    })), recursive = FALSE))
  } else {
    if (!length(n_per_group)) stop("'n_per_group' must be supplied when mode = 'find_power'.",
                                   call. = FALSE)
    do.call(rbind, unlist(lapply(as.integer(n_per_group), function(n) lapply(seq_along(alpha), function(i) {
      n_a <- max(2L, floor(round(n * (1 - dropout_rate), 8)))
      cbind(base(i), N_Per_Group = n, N_Analysable = n_a, Expected_Events = n_a * per_pair,
            Dropout_Rate = dropout_rate, Achieved_Power = power_at(n_a, a_eff[i]))
    })), recursive = FALSE))
  }
  rownames(scenario_table) <- NULL

  n_seq <- 2:80
  power_curve_data <- do.call(rbind, lapply(seq_along(alpha), function(i) {
    n_a <- pmax(2L, floor(round(n_seq * (1 - dropout_rate), 8)))
    data.frame(Hazard_Ratio = hazard_ratio, N_Per_Group = n_seq, N_Analysable = n_a,
               Power = power_at(n_a, a_eff[i]), Alpha = alpha[i],
               Alpha_Per_Comparison = a_eff[i], stringsAsFactors = FALSE)
  }))

  method_note <- paste0(
    "Log-rank test, by Schoenfeld's approximation: each treated arm against the control, ",
    "at alpha ", paste(signif(a_eff, 3), collapse = " / "), " per comparison",
    switch(p_adjust_method,
      bonferroni = paste0(" (Bonferroni over ", m, " comparisons)"),
      dunnett    = paste0(" (Dunnett over ", m, " comparisons)"),
      none       = " (no multiplicity adjustment)"),
    ". Expected events by the end: ", round(100 * p_c), " % of controls and ",
    round(100 * p_t), " % of treated animals (proportional hazards). The approximation ",
    "is asymptotic; with few events the real power is usually somewhat lower, so treat ",
    "the N as a minimum.")

  list(scenario_table = scenario_table, power_curve_data = power_curve_data,
       method_note = method_note, hazard_ratio = hazard_ratio,
       event_probs = c(control = p_c, treated = p_t), mode = mode)
}
