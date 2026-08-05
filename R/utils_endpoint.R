# Copyright (c) 2026 mouseExperiment Contributors
# Licensed under the MIT License - see LICENSE file

# Endpoint-day estimands.
#
# CODE_REVIEW.md R3.5 / R3.3 / G.3 — five functions computed their headline
# metric as a raw mean at the global last study day:
#
#   max_day <- max(wd$Day); final <- wd[wd$Day == max_day, ]
#   ctrl_mean <- mean(final$Volume[final$Treatment == reference_group])
#
# Animals leave these studies *because their tumours got large*, so conditioning
# on survival to max_day selects the slowest-growing animals in every arm — and
# most severely in the control arm, which loses animals earliest. Every TGI
# formed against that denominator is biased downward (efficacy understated), and
# when no control animal reaches max_day the denominator is NaN and every TGI
# silently becomes NaN.
#
# The maintainer's requirement (G.3) is to incorporate the dropout, not discard
# or truncate around it. The mechanism makes that achievable: euthanasia
# triggered by a *measured* volume crossing a protocol threshold is missingness
# that depends only on observed data, i.e. Missing At Random — under which
# likelihood-based methods are valid without modelling the dropout process.
#
# So the default here fits a log-scale mixed model to ALL observations from ALL
# animals and reads off each arm's marginal mean at day T. Nothing is discarded,
# nothing is truncated, and no synthetic data rows are fabricated: the
# extrapolation is the model's, and its uncertainty is the model's uncertainty.
# Back-transforming the log-scale EMM with exp() yields a geometric mean, which
# is also the right centre for a log-normal quantity and resolves the
# arithmetic/geometric inconsistency noted in R3.6.

#' Endpoint estimand methods
#' @noRd
#' @keywords internal
ME_ENDPOINT_METHODS <- c("model", "last_obs", "survivors")

#' Per-arm endpoint volumes under a chosen estimand
#'
#' @param df Long data frame.
#' @param id_column,treatment_column,time_column,volume_column Column names.
#' @param cage_column Optional; used only for the composite mouse key.
#' @param endpoint_day Day at which to evaluate. `NULL` uses the maximum
#'   observed day.
#' @param endpoint_method One of:
#'   \describe{
#'     \item{`"model"`}{(default) group geometric means at `endpoint_day` from a
#'       log-scale LMM fitted to every observation. Uses all animals.}
#'     \item{`"last_obs"`}{each animal's own last observation, as
#'       `bayesian_dose_response()` does since v0.4.14. Uses all animals but
#'       evaluates them at *different days*, so it is not an estimate of volume
#'       at `endpoint_day`: an animal removed on day 12 contributes its day-12
#'       volume. On a simulated study with volume-triggered euthanasia this was
#'       *more* biased than `"survivors"`, not less, because it understates the
#'       control arm most. Retained as a fallback for when the model cannot be
#'       fitted, not as an equal alternative.}
#'     \item{`"survivors"`}{raw mean among animals observed at `endpoint_day` —
#'       the pre-0.8.0 behaviour, retained for reproducibility. Warns.}
#'   }
#' @return A list with `group_means` (data frame: Treatment, Mean_Volume, N,
#'   plus SE / CI bounds under `"model"`), `per_mouse` (Treatment, MouseKey,
#'   Volume, Day_Used) for the resampling-based methods, `n_at_risk` (animals
#'   observed at `endpoint_day`, per arm), `endpoint_day`, `method`, and
#'   `attrition` (a one-row summary of how many animals were lost).
#' @noRd
#' @keywords internal
endpoint_volumes <- function(df,
                             id_column = "ID",
                             treatment_column = "Treatment",
                             time_column = "Day",
                             volume_column = "Volume",
                             cage_column = NULL,
                             endpoint_day = NULL,
                             endpoint_method = c("model", "last_obs", "survivors")) {
  endpoint_method <- match.arg(endpoint_method)

  d <- data.frame(
    MouseKey  = if (!is.null(cage_column) && cage_column %in% names(df)) {
      make_mouse_key(as.character(df[[treatment_column]]),
                     as.character(df[[id_column]]),
                     as.character(df[[cage_column]]))
    } else {
      make_mouse_key(as.character(df[[treatment_column]]),
                     as.character(df[[id_column]]))
    },
    Treatment = as.character(df[[treatment_column]]),
    Day       = as.numeric(df[[time_column]]),
    Volume    = as.numeric(df[[volume_column]]),
    stringsAsFactors = FALSE
  )
  d <- d[is.finite(d$Day) & is.finite(d$Volume), , drop = FALSE]
  if (nrow(d) == 0L) stop("No usable volume observations.", call. = FALSE)

  ep_day <- if (is.null(endpoint_day)) max(d$Day) else as.numeric(endpoint_day)

  # Attrition bookkeeping — the numbers that make survivor selection visible.
  all_mice   <- unique(d[, c("MouseKey", "Treatment")])
  at_risk    <- unique(d[d$Day >= ep_day, c("MouseKey", "Treatment")])
  n_total    <- table(all_mice$Treatment)
  n_at_risk  <- table(factor(at_risk$Treatment, levels = names(n_total)))
  attrition  <- data.frame(
    Treatment    = names(n_total),
    N_Enrolled   = as.integer(n_total),
    N_At_Endpoint = as.integer(n_at_risk),
    Pct_Lost     = round(100 * (1 - as.integer(n_at_risk) / as.integer(n_total)), 1),
    stringsAsFactors = FALSE
  )

  if (any(attrition$N_At_Endpoint < attrition$N_Enrolled)) {
    lost <- attrition[attrition$N_At_Endpoint < attrition$N_Enrolled, ]
    msg <- paste(sprintf("%s: %d/%d", lost$Treatment,
                         lost$N_At_Endpoint, lost$N_Enrolled), collapse = "; ")
    if (endpoint_method == "survivors") {
      warning("endpoint_method = 'survivors' conditions on being observed at ",
              "day ", ep_day, ", but animals were lost before then (", msg,
              "). Animals leave because their tumours grew, so this selects the ",
              "slowest growers — most severely in the control arm — and biases ",
              "TGI downward. Use endpoint_method = 'model' (default) or ",
              "'last_obs'.", call. = FALSE)
    } else {
      message("Animals lost before day ", ep_day, " (", msg,
              "); all of them still contribute under endpoint_method = '",
              endpoint_method, "'.")
    }
  }

  per_mouse <- NULL
  group_means <- NULL

  if (endpoint_method == "model") {
    gm <- model_endpoint_means(d, ep_day)
    if (is.null(gm)) {
      warning("Model-based endpoint means could not be fitted; falling back to ",
              "endpoint_method = 'last_obs'.", call. = FALSE)
      endpoint_method <- "last_obs"
    } else {
      group_means <- gm
    }
  }

  if (endpoint_method %in% c("last_obs", "survivors")) {
    per_mouse <- if (endpoint_method == "last_obs") {
      do.call(rbind, lapply(split(d, d$MouseKey, drop = TRUE), function(s) {
        s <- s[order(s$Day), ]
        # Each animal's last observation at or before the endpoint day.
        s <- s[s$Day <= ep_day, , drop = FALSE]
        if (nrow(s) == 0L) return(NULL)
        data.frame(MouseKey = s$MouseKey[1], Treatment = s$Treatment[1],
                   Volume = s$Volume[nrow(s)], Day_Used = s$Day[nrow(s)],
                   stringsAsFactors = FALSE)
      }))
    } else {
      sv <- d[d$Day == ep_day, , drop = FALSE]
      data.frame(MouseKey = sv$MouseKey, Treatment = sv$Treatment,
                 Volume = sv$Volume, Day_Used = sv$Day,
                 stringsAsFactors = FALSE)
    }
    group_means <- do.call(rbind, lapply(
      split(per_mouse, per_mouse$Treatment, drop = TRUE),
      function(g) data.frame(
        Treatment   = g$Treatment[1],
        Mean_Volume = mean(g$Volume, na.rm = TRUE),
        N           = sum(is.finite(g$Volume)),
        stringsAsFactors = FALSE)
    ))
    rownames(group_means) <- NULL
  }

  list(group_means = group_means, per_mouse = per_mouse,
       attrition = attrition, endpoint_day = ep_day,
       method = endpoint_method)
}

#' Group geometric mean volumes at a given day from a log-scale LMM
#'
#' Fits `log(Volume) ~ Treatment * Day + (1 | MouseKey)` to every observation
#' and returns each arm's marginal mean at `ep_day`, back-transformed. Because
#' the model uses all data from all animals, an animal euthanised on day 20 still
#' informs its arm's intercept and slope at day 35 — which is what makes this
#' valid under the MAR dropout mechanism these studies have, and why it does not
#' need the synthetic-row imputation removed in R3.9.
#'
#' @param d Data frame with `MouseKey`, `Treatment`, `Day`, `Volume`.
#' @param ep_day Day at which to marginalise.
#' @return Data frame with Treatment, Mean_Volume (geometric), SE_log,
#'   Lower_CL, Upper_CL, N; or `NULL` if the model or emmeans step fails.
#' @noRd
#' @keywords internal
model_endpoint_means <- function(d, ep_day) {
  pos <- d$Volume[is.finite(d$Volume) & d$Volume > 0]
  if (length(pos) == 0L) return(NULL)

  d$.logv <- log(pmax(d$Volume, min(pos) / 2))
  d$Treatment <- factor(d$Treatment)

  # A random intercept per animal needs at least two animals per arm and more
  # than one timepoint; fall back to NULL (caller uses last_obs) otherwise.
  if (length(unique(d$Day)) < 2L || nrow(d) < 6L) return(NULL)

  fit <- tryCatch(
    withCallingHandlers(
      lme4::lmer(.logv ~ Treatment * Day + (1 | MouseKey), data = d,
                 control = lme4::lmerControl(check.nobs.vs.nlev = "ignore",
                                             check.nobs.vs.nRE  = "ignore")),
      warning = function(w) invokeRestart("muffleWarning")
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)

  emm <- tryCatch(
    summary(emmeans::emmeans(fit, specs = "Treatment",
                             at = list(Day = ep_day))),
    error = function(e) NULL
  )
  if (is.null(emm)) return(NULL)

  n_by <- table(unique(d[, c("MouseKey", "Treatment")])$Treatment)

  data.frame(
    Treatment   = as.character(emm$Treatment),
    # exp() of a log-scale marginal mean is a geometric mean — the right centre
    # for a log-normal quantity, and consistent with the modelling scale.
    Mean_Volume = exp(emm$emmean),
    SE_log      = emm$SE,
    Lower_CL    = exp(emm$lower.CL),
    Upper_CL    = exp(emm$upper.CL),
    N           = as.integer(n_by[as.character(emm$Treatment)]),
    stringsAsFactors = FALSE
  )
}

#' TGI from a group-means table, with the reference arm pinned at 0
#'
#' @param group_means Output of `endpoint_volumes()$group_means`.
#' @param reference_group Control arm name.
#' @return `group_means` with a `TGI` column added.
#' @noRd
#' @keywords internal
endpoint_tgi <- function(group_means, reference_group) {
  ctrl <- group_means$Mean_Volume[group_means$Treatment == reference_group]
  if (length(ctrl) != 1L || !is.finite(ctrl) || ctrl <= 0) {
    stop("Reference group '", reference_group, "' has no usable endpoint ",
         "volume, so TGI is undefined.", call. = FALSE)
  }
  group_means$TGI <- (1 - group_means$Mean_Volume / ctrl) * 100
  group_means$TGI[group_means$Treatment == reference_group] <- 0
  group_means
}

#' Per-animal baseline at each animal's own first observation
#'
#' CODE_REVIEW.md R15.2. Four toxicity functions computed a baseline as
#' `data[data$Day == min(data$Day), ]` — the *global* earliest study day. Any
#' animal without an observation on that exact day is then either dropped by the
#' subsequent merge (therapeutic window) or carries an `NA` baseline that
#' propagates through every percentage derived from it.
#'
#' That is not an edge case. Staggered enrolment, a missed first weighing, or an
#' animal added after the study opened all produce it, and the bias has a
#' direction: the excluded animals are removed from the toxicity denominator, so
#' the reported weight loss is **understated** and the therapeutic window looks
#' better than it is. A worked example dropped the two most-toxic animals in an
#' arm and reported 10.0 % mean loss where the true figure across all six was
#' 16.7 %.
#'
#' Each animal's own first observation is the right baseline anyway: percentage
#' weight change is a within-animal quantity, and anchoring it to a day the animal
#' was not measured on is meaningless even when the row happens to exist.
#'
#' @param df Long-format data.
#' @param key_cols Character vector identifying an animal (e.g. `"MouseKey"`, or
#'   `c("ID", "Treatment")`).
#' @param value_col Name of the measurement column.
#' @param day_col Name of the time column.
#' @param out_name Name for the returned baseline column.
#' @return A data frame of `key_cols` plus `out_name`, one row per animal.
#'   Duplicate rows on an animal's first day are averaged, matching the previous
#'   behaviour.
#' @noRd
#' @keywords internal
me_per_mouse_baseline <- function(df, key_cols, value_col, day_col = "Day",
                                  out_name = "Baseline_Weight") {
  keep <- stats::complete.cases(df[, c(key_cols, day_col), drop = FALSE])
  d <- df[keep, , drop = FALSE]
  if (!nrow(d)) {
    out <- d[, key_cols, drop = FALSE]
    out[[out_name]] <- numeric(0)
    return(out)
  }
  key <- do.call(paste, c(lapply(key_cols, function(k) as.character(d[[k]])),
                          list(sep = "\r")))
  first_day <- stats::ave(as.numeric(d[[day_col]]), key,
                          FUN = function(x) min(x, na.rm = TRUE))
  at_first <- d[as.numeric(d[[day_col]]) == first_day, , drop = FALSE]

  fml <- stats::as.formula(paste(value_col, "~", paste(key_cols, collapse = " + ")))
  out <- stats::aggregate(fml, data = at_first, FUN = mean, na.rm = TRUE)
  names(out)[ncol(out)] <- out_name
  out
}
