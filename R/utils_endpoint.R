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
#' @param endpoint_day Day at which to evaluate. `NULL` uses the last day on
#'   which every arm in `arms` is evaluable (see [evaluable_days()]); an
#'   explicit day must be evaluable (R20.1, R20-K).
#' @param endpoint_method One of:
#'   \describe{
#'     \item{`"model"`}{(default) each arm's geometric mean at `endpoint_day`
#'       from `me_endpoint_model()`, fitted to every observation of every
#'       animal. Its intervals come from draws of that model's fixed effects.}
#'     \item{`"last_obs"`}{each animal's own last observation at or before the
#'       endpoint day. Animals are evaluated at *different days*; on a simulated
#'       study with volume-triggered euthanasia this was *more* biased than
#'       `"survivors"`. Retained as a fallback for when the model cannot be
#'       fitted, not as an equal alternative.}
#'     \item{`"survivors"`}{raw mean among animals observed at `endpoint_day`,
#'       the pre-0.8.0 behaviour, retained for reproducibility. Warns when
#'       animals were lost.}
#'   }
#' @param arms Arms that must be evaluable, and whose means are returned.
#'   `NULL` means every arm in `df`.
#' @param model A fitted `me_endpoint_model()` object to reuse, e.g. across
#'   the days of an over-time analysis. `NULL` fits one.
#' @return A list with `group_means` (data frame: Treatment, Mean_Volume, N,
#'   N_On_Study, plus SE_log / CI bounds under `"model"`), `per_mouse`
#'   (Treatment, MouseKey, Volume, Day_Used) for the resampling-based methods,
#'   `attrition`, `endpoint_day`, `method`, `evaluability` (see
#'   [evaluable_days()]) and, under `"model"`, `model`.
#' @noRd
#' @keywords internal
endpoint_volumes <- function(df,
                             id_column = "ID",
                             treatment_column = "Treatment",
                             time_column = "Day",
                             volume_column = "Volume",
                             cage_column = NULL,
                             endpoint_day = NULL,
                             endpoint_method = c("model", "last_obs", "survivors"),
                             arms = NULL,
                             model = NULL) {
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

  # R20.1 / R20-K: evaluate only where every arm has enough animals on study.
  ev     <- me_evaluability(d, arms)
  ep_day <- me_resolve_eval_day(ev, endpoint_day)
  arms   <- ev$arms

  # Attrition bookkeeping -- the numbers that make survivor selection visible.
  da         <- d[d$Treatment %in% arms, , drop = FALSE]
  all_mice   <- unique(da[, c("MouseKey", "Treatment")])
  at_risk    <- unique(da[da$Day >= ep_day, c("MouseKey", "Treatment")])
  n_total    <- table(factor(all_mice$Treatment, levels = arms))
  n_at_risk  <- table(factor(at_risk$Treatment, levels = arms))
  attrition  <- data.frame(
    Treatment     = arms,
    N_Enrolled    = as.integer(n_total),
    N_At_Endpoint = as.integer(n_at_risk),
    Pct_Lost      = round(100 * (1 - as.integer(n_at_risk) / as.integer(n_total)), 1),
    stringsAsFactors = FALSE
  )

  if (endpoint_method == "survivors" &&
      any(attrition$N_At_Endpoint < attrition$N_Enrolled)) {
    lost <- attrition[attrition$N_At_Endpoint < attrition$N_Enrolled, ]
    msg <- paste(sprintf("%s: %d/%d", lost$Treatment,
                         lost$N_At_Endpoint, lost$N_Enrolled), collapse = "; ")
    warning("endpoint_method = 'survivors' conditions on being observed at ",
            "day ", ep_day, ", but animals were lost before then (", msg,
            "). Animals leave because their tumours grew, so this selects the ",
            "slowest growers -- most severely in the control arm -- and biases ",
            "TGI downward. Use endpoint_method = 'model' (default) or ",
            "'last_obs'.", call. = FALSE)
  }

  per_mouse <- NULL
  group_means <- NULL

  if (endpoint_method == "model") {
    em <- if (!is.null(model)) model else me_endpoint_model(d)
    if (is.null(em)) {
      warning("Model-based endpoint means could not be fitted; falling back to ",
              "endpoint_method = 'last_obs'.", call. = FALSE)
      endpoint_method <- "last_obs"
    } else {
      lm_ <- me_endpoint_logmeans(em, arms, ep_day)
      z <- stats::qnorm(0.975)
      group_means <- data.frame(
        Treatment   = arms,
        # exp() of a log-scale mean is a geometric mean -- the right centre
        # for a log-normal quantity, and consistent with the modelling scale.
        Mean_Volume = exp(lm_$log_mean),
        SE_log      = lm_$se_log,
        Lower_CL    = exp(lm_$log_mean - z * lm_$se_log),
        Upper_CL    = exp(lm_$log_mean + z * lm_$se_log),
        N           = attrition$N_Enrolled,
        N_On_Study  = attrition$N_At_Endpoint,
        stringsAsFactors = FALSE
      )
    }
  }

  if (endpoint_method %in% c("last_obs", "survivors")) {
    em <- NULL
    per_mouse <- if (endpoint_method == "last_obs") {
      do.call(rbind, lapply(split(da, da$MouseKey, drop = TRUE), function(s) {
        s <- s[order(s$Day), ]
        # Each animal's last observation at or before the endpoint day.
        s <- s[s$Day <= ep_day, , drop = FALSE]
        if (nrow(s) == 0L) return(NULL)
        data.frame(MouseKey = s$MouseKey[1], Treatment = s$Treatment[1],
                   Volume = s$Volume[nrow(s)], Day_Used = s$Day[nrow(s)],
                   stringsAsFactors = FALSE)
      }))
    } else {
      sv <- da[da$Day == ep_day, , drop = FALSE]
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
    group_means$N_On_Study <- attrition$N_At_Endpoint[
      match(group_means$Treatment, attrition$Treatment)]
  }

  list(group_means = group_means, per_mouse = per_mouse,
       attrition = attrition, endpoint_day = ep_day,
       method = endpoint_method, evaluability = ev,
       model = if (endpoint_method == "model") em else NULL)
}

#' The endpoint model: log volume on a per-arm spline in time, random slopes
#'
#' CODE_REVIEW.md R20.1 / R20.29. The previous model,
#' `log(V) ~ Treatment * Day + (1 | animal)`, forced log volume to be linear in
#' time within each arm and gave every animal the arm's slope. Real xenograft
#' growth decelerates on the log scale, so the straight line overshot wherever
#' an arm's data thinned out; on the Combo demo it put the control at
#' 28,544 mm3 on day 32 against about 3,000 observed. Leaving out per-animal
#' slopes also biased TGI by about 6 points under informative dropout.
#'
#' This model:
#' - gives each arm its own natural spline in time (3 degrees of freedom),
#'   or a straight line when an arm has fewer than 4 measured days;
#' - uses correlated per-animal random slopes, falling back to uncorrelated
#'   slopes and then to a random intercept only when a fit fails (R20-K);
#' - leaves out zero volumes measured before an animal's first positive
#'   volume (pre-palpable tumours): flooring them at half the smallest volume
#'   invented values far below the data (R20.1). A zero after a positive
#'   measurement (a regression) is set to the smallest positive volume in
#'   the study, the detection limit.
#'
#' It is evaluated only on evaluable days (see [evaluable_days()]).
#'
#' @param d Data frame with `MouseKey`, `Treatment`, `Day` and `Volume`.
#' @return An `me_endpoint_model` list, or `NULL` when no model can be fitted.
#' @noRd
#' @keywords internal
me_endpoint_model <- function(d) {
  d <- d[is.finite(d$Day) & is.finite(d$Volume), , drop = FALSE]
  pos <- d$Volume > 0
  if (!any(pos)) return(NULL)
  det_limit <- min(d$Volume[pos])
  first_pos <- stats::ave(ifelse(pos, d$Day, Inf), d$MouseKey, FUN = min)
  prepalp   <- !pos & d$Day < first_pos
  fd <- d[!prepalp, , drop = FALSE]
  n_floored <- sum(fd$Volume <= 0)
  fd$.logv <- log(pmax(fd$Volume, det_limit))

  lev <- sort(unique(as.character(d$Treatment)))
  fd$Treatment <- factor(as.character(fd$Treatment), levels = lev)
  if (any(table(fd$Treatment) == 0L)) return(NULL)
  if (length(unique(fd$Day)) < 2L || nrow(fd) < 6L) return(NULL)

  n_days <- tapply(fd$Day, fd$Treatment, function(x) length(unique(x)))
  spline <- all(n_days >= 4L)
  if (spline) {
    basis <- splines::ns(fd$Day, df = 3L)
    B <- matrix(as.numeric(basis), nrow = nrow(fd))
  } else {
    basis <- NULL
    B <- matrix(fd$Day, ncol = 1L)
  }
  cols <- paste0(".tb", seq_len(ncol(B)))
  colnames(B) <- cols
  fd <- cbind(fd, as.data.frame(B))
  day_mean <- mean(fd$Day)
  fd$.day_c <- fd$Day - day_mean

  rhs_fixed <- paste("Treatment * (", paste(cols, collapse = " + "), ")")
  re_terms <- c(correlated   = "(.day_c | MouseKey)",
                uncorrelated = "(.day_c || MouseKey)",
                intercept    = "(1 | MouseKey)")
  ctrl <- lme4::lmerControl(check.nobs.vs.nlev = "ignore",
                            check.nobs.vs.nRE  = "ignore")
  fit <- NULL; used <- NA_character_; fallback <- character(0)
  for (nm in names(re_terms)) {
    msgs <- character(0)
    f <- tryCatch(
      withCallingHandlers(
        lme4::lmer(stats::as.formula(paste(".logv ~", rhs_fixed, "+", re_terms[[nm]])),
                   data = fd, REML = TRUE, control = ctrl),
        warning = function(w) {
          msgs <<- c(msgs, conditionMessage(w))
          invokeRestart("muffleWarning")
        },
        message = function(m) invokeRestart("muffleMessage")),
      error = function(e) {
        msgs <<- c(msgs, conditionMessage(e))
        NULL
      })
    failed <- is.null(f) ||
      any(grepl("failed to converge|unable to evaluate|unidentifiable|degenerate",
                msgs, ignore.case = TRUE))
    if (!failed) {
      fit <- f; used <- nm
      break
    }
    fallback <- c(fallback, sprintf("%s random effects: %s", nm,
                                    if (length(msgs)) msgs[1L] else "failed"))
  }
  if (is.null(fit)) return(NULL)

  beta <- lme4::fixef(fit)
  structure(list(
    fit       = fit,
    rhs       = stats::as.formula(paste("~", rhs_fixed)),
    basis     = basis,
    spline    = spline,
    cols      = cols,
    levels    = lev,
    beta      = beta,
    V         = as.matrix(stats::vcov(fit))[names(beta), names(beta), drop = FALSE],
    structure = used,
    fallback  = fallback,
    detection_limit = det_limit,
    n_prepalpable   = sum(prepalp),
    n_floored       = n_floored,
    day_range       = range(fd$Day),
    day_mean        = day_mean
  ), class = "me_endpoint_model")
}

#' Fixed-effect design rows for one arm at the given days
#' @noRd
#' @keywords internal
me_endpoint_X <- function(em, arm, t) {
  nd <- data.frame(Treatment = factor(rep(arm, length(t)), levels = em$levels))
  B <- if (em$spline) {
    matrix(as.numeric(stats::predict(em$basis, t)), nrow = length(t))
  } else {
    matrix(t, ncol = 1L)
  }
  colnames(B) <- em$cols
  nd <- cbind(nd, as.data.frame(B))
  X <- stats::model.matrix(em$rhs, nd)
  X[, names(em$beta), drop = FALSE]
}

#' Each arm's log-scale mean and standard error at one day
#' @noRd
#' @keywords internal
me_endpoint_logmeans <- function(em, arms, t) {
  X <- do.call(rbind, lapply(arms, function(a) me_endpoint_X(em, a, t)))
  data.frame(
    Treatment = arms,
    log_mean  = as.numeric(X %*% em$beta),
    se_log    = sqrt(pmax(rowSums((X %*% em$V) * X), 0)),
    stringsAsFactors = FALSE)
}

#' Draws of the endpoint model's fixed effects
#'
#' Parametric draws from the estimated sampling distribution of the fixed
#' effects, beta ~ N(beta-hat, V). Any function of the arms' means, such as
#' TGI, the Bliss excess or an AUC, gets its interval by evaluating it on each
#' draw. The interval then describes the reported estimate itself, which the
#' previous bootstrap of last observations did not (R20.2), and needs no
#' refits.
#'
#' @param em An `me_endpoint_model`.
#' @param n Number of draws.
#' @param seed Optional seed; the caller's RNG state is restored.
#' @return An n x p matrix.
#' @noRd
#' @keywords internal
me_beta_draws <- function(em, n, seed = NULL) {
  me_with_seed(seed, {
    p <- length(em$beta)
    L <- tryCatch(chol(em$V), error = function(e) NULL)
    if (is.null(L)) {
      eg <- eigen(em$V, symmetric = TRUE)
      L <- t(eg$vectors %*% diag(sqrt(pmax(eg$values, 0)), p))
    }
    Z <- matrix(stats::rnorm(n * p), n, p)
    sweep(Z %*% L, 2L, em$beta, "+")
  })
}

#' Evaluate an expression under a seed, restoring the caller's RNG state
#' @noRd
#' @keywords internal
me_with_seed <- function(seed, expr) {
  if (!is.null(seed)) {
    old <- if (exists(".Random.seed", envir = .GlobalEnv)) {
      get(".Random.seed", envir = .GlobalEnv)
    } else NULL
    on.exit({
      if (!is.null(old)) assign(".Random.seed", old, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }
  expr
}

#' Log-scale arm means on each draw
#'
#' @return A draws x arms matrix of log-scale means at day `t`.
#' @noRd
#' @keywords internal
me_draw_logmeans <- function(em, draws, arms, t) {
  X <- do.call(rbind, lapply(arms, function(a) me_endpoint_X(em, a, t)))
  out <- draws %*% t(X)
  colnames(out) <- arms
  out
}

#' Summary of how the endpoint model was fitted, for results and displays
#' @noRd
#' @keywords internal
me_endpoint_model_info <- function(em) {
  if (is.null(em)) return(NULL)
  list(
    formula = paste0("log(volume) ~ Treatment x ",
                     if (em$spline) "natural spline in day (3 df)" else "day",
                     " + ",
                     switch(em$structure,
                            correlated   = "(day | animal)",
                            uncorrelated = "(day || animal)",
                            intercept    = "(1 | animal)")),
    time_basis      = if (em$spline) "natural spline, 3 df" else "linear",
    random_effects  = em$structure,
    fallback        = em$fallback,
    detection_limit = em$detection_limit,
    n_prepalpable_excluded = em$n_prepalpable,
    n_zero_after_positive  = em$n_floored
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
