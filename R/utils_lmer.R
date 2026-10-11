# A failed gradient check, confirmed or dismissed -- CODE_REVIEW.md R20-T
#
# lme4 checks the gradient at the optimum and warns "Model failed to converge
# with max|grad| = ..." when it exceeds 0.002. The check is often a false
# alarm (?lme4::convergence), and whether it fires moves with the lme4
# version, the platform and trivial changes to the data. On the shipped data
# the body-weight model warned; the warning sent it to a random intercept, and
# the Drug_B rate difference went from -0.0111 to -0.0068 g/day. Adding 1.5 g
# to one arm, which cannot change the slopes, warned in CI and not locally.
# Every one of those fits had reached the same optimum (log-likelihood
# -451.2811, the same slopes).
#
# lme4's documented checks are a restart from the estimates and a fit with a
# different optimizer. A clean restart, or a clean fit with the other
# optimizer, is used. When the other optimizer warns too but reaches the same
# optimum, the warning was a false alarm. Otherwise the caller's fallback
# proceeds as before, and the note says why.

#' Confirm or dismiss an lmer gradient warning
#'
#' @param fit A fitted `merMod`, or `NULL`.
#' @param msgs Warnings recorded while fitting it.
#' @param env Where the fit's call can be evaluated again (the caller's frame).
#' @return `list(fit, msgs, note)`. When the warning is dismissed, `msgs`
#'   loses it; `note` says how it was checked, or why it was upheld. With no
#'   gradient warning, `note` is NULL.
#' @noRd
#' @keywords internal
me_lmer_restart <- function(fit, msgs, env = parent.frame()) {
  grad_warning <- grepl("failed to converge with max|grad|", msgs, fixed = TRUE)
  if (is.null(fit) || !any(grad_warning)) return(list(fit = fit, msgs = msgs, note = NULL))

  refit_quietly <- function(cl) {
    w <- character()
    m <- tryCatch(
      withCallingHandlers(eval(cl, env),
        warning = function(x) {
          w <<- c(w, conditionMessage(x))
          invokeRestart("muffleWarning")
        },
        message = function(m) invokeRestart("muffleMessage")),
      error = function(e) NULL)
    list(fit = m, msgs = w)
  }
  clean <- function(r) !is.null(r$fit) && !any(grepl("converge", r$msgs, ignore.case = TRUE))
  # Values, not expressions: update() would insert getME(fit, ...), which the
  # caller's frame cannot evaluate.
  base_call <- stats::update(fit, evaluate = FALSE)

  # 1. Restart from the estimates.
  cl <- base_call
  cl$start <- list(theta = lme4::getME(fit, "theta"))
  r <- refit_quietly(cl)
  if (clean(r)) {
    return(list(fit = r$fit, msgs = c(msgs[!grad_warning], r$msgs),
                note = "lme4's gradient check failed; a restart from the estimates converged."))
  }

  # 2. Another optimizer, from lme4's default start. If it converges cleanly,
  #    its fit is used; if it reaches the same optimum as the first fit, the
  #    first warning was a false alarm.
  ctrl <- if (is.null(base_call$control)) lme4::lmerControl() else eval(base_call$control, env)
  ctrl$optimizer <- if (identical(fit@optinfo$optimizer, "bobyqa")) "nloptwrap" else "bobyqa"
  cl <- base_call
  cl$control <- ctrl
  r <- refit_quietly(cl)
  if (clean(r)) {
    return(list(fit = r$fit, msgs = c(msgs[!grad_warning], r$msgs), note = paste0(
      "lme4's gradient check failed; a fit with the ", ctrl$optimizer,
      " optimizer converged and is used.")))
  }
  # Practically the same optimum (lme4's allFit() guidance): log-likelihoods
  # within 0.01 and fixed effects within 0.1 SE. Tighter limits (0.001, 0.01
  # SE) still rejected the shipped data's two fits on Linux. The fit with the
  # higher log-likelihood is kept.
  if (!is.null(r$fit)) {
    se    <- sqrt(diag(as.matrix(stats::vcov(fit))))
    d_ll  <- as.numeric(stats::logLik(r$fit)) - as.numeric(stats::logLik(fit))
    d_b   <- max(abs(lme4::fixef(r$fit) - lme4::fixef(fit)) / se)
    if (abs(d_ll) < 0.01 && d_b <= 0.1) {
      return(list(fit = if (d_ll > 0) r$fit else fit, msgs = msgs[!grad_warning],
                  note = paste0(
        "lme4's gradient check failed, but a fit with the ", ctrl$optimizer,
        " optimizer reached the same optimum, so the warning was treated as a ",
        "false alarm (?lme4::convergence).")))
    }
    return(list(fit = fit, msgs = msgs, note = sprintf(paste0(
      "lme4's gradient check failed and was upheld: a restart warned again, and ",
      "the %s optimizer reached a different optimum (log-likelihood %+.4f, ",
      "fixed effects up to %.2f SE apart)."), ctrl$optimizer, d_ll, d_b)))
  }
  list(fit = fit, msgs = msgs, note = paste0(
    "lme4's gradient check failed and was upheld: a restart warned again, and ",
    "the ", ctrl$optimizer, " optimizer could not fit the model."))
}

# Starting values from a fit on centred time -- CODE_REVIEW.md R20-W
#
# With correlated random slopes on uncentred time the intercept and slope are
# strongly correlated, and the optimiser could stop short: counting the same
# data from implantation instead of randomisation (+14 days) gave
# log-likelihood 366.45 against 365.93, which moved the standard errors. The
# model is fitted first with time centred at `shift`, where it is well
# conditioned. Its random-effect covariance is mapped back to the original
# origin -- b0 = b0' - shift * b1, so Sigma = A Sigma' A' with
# A = [[1, -shift], [0, 1]] -- and returned as the starting theta for the
# model as written, which therefore starts at the optimum. The model and
# everything reported from it are unchanged; only where the optimiser starts.

#' Starting theta for a correlated random-slope model, from centred time
#'
#' @param formula The model formula (time enters as `time_column`).
#' @param data The data.
#' @param time_column,id_column The time variable and the animal grouping.
#' @param shift Where time is centred (the mean study day).
#' @param control An `lmerControl()`.
#' @return A theta vector for `lmer(start = list(theta = ...))`, or NULL.
#' @noRd
#' @keywords internal
me_centred_start <- function(formula, data, time_column, id_column, shift, control) {
  if (!is.finite(shift) || shift == 0) return(NULL)
  dc <- data
  dc[[time_column]] <- dc[[time_column]] - shift
  fc <- tryCatch(suppressWarnings(suppressMessages(
    lme4::lmer(formula, data = dc, control = control))), error = function(e) NULL)
  if (is.null(fc)) return(NULL)
  th <- lme4::getME(fc, "theta")
  S  <- tryCatch(as.matrix(lme4::VarCorr(fc)[[id_column]]), error = function(e) NULL)
  if (is.null(S) || !identical(rownames(S), c("(Intercept)", time_column))) return(NULL)
  A  <- matrix(c(1, 0, -shift, 1), 2L)
  L  <- tryCatch(t(chol(A %*% S %*% t(A))) / stats::sigma(fc), error = function(e) NULL)
  if (is.null(L)) return(NULL)
  idx <- startsWith(names(th), paste0(id_column, "."))
  if (sum(idx) != 3L) return(NULL)
  th[idx] <- L[lower.tri(L, diag = TRUE)]
  unname(th)
}
