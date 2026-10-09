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
# optimum, the warning was a false alarm and the first fit is kept. Otherwise
# the caller's fallback proceeds as before.

#' Confirm or dismiss an lmer gradient warning
#'
#' @param fit A fitted `merMod`, or `NULL`.
#' @param msgs Warnings recorded while fitting it.
#' @param env Where the fit's call can be evaluated again (the caller's frame).
#' @return `list(fit, msgs, note)`. When the warning is dismissed, `msgs`
#'   loses it and `note` says how it was checked; otherwise `note` is NULL.
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
  if (!is.null(r$fit)) {
    se <- sqrt(diag(as.matrix(stats::vcov(fit))))
    same_ll   <- abs(as.numeric(stats::logLik(r$fit)) - as.numeric(stats::logLik(fit))) < 1e-3
    same_beta <- all(abs(lme4::fixef(r$fit) - lme4::fixef(fit)) <= 0.01 * se)
    if (same_ll && same_beta) {
      return(list(fit = fit, msgs = msgs[!grad_warning], note = paste0(
        "lme4's gradient check failed, but a fit with the ", ctrl$optimizer,
        " optimizer reached the same optimum, so the warning was treated as a ",
        "false alarm (?lme4::convergence).")))
    }
  }
  list(fit = fit, msgs = msgs, note = NULL)
}
