# =============================================================================
# How early deaths move the model-based AUC (CODE_REVIEW.md R20-X addendum)
#
# Two arms of 10, log-linear growth with per-animal intercepts and slopes,
# measured every 3 days to day 30. Each study is generated once, and every
# dropout mechanism is applied to the same treated-arm animals, so the shift
# each mechanism causes is measured on paired data. The truth is the area under
# the treated arm's population geometric-mean curve over the AUC window.
#
# Run from the package root:  Rscript dev/sim_auc_dropout.R  (about 10 minutes)
# Results on 2026-10-22 (200 studies) are in CODE_REVIEW.md R20-X.
# =============================================================================

suppressMessages(devtools::load_all(".", quiet = TRUE))
days <- seq(0, 30, by = 3); n <- 10; a <- log(100)
b <- c(Control = 0.12, Drug = 0.07); sd_u <- 0.25; sd_v <- 0.03; sd_e <- 0.15; L <- 1000
true_auc <- function(b, T) (exp(a + b * T) - exp(a)) / b      # area under the population geometric-mean curve
sim_arm <- function(arm) {
  u <- rnorm(n, 0, sd_u); v <- rnorm(n, 0, sd_v)
  d <- expand.grid(i = 1:n, Day = days)
  d$logV <- a + u[d$i] + (b[[arm]] + v[d$i]) * d$Day + rnorm(nrow(d), 0, sd_e)
  d$Treatment <- arm; d$ID <- paste0(arm, d$i); d$v <- v[d$i]
  d[order(d$i, d$Day), ]
}
# Last day each drug animal is measured, by scenario.
last_day <- function(x, scen) switch(scen,
  none      = max(x$Day),
  # removed at the volume limit, after the measurement that crossed it (MAR)
  mar_limit = { k <- which(exp(x$logV) > L); if (length(k)) x$Day[k[1]] else max(x$Day) },
  # 30 % die at random, unrelated to the tumour (MCAR)
  mcar      = if (runif(1) < 0.3) sample(days[days >= 6 & days <= 27], 1) else max(x$Day),
  # the most aggressive 30 % die on days 9-18 whatever their volume (latent growth rate)
  mnar_aggr = if (x$v[1] > qnorm(0.7, 0, sd_v)) sample(days[days >= 9 & days <= 18], 1) else max(x$Day),
  # die just before the visit at which they would cross the limit (never measured large)
  mnar_next = { k <- which(exp(x$logV) > L); if (length(k)) x$Day[max(1, k[1] - 1)] else max(x$Day) })
fit_auc <- function(dd) {
  res <- tryCatch(suppressWarnings(suppressMessages(tumor_growth_statistics(dd[, c("ID", "Treatment", "Day", "Volume")],
           cage_column = NULL, model_type = "auc", p_adjust_method = "holm", reference_group = "Control",
           plots = FALSE, include_diagnostics = FALSE, return_model = FALSE))), error = function(e) NULL)
  if (is.null(res)) return(NULL)
  T <- res$auc_analysis$window[["end"]]; mb <- res$auc_analysis$model_based
  list(T = T, model = mb$AUC[mb$Treatment == "Drug"])
}
naive_auc <- function(drg, T) {
  g <- aggregate(logV ~ Day, data = drg[drg$Day <= T, ], FUN = mean)
  sum(diff(g$Day) * (head(exp(g$logV), -1) + tail(exp(g$logV), -1)) / 2)
}
scens <- c("none", "mcar", "mar_limit", "mnar_next", "mnar_aggr")
# One study, every scenario applied to the same animals.
one_paired <- function() {
  ctl <- sim_arm("Control"); full <- sim_arm("Drug")
  # the random draws of each scenario's death times are made once per animal
  out <- list()
  for (s in scens) {
    drg <- do.call(rbind, lapply(split(full, full$ID), function(x) x[x$Day <= last_day(x, s), ]))
    dd <- rbind(ctl, drg); dd$Volume <- exp(dd$logV)
    f <- fit_auc(dd); if (is.null(f)) return(NULL)
    tr <- true_auc(b[["Drug"]], f$T)
    out[[s]] <- c(lost = mean(tapply(drg$Day, drg$ID, max) < 30), T = f$T,
                  model = f$model / tr - 1, naive = naive_auc(drg, f$T) / tr - 1)
  }
  out
}
set.seed(20261023)
R <- 200
runs <- Filter(Negate(is.null), replicate(R, one_paired(), simplify = FALSE))
tab <- do.call(rbind, lapply(scens, function(s) {
  m <- do.call(rbind, lapply(runs, `[[`, s)); m0 <- do.call(rbind, lapply(runs, `[[`, "none"))
  dm <- m[, "model"] - m0[, "model"]; dn <- m[, "naive"] - m0[, "naive"]
  data.frame(scenario = s, studies = nrow(m), lost = sprintf("%.0f%%", 100 * mean(m[, "lost"])),
             model_vs_truth = sprintf("%+.1f%%", 100 * mean(m[, "model"])),
             naive_vs_truth = sprintf("%+.1f%%", 100 * mean(m[, "naive"])),
             model_shift_from_deaths = sprintf("%+.1f%% (SE %.1f)", 100 * mean(dm), 100 * sd(dm) / sqrt(length(dm))),
             naive_shift_from_deaths = sprintf("%+.1f%% (SE %.1f)", 100 * mean(dn), 100 * sd(dn) / sqrt(length(dn))))
}))
print(tab, row.names = FALSE)
