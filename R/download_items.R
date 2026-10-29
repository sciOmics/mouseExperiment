# =============================================================================
# What each analysis tab's "Download Results" offers (CODE_REVIEW.md R20-D.X,
# R20-D.Y). Each function takes the module's stored result and returns a list
# of dl_item()s for what that result holds. Plots are passed in as functions
# and drawn only when the ZIP is written (see R/download_results.R).
# =============================================================================

dl_has <- function(x) !is.null(x) && (!is.data.frame(x) || nrow(x) > 0L)
dl_get <- function(x) function() x
dl_drawn <- function(f) function() tryCatch(f(), error = function(e) NULL)

#' Survival
#' @param res `rv$surv_results`.
#' @param km_plot,forest_plot Functions returning the Kaplan-Meier and forest
#'   plots as displayed.
#' @param posterior_plot A function returning the posterior distributions as
#'   displayed (R23.12), or NULL.
#' @keywords internal
surv_dl_items <- function(res, km_plot, forest_plot, posterior_plot = NULL) {
  if (is.null(res) || isTRUE(res$error)) return(NULL)
  bayes <- identical(res$model_type_used, "bayes")
  aft   <- identical(res$model_type_used, "aft")
  ph    <- res$ph_test
  draws <- res$posterior_draws
  Filter(Negate(is.null), list(
    if (!is.null(res$call_record)) script_dl_item(function(up) surv_script_bundle(res, up)),
    if (!is.null(surv_effects(res)))   # R22.2: a Bayesian run's are treatment_effects
      dl_item("results", "tables",
              if (bayes) "Treatment effects" else if (aft) "Time ratios" else "Hazard ratios",
              "survival_results", "table", function() surv_effects(res)$table,
              if (bayes) "Posterior time ratios by arm, with credible intervals"
              else if (aft) "Each arm's time ratio against the reference, with the correction selected"
              else "Each arm against the reference, with the correction selected"),
    # R23.11: the AFT model's checks.
    if (aft && dl_has(res$distribution_comparison))
      dl_item("aft_dists", "diagnostics", "AFT distributions (CSV)", "aft_distributions", "table",
              dl_get(res$distribution_comparison), "Each distribution's AIC, and the one used"),
    if (aft && dl_has(res$fitted_curves))
      dl_item("aft_fit", "diagnostics", "AFT fit against Kaplan-Meier", "aft_fit", "plot",
              function() aft_fit_plot(res), "Kaplan-Meier curves with the fitted AFT curves"),
    if (!is.null(res$km_fit))
      dl_item("risk", "tables", "Risk table", "risk_table", "table",
              function() surv_risk_table(res$km_fit),
              "Animals at risk, events and censored, by arm and day"),
    dl_item("plot_km", "plots", "Kaplan-Meier curves", "kaplan_meier", "plot",
            dl_drawn(km_plot), "Survival curves by arm, as displayed"),
    if (!is.null(surv_effects(res)))
      dl_item("plot_forest", "plots", "Forest plot", "forest_plot", "plot",
              dl_drawn(forest_plot),
              if (bayes) "Time ratios with their credible intervals"
              else if (aft) "Time ratios with their intervals" else "Hazard ratios with their intervals"),
    if (!is.null(ph) && !is.null(ph$table))
      dl_item("ph_table", "diagnostics", "Proportional-hazards test (CSV)", "ph_test",
              "table", function() {
                tab <- as.data.frame(ph$table)
                cbind(Covariate = rownames(tab), tab, row.names = NULL)
              }, "Schoenfeld-residual test of proportional hazards (cox.zph)"),
    if (!is.null(ph))
      dl_item("ph_plot", "diagnostics", "Schoenfeld residuals", "schoenfeld_residuals",
              "plot", function() function() plot(ph),
              "Schoenfeld residuals over time, by covariate"),
    if (bayes && dl_has(res$mcmc_diagnostics))
      dl_item("diag_mcmc", "diagnostics", "MCMC diagnostics (CSV)", "mcmc_diagnostics",
              "table", dl_get(res$mcmc_diagnostics), "Rhat and effective sample sizes"),
    if (bayes)
      dl_item("diag_influential", "diagnostics", "Influential observations (CSV)",
              "influential_observations", "table",
              function() bayes_influential_obs(res),
              "Observations with Pareto-k above the PSIS-LOO threshold"),
    if (bayes && dl_has(draws))
      dl_item("diag_trace", "diagnostics", "MCMC trace plots", "mcmc_trace", "plot",
              function() bayesplot::mcmc_trace(draws), "Each chain's draws, after warmup"),
    if (dl_has(res$pp_check_plot))
      dl_item("diag_pp_check", "diagnostics", "Posterior predictive check", "pp_check",
              "plot", dl_get(res$pp_check_plot), "Observed data against posterior predictions"),
    if (bayes && dl_has(res$prior_posterior_data) && is.function(posterior_plot))
      dl_item("plot_posterior", "plots", "Posterior distributions", "posterior_distributions",
              "plot", dl_drawn(posterior_plot),
              "Each treatment effect's posterior, with its highest-density interval, as displayed"),
    if (dl_has(res$survival_data))
      dl_item("data", "data", "Survival data", "survival_data", "table",
              dl_get(res$survival_data), "The data as analysed: one row per animal")
  ))
}

#' Toxicity (body weight, or the therapeutic window)
#' @param res `rv$tox_results`.
#' @param traj_plot,km_plot,window_plot Functions returning the plots.
#' @keywords internal
tox_dl_items <- function(res, traj_plot, km_plot, window_plot) {
  if (is.null(res) || isTRUE(res$error)) return(NULL)
  if (identical(res$model_type, "et")) {
    twm <- res$twm
    if (is.null(twm) || isTRUE(twm$error)) return(NULL)
    return(Filter(Negate(is.null), list(
      if (!is.null(res$call_records)) script_dl_item(function(up) tox_script_bundle(res, up)),
      if (dl_has(twm$window_table))
        dl_item("window", "tables", "Therapeutic window", "therapeutic_window", "table",
                dl_get(twm$window_table), "TGI and worst weight loss by arm, with tolerability"),
      if (dl_has(twm$weight_loss_data))
        dl_item("window_animals", "tables", "Per-animal weight loss", "weight_loss_per_animal",
                "table", dl_get(twm$weight_loss_data), "Each animal's worst weight loss"),
      dl_item("plot_window", "plots", "Therapeutic window", "therapeutic_window", "plot",
              dl_drawn(window_plot), "TGI against worst weight loss, by arm")
    )))
  }
  mm <- res$mixed_model
  mm_ok <- !is.null(mm) && !isTRUE(mm$error)
  wl <- res$weight_loss
  wl_ok <- !is.null(wl) && !isTRUE(wl$error)
  Filter(Negate(is.null), list(
    if (!is.null(res$call_records) && (mm_ok || wl_ok))
      script_dl_item(function(up) tox_script_bundle(res, up)),
    if (mm_ok && dl_has(mm$fixed_effects))
      dl_item("fixed", "tables", "Model coefficients", "body_weight_model", "table",
              dl_get(mm$fixed_effects), "Fixed effects of the body-weight model"),
    if (mm_ok && !is.null(mm$emmeans_table))
      dl_item("rates", "tables", "Weight change per day", "weight_change_rates", "table",
              function() as.data.frame(mm$emmeans_table), "Each arm's rate of weight change"),
    if (mm_ok && !is.null(mm$pairwise_comparisons))
      dl_item("rate_comparisons", "tables", "Rate comparisons", "rate_comparisons", "table",
              function() as.data.frame(mm$pairwise_comparisons),
              "Each arm's rate against the reference's"),
    if (mm_ok && !is.null(mm$level_comparisons))
      dl_item("level_comparisons", "tables", "Weight at the mean day", "level_comparisons",
              "table", function() as.data.frame(mm$level_comparisons),
              "Arms compared on weight at the mean study day"),
    if (wl_ok && dl_has(wl$event_data))
      dl_item("events", "tables", "Weight-loss events", "weight_loss_events", "table",
              dl_get(wl$event_data), "Each animal's time to the weight-loss threshold"),
    if (dl_has(res$trajectory_data))
      dl_item("plot_traj", "plots", "Body-weight trajectories", "body_weight_trajectories",
              "plot", dl_drawn(traj_plot), "Body weight over time, as displayed"),
    if (wl_ok)
      dl_item("plot_km", "plots", "Time to weight loss", "weight_loss_km", "plot",
              dl_drawn(km_plot), "Kaplan-Meier curves of time to the threshold"),
    if (mm_ok && dl_has(mm$diag_qq_plot))
      dl_item("diag_qq", "diagnostics", "Residual Q-Q", "residual_qq", "plot",
              dl_get(mm$diag_qq_plot), "Q-Q plot of the standardised residuals"),
    if (mm_ok && dl_has(mm$diag_resid_fitted_plot))
      dl_item("diag_resid", "diagnostics", "Residuals vs fitted", "residuals_vs_fitted",
              "plot", dl_get(mm$diag_resid_fitted_plot), "Residuals against fitted values"),
    if (mm_ok && dl_has(mm$diag_scale_location_plot))
      dl_item("diag_scale", "diagnostics", "Scale-location", "scale_location", "plot",
              dl_get(mm$diag_scale_location_plot), "Spread of the residuals against fitted values"),
    if (mm_ok && dl_has(mm$diag_re_qq_plot))
      dl_item("diag_re_qq", "diagnostics", "Random-effects Q-Q", "random_effects_qq", "plot",
              dl_get(mm$diag_re_qq_plot), "Q-Q plot of the per-animal random effects"),
    if (mm_ok && dl_has(mm$diag_cooks_distance))
      dl_item("diag_cooks", "diagnostics", "Cook's distance (CSV)", "cooks_distance", "table",
              dl_get(mm$diag_cooks_distance), "Influence of each animal on the fit"),
    if (mm_ok && dl_has(mm$diag_dfbetas))
      dl_item("diag_dfbetas", "diagnostics", "DFBETAS (CSV)", "dfbetas", "table",
              dl_get(mm$diag_dfbetas), "Influence of each animal on each coefficient"),
    if (dl_has(res$trajectory_data))
      dl_item("data", "data", "Body-weight data", "body_weight_data", "table",
              dl_get(res$trajectory_data),
              "The data as analysed: one row per weighing, after any tumour-mass correction")
  ))
}

#' Drug Synergy
#' @param res `rv$syn_results`.
#' @param over_time Whether the run covered every evaluable day.
#' @param main_plot,trend_plot,tgi_plot Functions returning the plots last drawn.
#' @param diag The diagnostics' source (the displayed result, or the last day's).
#' @param plot_record A function returning the Bliss plot's record (F9), or NULL.
#' @keywords internal
syn_dl_items <- function(res, over_time, main_plot, trend_plot, tgi_plot, diag, plot_record = NULL) {
  if (is.null(res) || isTRUE(res$error)) return(NULL)
  metrics <- if (isTRUE(over_time) && dl_has(res$synergy_summary)) res$synergy_summary
             else res$synergy_metrics
  shown <- "Drawn with the Plots tab's current settings (R23.14)."
  Filter(Negate(is.null), list(
    if (!is.null(res$call_record))
      script_dl_item(function(up) syn_script_bundle(res, up, plot = if (is.function(plot_record)) plot_record())),
    if (!is.null(metrics))
      dl_item("metrics", "tables", "Synergy metrics", "synergy_metrics", "table",
              function() as.data.frame(metrics),
              if (isTRUE(over_time)) "Synergy metrics for each evaluable day"
              else "Bliss expectation, observed TGI and the difference"),
    if (is.data.frame(res$summary) && nrow(res$summary))
      dl_item("tgi", "tables", "TGI by arm", "tgi_by_arm", "table", dl_get(res$summary),
              "Each arm's tumour growth inhibition"),
    if (dl_has(res$synergy_ci))
      dl_item("ci", "tables", "Intervals", "synergy_intervals", "table",
              dl_get(res$synergy_ci), "95% intervals of the synergy metrics"),
    if (dl_has(res$attrition))
      dl_item("attrition", "tables", "Attrition", "attrition", "table",
              dl_get(res$attrition), "Animals on study by arm and day"),
    dl_item("plot_main", "plots", "Synergy plot", "synergy_plot", "plot",
            dl_drawn(main_plot), paste("Observed against expected TGI.", shown)),
    if (isTRUE(over_time))
      dl_item("plot_trend", "plots", "Synergy over time", "synergy_trend", "plot",
              dl_drawn(trend_plot), paste("The Bliss difference by day.", shown)),
    # The TGI plot is drawn for the analysis over every evaluable day only.
    if (isTRUE(over_time))
      dl_item("plot_tgi", "plots", "TGI by arm", "tgi_by_arm", "plot",
              dl_drawn(tgi_plot), paste("Each arm's TGI over time.", shown)),
    if (dl_has(diag$diag_group_qq_plot))
      dl_item("diag_qq", "diagnostics", "Residuals by arm (Q-Q)", "residuals_by_arm_qq",
              "plot", dl_get(diag$diag_group_qq_plot), "Q-Q plots of the residuals, by arm"),
    if (dl_has(diag$diag_fit_plot %||% diag$diag_group_boxplot))
      dl_item("diag_fit", "diagnostics",
              if (dl_has(diag$diag_fit_plot)) "Fitted curves" else "Volumes by arm",
              "fit_by_arm", "plot", dl_get(diag$diag_fit_plot %||% diag$diag_group_boxplot),
              "Each arm's data with the model's curve, or the volumes used")
  ))
}

#' Dose Response
#' @param res `rv$dr_results`.
#' @param scatter_plot,tgi_plot Functions returning the plots last drawn.
#' @keywords internal
dr_dl_items <- function(res, scatter_plot, tgi_plot) {
  if (is.null(res) || isTRUE(res$error)) return(NULL)
  st <- res$statistics
  shown <- "Drawn with the Plots tab's current settings (R23.14)."
  Filter(Negate(is.null), list(
    if (!is.null(res$call_record)) script_dl_item(function(up) dr_script_bundle(res, up)),
    if (dl_has(res$summary_table))
      dl_item("summary", "tables", "Summary by dose", "dose_summary", "table",
              function() as.data.frame(res$summary_table), "Tumour volume by dose group"),
    if (!is.null(st))
      dl_item("stats", "tables", "Dose-response statistics", "dose_response_statistics",
              "table", function() dr_stats_table(res),
              "Trend tests, EC50 with its interval, and the curve's parameters"),
    if (dl_has(res$tgi_table))
      dl_item("tgi", "tables", "TGI by dose", "tgi_by_dose", "table",
              function() as.data.frame(res$tgi_table), "Each dose's TGI with its interval"),
    if (!is.null(st$tukey_results))
      dl_item("tukey", "tables", "Pairwise comparisons", "pairwise_comparisons", "table",
              function() as.data.frame(st$tukey_results), "Comparisons between doses"),
    dl_item("plot_scatter", "plots", "Dose-response curve", "dose_response", "plot",
            dl_drawn(scatter_plot), paste("Volume against dose, with the fitted curve.", shown)),
    dl_item("plot_tgi", "plots", "TGI by dose", "tgi_by_dose", "plot",
            dl_drawn(tgi_plot), paste("Each dose's TGI.", shown)),
    if (!is.null(st$dr_lack_of_fit))
      dl_item("diag_lof", "diagnostics", "Lack-of-fit test (CSV)", "lack_of_fit", "table",
              function() as.data.frame(st$dr_lack_of_fit),
              "Whether the curve fits the dose-group means"),
    if (dl_has(st$dr_residuals_plot))
      dl_item("diag_resid", "diagnostics", "Residuals vs dose", "residuals_vs_dose", "plot",
              dl_get(st$dr_residuals_plot), "The curve's residuals against dose"),
    if (dl_has(res$analysis_data))
      dl_item("data", "data", "Analysed data", "analysed_data", "table",
              function() as.data.frame(res$analysis_data),
              "The data as analysed: one row per animal at the evaluation day")
  ))
}

#' Power Analysis
#' @param res `rv$ap_results`.
#' @param curve_plot Function returning the power-curve plot last drawn.
#' @keywords internal
pw_dl_items <- function(res, curve_plot) {
  if (is.null(res) || isTRUE(res$error)) return(NULL)
  lmm <- identical(res$mode, "apriori_lmm")
  Filter(Negate(is.null), list(
    if (length(res$call_records)) script_dl_item(function(up) pw_script_bundle(res)),
    if (dl_has(res$scenario_table))
      dl_item("scenarios", "tables", if (lmm) "Simulated power" else "Scenarios",
              if (lmm) "lmm_power" else "scenarios", "table",
              function() as.data.frame(res$scenario_table),
              if (lmm) "Power by N per group, from the simulations"
              else "Required N or achieved power for each scenario"),
    if (dl_has(res$sensitivity_table))
      dl_item("sensitivity", "tables", "SD sensitivity", "sd_sensitivity", "table",
              function() as.data.frame(res$sensitivity_table),
              "How the required N changes with the assumed SD"),
    if (dl_has(res$power_curve_data))
      dl_item("curve_data", "tables", "Power-curve data", "power_curve_data", "table",
              function() as.data.frame(res$power_curve_data), "The points behind the power curves"),
    dl_item("plot_curves", "plots", "Power curves", "power_curves", "plot",
            dl_drawn(curve_plot), "Power against N per group, drawn with the current settings.")
  ))
}
