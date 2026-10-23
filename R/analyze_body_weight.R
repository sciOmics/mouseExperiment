#' Analyze Body Weight via Mixed Modeling
#'
#' Fits a linear mixed-effects model to longitudinal body weight data,
#' optionally adjusting for estimated tumor weight, sex, initial mass,
#' and tumor volume as covariates.
#'
#' @param df Data frame with longitudinal data.
#' @param weight_column Name of the body weight column (grams).
#' @param time_column Name of the time/day column.
#' @param treatment_column Name of the treatment group column.
#' @param id_column Name of the mouse/subject ID column.
#' @param volume_column Name of the tumor volume column (mm³). NULL to skip tumor adjustment.
#' @param sex_column Name of the sex column. NULL to omit.
#' @param cage_column Name of the cage column. NULL to omit.
#' @param adjust_tumor_weight Logical; subtract estimated tumor weight from body mass.
#' @param volume_units Units of the volume column, \code{"mm3"} or \code{"cm3"}.
#'   Required when \code{adjust_tumor_weight = TRUE} and a volume column is used,
#'   because volume is converted to mass and subtracted from body weight: the
#'   wrong unit scales that correction by 1000. Units were inferred from the data
#'   before v0.25.0, which read small-tumour mm3 studies as cm3 (CODE_REVIEW.md
#'   R20.83); the data are now only checked against the declared unit, with a
#'   warning on disagreement and an error when the implied tumour mass exceeds
#'   half the body weight.
#' @param tumor_density Density in g/cm³ for tumor weight estimation (default 1.0).
#' @param covariates Character vector of optional covariates: "volume", "sex", "initial_mass".
#' @param estimation Character; "REML" (default) or "ML".
#' @param comparison_family Which comparisons to report and adjust over:
#'   "vs_reference" (default), "all_pairs", or "custom" (with
#'   \code{custom_contrasts}). The multiplicity adjustment covers exactly this
#'   family.
#' @param custom_contrasts Named list of contrast coefficient vectors over the
#'   treatment levels; required when \code{comparison_family = "custom"}.
#' @param p_adjust_method Multiplicity adjustment: "bonferroni" (default),
#'   "holm", "fdr", "dunnett" (exact many-to-one, requires "vs_reference"),
#'   "tukey" (requires "all_pairs"), or "none".
#' @param reference_group Name of the control/reference group. NULL auto-selects.
#' @return A list with components: model, fixed_effects, random_effects,
#'   emmeans_table (group marginal means at the mean study day),
#'   pairwise_comparisons (the requested contrast family with adjusted
#'   p-values), comparison_family, p_adjust_method_used, model_info,
#'   weight_data, summary_text.
#' @examples
#' bw <- analyze_body_weight(master_synthetic_data, volume_column = "Volume",
#'                           volume_units = "mm3", cage_column = "Cage",
#'                           reference_group = "Vehicle")
#' names(bw)
#' @export
analyze_body_weight <- function(df,
                                weight_column    = "Weight",
                                time_column      = "Day",
                                treatment_column = "Treatment",
                                id_column        = "ID",
                                volume_column    = NULL,
                                sex_column       = NULL,
                                cage_column      = NULL,
                                adjust_tumor_weight = TRUE,
                                tumor_density    = 1.0,
                                volume_units     = NULL,
                                covariates       = c("volume"),
                                estimation       = c("REML", "ML"),
                                comparison_family = c("vs_reference", "all_pairs", "custom"),
                                custom_contrasts = NULL,
                                p_adjust_method  = c("bonferroni", "holm", "fdr",
                                                     "dunnett", "tukey", "none"),
                                reference_group  = NULL) {

  estimation <- match.arg(estimation)
  comparison_family <- match.arg(comparison_family)
  p_adjust_method   <- match.arg(p_adjust_method)

  # CODE_REVIEW.md R3.12 / G.1 — same comparison-family contract as
  # tumor_growth_statistics(); both paths here fit a joint model, so Tukey and
  # Dunnett are available.
  comparison_spec <- resolve_comparison_spec(
    comparison_family, p_adjust_method, custom_contrasts, supports_joint = TRUE
  )

  # --- Validate required columns ---
  required <- c(weight_column, time_column, treatment_column, id_column)
  missing_cols <- setdiff(required, names(df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # --- Prepare working data ---
  wd <- data.frame(
    ID        = as.factor(df[[id_column]]),
    Treatment = as.factor(df[[treatment_column]]),
    Day       = as.numeric(df[[time_column]]),
    Weight    = as.numeric(df[[weight_column]]),
    stringsAsFactors = FALSE
  )

  # Tumor volume (for adjustment and/or covariate). R20.22: volume is filled
  # in on weighing days without a calliper measurement, so those weighings are
  # kept rather than dropped.
  has_volume <- !is.null(volume_column) && volume_column %in% names(df)
  if (has_volume) {
    fill_key <- make_mouse_key(
      as.character(df[[treatment_column]]), as.character(df[[id_column]]),
      if (!is.null(cage_column) && cage_column %in% names(df))
        as.character(df[[cage_column]]) else "")
    wd$Volume <- me_fill_volume(fill_key, wd$Day, as.numeric(df[[volume_column]]))
  }

  # Tumor weight adjustment
  if (adjust_tumor_weight && has_volume) {
    # CODE_REVIEW.md R3.30 — units are resolved explicitly (auto-detected when
    # NULL) rather than assuming mm³, and the resulting mass is range-checked
    # against body weight so a 1000x unit error cannot pass silently.
    volume_units <- resolve_volume_units(wd$Volume, volume_units)
    wd$Tumor_Weight <- volume_to_mass(wd$Volume, tumor_density, volume_units)
    check_tumor_mass_plausible(wd$Tumor_Weight, wd$Weight, volume_units)
    wd$Net_Weight   <- wd$Weight - wd$Tumor_Weight
    response_col    <- "Net_Weight"
  } else {
    wd$Net_Weight   <- wd$Weight
    response_col    <- "Net_Weight"
  }

  # Sex covariate
  has_sex <- !is.null(sex_column) && sex_column %in% names(df)
  if (has_sex) {
    wd$Sex <- as.factor(df[[sex_column]])
  }

  # Cage
  has_cage <- !is.null(cage_column) && cage_column %in% names(df)
  if (has_cage) {
    wd$Cage <- as.factor(df[[cage_column]])
    warn_reused_cage_labels(df, cage_column, treatment_column)   # R22.3
  }

  # Initial mass per mouse (baseline weight at earliest timepoint).
  #
  # CODE_REVIEW.md R3.11 — aggregating by ID alone collapses mice that share a
  # numeric ear-tag ID across treatment groups or cages (the normal case), so
  # they all received one another's baseline. Key on the composite mouse key
  # like the rest of the package, and filter to each mouse's own earliest day
  # rather than relying on row order (the Round 1 1.7 defect class).
  wd$.MouseKey <- if (has_cage) {
    make_mouse_key(as.character(wd$Treatment), as.character(wd$ID),
                   as.character(wd$Cage))
  } else {
    make_mouse_key(as.character(wd$Treatment), as.character(wd$ID))
  }
  wd <- wd[order(wd$.MouseKey, wd$Day), ]
  # CODE_REVIEW.md T1 / R20.3 -- the random effects were (1 + Day | ID), so a
  # mouse "1" in every arm became one animal (null familywise rejection 0.71,
  # SE ratio 0.30, n_subjects 8 for 24 mice). Group by the composite key.
  wd$Animal <- factor(wd$.MouseKey)


  # Remove rows with missing response
  wd <- wd[!is.na(wd$Net_Weight) & !is.na(wd$Day), ]
  if (nrow(wd) == 0) stop("No valid weight observations after cleaning.")

  # R20.44: the baseline was computed before this filter, from each animal's
  # first row. Pre-implant weighings (no volume, so no net weight) made it fail
  # with "no rows to aggregate", and a missing first weight gave the animal a
  # missing baseline, so lmer silently dropped all of its rows. It is now each
  # animal's first available weight, and only when it is a covariate.
  if ("initial_mass" %in% covariates) {
    baseline <- me_per_mouse_baseline(wd, ".MouseKey", "Net_Weight", "Day",
                                      out_name = "Initial_Mass")
    wd <- merge(wd, baseline, by = ".MouseKey", all.x = TRUE)
    wd <- wd[order(wd$.MouseKey, wd$Day), ]
  }

  # R20.46: shared reference rule. A reference_group that matched no arm used
  # to be ignored silently, and the first level used instead.
  reference_group <- me_resolve_reference(wd$Treatment, reference_group,
                                          "analyze_body_weight")
  wd$Treatment <- stats::relevel(factor(as.character(wd$Treatment)),
                                 ref = reference_group)
  # The level the vs_reference contrasts are taken against.
  ref_level <- reference_group

  # --- Build formula ---
  # CODE_REVIEW.md R3.10 — `adjust_tumor_weight` and a "volume" covariate are
  # two mutually exclusive ways of handling the same confounder, and the
  # shipped defaults applied both: Net_Weight already has a deterministic
  # linear function of Volume subtracted from it, so entering Volume as a
  # predictor of that response makes its coefficient uninterpretable (it
  # absorbs the residual of a quantity already removed) and adjusts the
  # treatment effect for tumour burden twice, in two different functional
  # forms. Subtract the mass *or* condition on it, never both.
  if (adjust_tumor_weight && has_volume && "volume" %in% covariates) {
    covariates <- setdiff(covariates, "volume")
    message(
      "Tumour burden is already removed from the response via ",
      "adjust_tumor_weight = TRUE; dropping the redundant 'volume' covariate. ",
      "To condition on tumour volume instead of subtracting its mass, set ",
      "adjust_tumor_weight = FALSE."
    )
  }

  fixed_terms <- "Treatment * Day"
  if ("volume" %in% covariates && has_volume) {
    fixed_terms <- paste(fixed_terms, "+ Volume")
  }
  if ("sex" %in% covariates && has_sex) {
    fixed_terms <- paste(fixed_terms, "+ Sex")
  }
  if ("initial_mass" %in% covariates) {
    fixed_terms <- paste(fixed_terms, "+ Initial_Mass")
  }

  # --- Fit model with auto-fallback ---
  model <- NULL
  model_simplified <- FALSE

  # Random-effects spec: per-mouse intercept + slope, plus cage random
  # intercept when a cage column was supplied. CODE_REVIEW.md J.11 — the
  # cage column was previously attached to the data frame but never
  # appeared in the formula, the same silent-ignore bug class as
  # Round 1 1.1 (handle_cage_effects in tumor_growth_statistics).
  re_full   <- if (has_cage) "(1 + Day | Animal) + (1 | Cage)" else "(1 + Day | Animal)"
  re_simple <- if (has_cage) "(1 | Animal) + (1 | Cage)"        else "(1 | Animal)"

  # R20.73: tryCatch() abandoned the full model on ANY warning, and the summary
  # then said "simplified due to sample size" whatever the cause. Only a
  # failure, a singular fit or non-convergence now leads to the simpler model;
  # the reason is recorded, and other warnings are passed on.
  fit_lmer <- function(f) {
    w <- character()
    m <- withCallingHandlers(
      tryCatch(lme4::lmer(f, data = wd, REML = (estimation == "REML")),
               error = function(e) { w <<- c(w, conditionMessage(e)); NULL }),
      warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
    # R20-T: a failed gradient check is retried from the estimates before
    # the random slope is given up; under lme4 2.0.6 the shipped data's
    # model warned spuriously and lost its slopes.
    rs <- me_lmer_restart(m, w)
    list(model = rs$fit, warnings = rs$msgs, note = rs$note)
  }
  formula_full <- stats::as.formula(paste(response_col, "~", fixed_terms, "+", re_full))
  full <- fit_lmer(formula_full)
  simplify_reason <- if (is.null(full$model)) {
    "the random-slope model could not be fitted"
  } else if (lme4::isSingular(full$model)) {
    "the random-slope model was singular (slope variance estimated at zero)"
  } else if (any(grepl("converge", full$warnings, ignore.case = TRUE))) {
    "the random-slope model did not converge"
  } else NULL
  if (is.null(simplify_reason)) {
    model <- full$model
    for (w in full$warnings) warning(w, call. = FALSE)
  } else {
    formula_simple <- stats::as.formula(paste(response_col, "~", fixed_terms, "+", re_simple))
    simple <- fit_lmer(formula_simple)
    if (is.null(simple$model)) {
      stop("Mixed model failed to converge: ", paste(simple$warnings, collapse = "; "))
    }
    model <- simple$model
    model_simplified <- TRUE
  }

  # --- Extract results ---
  model_summary <- summary(model)
  fe <- as.data.frame(stats::coef(model_summary))
  fe$Term <- rownames(fe)
  fe <- fe[, c("Term", setdiff(names(fe), "Term"))]
  rownames(fe) <- NULL

  re <- as.data.frame(lme4::VarCorr(model))

  # R20.44 (decided 2026-10-09): the primary comparison is the rate of weight
  # change (g/day), each arm against the reference, from the Treatment x Day
  # term. It used to be the weight level at the mean study day, which carried
  # any baseline imbalance (+1.5 g at baseline showed as a +1 g "effect") and
  # had power 0.14 where the slope contrast has 0.84. The level comparison is
  # kept as level_means / level_comparisons.
  trend_obj <- tryCatch(emmeans::emtrends(model, ~ Treatment, var = "Day"),
                        error = function(e) NULL)
  emm <- if (!is.null(trend_obj)) {
    t <- as.data.frame(trend_obj)
    names(t)[names(t) == "Day.trend"] <- "Weight_Change_Per_Day"
    t
  } else NULL
  pairwise <- if (!is.null(trend_obj)) {
    bw_pairwise_table(trend_obj, comparison_spec, ref_level, custom_contrasts)
  } else NULL

  # Weight level at the mean study day (R3.12), secondary.
  emm_obj <- tryCatch(
    emmeans::emmeans(model, ~ Treatment, at = list(Day = mean(wd$Day))),
    error = function(e) NULL
  )
  level_means <- if (!is.null(emm_obj)) as.data.frame(emm_obj) else NULL
  level_pairwise <- if (!is.null(emm_obj)) {
    bw_pairwise_table(emm_obj, comparison_spec, ref_level, custom_contrasts)
  } else NULL

  n_obs_model  <- stats::nobs(model)
  n_subj_model <- as.integer(lme4::ngrps(model)[["Animal"]])

  # --- Summary text ---
  lines <- c(
    "=== BODY WEIGHT MIXED MODEL ===",
    "",
    sprintf("Response: %s", if (adjust_tumor_weight && has_volume) "Net Weight (body - tumor)" else "Body Weight"),
    sprintf("Estimation: %s", estimation),
    sprintf("Random effects: %s", if (model_simplified) "(1 | Animal) [simplified]" else "(1 + Day | Animal)"),
    sprintf("Fixed effects: %s", fixed_terms),
    sprintf("Observations: %d  |  Subjects: %d  |  Groups: %d",
            n_obs_model, n_subj_model, length(levels(wd$Treatment))),
    "Comparison: rate of weight change (g/day) vs the reference, from Treatment x Day",
    ""
  )

  if (model_simplified) {
    lines <- c(lines,
      paste0("NOTE: random intercept only, because ", simplify_reason, "."),
      "")
  }

  summary_text <- paste(lines, collapse = "\n")

  # ── Diagnostic plots (LMM path only) ───────────────────────────────────────
  # CODE_REVIEW.md J.12 — the Bayesian counterpart returns pp_check_plot,
  # mcmc_trace_plot, etc. The frequentist version returned just model +
  # summary_text; users had no way to check residual normality or
  # heteroscedasticity. Add QQ and residuals-vs-fitted plots.
  # CODE_REVIEW.md J.12 — Bayesian counterpart returns pp_check_plot,
  # mcmc_trace_plot, etc. v0.4.8 — use the shared helper so TG/BW/AUC
  # paths produce identical diagnostic shape.
  rd <- build_residual_diagnostic_plots(model, title_prefix = "Body-weight LMM")
  diag_qq_plot             <- rd$diag_qq_plot
  diag_resid_fitted_plot   <- rd$diag_resid_fitted_plot
  diag_scale_location_plot <- rd$diag_scale_location_plot
  diag_re_qq_plot          <- build_random_effects_qq_plot(
    model, title_prefix = "Body-weight LMM")
  # CODE_REVIEW.md DIAGNOSTICS gap (13) — LMM influence diagnostics.
  lmm_infl <- build_lmm_influence(model, groups = "Animal")

  list(
    reference_group = reference_group,   # R20.46
    model          = model,
    fixed_effects  = fe,
    random_effects = re,
    estimand       = "difference in rate of weight change (g/day) vs the reference",
    emmeans_table  = emm,               # weight change per day by arm (R20.44)
    pairwise_comparisons = pairwise,    # slope differences vs the reference
    level_means       = level_means,    # weight at the mean study day
    level_comparisons = level_pairwise,
    comparison_family    = comparison_spec$family,
    p_adjust_method_used = comparison_spec$p_adjust_method,
    model_info     = list(
      estimation       = estimation,
      response         = response_col,
      fixed_formula    = fixed_terms,
      model_simplified = model_simplified,
      simplify_reason  = if (model_simplified) simplify_reason else NA_character_,
      convergence_note = full$note,   # R20-T: how a gradient warning was checked
      adjust_tumor     = adjust_tumor_weight && has_volume,
      tumor_density    = tumor_density,
      cage_in_model    = has_cage,
      n_obs            = n_obs_model,     # R20.44: what the model used
      n_subjects       = n_subj_model,
      n_groups         = length(levels(wd$Treatment))
    ),
    diag_qq_plot             = diag_qq_plot,
    diag_resid_fitted_plot   = diag_resid_fitted_plot,
    diag_scale_location_plot = diag_scale_location_plot,
    diag_re_qq_plot          = diag_re_qq_plot,
    diag_cooks_distance      = if (!is.null(lmm_infl)) lmm_infl$cooks_distance,
    diag_dfbetas             = if (!is.null(lmm_infl)) lmm_infl$dfbetas,
    weight_data    = wd,
    summary_text   = summary_text
  )
}
