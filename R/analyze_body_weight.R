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
#' @param tumor_density Density in g/cm³ for tumor weight estimation (default 1.0).
#' @param covariates Character vector of optional covariates: "volume", "sex", "initial_mass".
#' @param estimation Character; "REML" (default) or "ML".
#' @param model_type Character; "lmm" (default — linear mixed model via lme4) or
#'   "gam" (generalized additive mixed model via gamm4 with a group-specific
#'   smoother on Day, preferred when weight trajectories are non-monotonic).
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
                                model_type       = c("lmm", "gam"),
                                comparison_family = c("vs_reference", "all_pairs", "custom"),
                                custom_contrasts = NULL,
                                p_adjust_method  = c("bonferroni", "holm", "fdr",
                                                     "dunnett", "tukey", "none"),
                                reference_group  = NULL) {

  estimation <- match.arg(estimation)
  model_type <- match.arg(model_type)
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

  # Tumor volume (for adjustment and/or covariate)
  has_volume <- !is.null(volume_column) && volume_column %in% names(df)
  if (has_volume) {
    wd$Volume <- as.numeric(df[[volume_column]])
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

  first_day <- stats::aggregate(Day ~ .MouseKey, data = wd, FUN = min)
  names(first_day)[2] <- ".FirstDay"
  bl_rows <- merge(wd[, c(".MouseKey", "Day", "Net_Weight")], first_day,
                   by = ".MouseKey")
  bl_rows <- bl_rows[bl_rows$Day == bl_rows$.FirstDay, , drop = FALSE]
  baseline <- stats::aggregate(Net_Weight ~ .MouseKey, data = bl_rows,
                               FUN = mean)
  names(baseline)[2] <- "Initial_Mass"
  wd <- merge(wd, baseline, by = ".MouseKey", all.x = TRUE)

  # Remove rows with missing response
  wd <- wd[!is.na(wd$Net_Weight) & !is.na(wd$Day), ]
  if (nrow(wd) == 0) stop("No valid weight observations after cleaning.")

  # Set reference group
  if (!is.null(reference_group) && reference_group %in% levels(wd$Treatment)) {
    wd$Treatment <- stats::relevel(wd$Treatment, ref = reference_group)
  }
  # The level the vs_reference contrasts are taken against.
  ref_level <- levels(wd$Treatment)[1]

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

  # GAM path — fit via gamm4 with a group-specific smoother on Day and
  # return early with a compatible result-list shape.
  if (model_type == "gam") {
    if (!requireNamespace("gamm4", quietly = TRUE)) {
      stop(
        "Package 'gamm4' is required for model_type = 'gam'.\n",
        "Install it with: install.packages('gamm4')"
      )
    }
    n_days <- length(unique(wd$Day))
    k_val  <- max(3L, min(10L, n_days - 1L))

    gam_fixed <- paste0(
      response_col, " ~ Treatment + s(Day, by = Treatment, k = ", k_val, ")"
    )
    if ("volume" %in% covariates && has_volume) {
      gam_fixed <- paste(gam_fixed, "+ Volume")
    }
    if ("sex" %in% covariates && has_sex) {
      gam_fixed <- paste(gam_fixed, "+ Sex")
    }
    if ("initial_mass" %in% covariates) {
      gam_fixed <- paste(gam_fixed, "+ Initial_Mass")
    }
    gam_random <- if (has_cage) {
      stats::as.formula("~ (1 | Cage) + (1 | ID)")
    } else {
      stats::as.formula("~ (1 | ID)")
    }

    gam_fit <- tryCatch(
      gamm4::gamm4(
        formula = stats::as.formula(gam_fixed),
        random  = gam_random,
        data    = wd,
        REML    = (estimation == "REML")
      ),
      error = function(e) {
        stop("gamm4 fit failed: ", conditionMessage(e))
      }
    )

    # CODE_REVIEW.md R3.4 — gamm4's $gam component is a stub whose class vector
    # lacks c("glm","lm") and whose $call is NULL, so emmeans::recover_data.gam
    # rejects it with "Can't handle an object of class 'NULL'". v0.4.11 patched
    # this for the tumour-growth path inside tgs_fit_gamm4_model(), but this
    # function fits gamm4 inline and never got the patch — so the emmeans call
    # below always errored, the tryCatch turned it into NULL, and the
    # body-weight GAM path silently returned an empty marginal-means table,
    # indistinguishable from "no effect".
    gam_fit <- patch_gamm4_stub(gam_fit)

    emm_obj <- tryCatch(
      emmeans::emmeans(gam_fit$gam, ~ Treatment,
                       at = list(Day = mean(wd$Day))),
      error = function(e) {
        warning("emmeans on the fitted GAMM failed: ", conditionMessage(e),
                call. = FALSE)
        NULL
      }
    )
    emm <- if (!is.null(emm_obj)) as.data.frame(emm_obj) else NULL

    # CODE_REVIEW.md R3.12 — the function previously returned group marginal
    # means and nothing inferential, so it could not answer "did this arm lose
    # more weight than control", the primary toxicity question.
    pairwise <- if (!is.null(emm_obj)) {
      bw_pairwise_table(emm_obj, comparison_spec, ref_level, custom_contrasts)
    } else NULL

    return(list(
      model          = gam_fit,
      fixed_effects  = data.frame(
        Term = "Smooth: s(Day, by = Treatment)",
        Note = paste0("k = ", k_val, " (auto-chosen)"),
        stringsAsFactors = FALSE
      ),
      random_effects = tryCatch(
        as.data.frame(lme4::VarCorr(gam_fit$mer)),
        error = function(e) NULL
      ),
      emmeans_table  = emm,
      pairwise_comparisons = pairwise,
      comparison_family    = comparison_spec$family,
      p_adjust_method_used = comparison_spec$p_adjust_method,
      model_info     = list(
        estimation       = estimation,
        response         = response_col,
        fixed_formula    = gam_fixed,
        model_type       = "gam",
        model_simplified = FALSE,
        adjust_tumor     = adjust_tumor_weight && has_volume,
        tumor_density    = tumor_density,
        smoother_k       = k_val,
        n_obs            = nrow(wd),
        n_subjects       = length(unique(wd$ID)),
        n_groups         = length(levels(wd$Treatment))
      ),
      weight_data    = wd,
      summary_text   = paste(c(
        "=== BODY WEIGHT GAMM ===",
        "",
        sprintf("Response: %s",
                if (adjust_tumor_weight && has_volume) "Net Weight (body - tumor)" else "Body Weight"),
        sprintf("Estimation: %s", estimation),
        sprintf("Random effects: %s", deparse1(gam_random)),
        sprintf("Fixed effects: %s (smoother basis k = %d)", gam_fixed, k_val),
        sprintf("Observations: %d  |  Subjects: %d  |  Groups: %d",
                nrow(wd), length(unique(wd$ID)), length(levels(wd$Treatment)))
      ), collapse = "\n")
    ))
  }

  # Random-effects spec: per-mouse intercept + slope, plus cage random
  # intercept when a cage column was supplied. CODE_REVIEW.md J.11 — the
  # cage column was previously attached to the data frame but never
  # appeared in the formula, the same silent-ignore bug class as
  # Round 1 1.1 (handle_cage_effects in tumor_growth_statistics).
  re_full   <- if (has_cage) "(1 + Day | ID) + (1 | Cage)" else "(1 + Day | ID)"
  re_simple <- if (has_cage) "(1 | ID) + (1 | Cage)"        else "(1 | ID)"

  # Try random slope + intercept first
  formula_full <- stats::as.formula(
    paste(response_col, "~", fixed_terms, "+", re_full)
  )
  tryCatch({
    model <- lme4::lmer(formula_full, data = wd, REML = (estimation == "REML"))
    # Check for singular fit
    if (lme4::isSingular(model)) {
      model <- NULL
    }
  }, error = function(e) {
    model <<- NULL
  }, warning = function(w) {
    if (grepl("singular|converge", w$message, ignore.case = TRUE)) {
      model <<- NULL
    }
  })

  # Fallback to random intercept only (still preserves cage RE when present)
  if (is.null(model)) {
    formula_simple <- stats::as.formula(
      paste(response_col, "~", fixed_terms, "+", re_simple)
    )
    tryCatch({
      model <- lme4::lmer(formula_simple, data = wd, REML = (estimation == "REML"))
      model_simplified <- TRUE
    }, error = function(e) {
      stop("Mixed model failed to converge: ", conditionMessage(e))
    })
  }

  # --- Extract results ---
  model_summary <- summary(model)
  fe <- as.data.frame(stats::coef(model_summary))
  fe$Term <- rownames(fe)
  fe <- fe[, c("Term", setdiff(names(fe), "Term"))]
  rownames(fe) <- NULL

  re <- as.data.frame(lme4::VarCorr(model))

  # Marginal means per treatment.
  # CODE_REVIEW.md R3.12 — marginalise explicitly at the mean study day so this
  # matches how tumor_growth_statistics() defines a group's adjusted mean; the
  # bare `~ Treatment` relied on emmeans' default reference grid, so the two
  # functions documented the same quantity two different ways.
  emm_obj <- tryCatch(
    emmeans::emmeans(model, ~ Treatment, at = list(Day = mean(wd$Day))),
    error = function(e) NULL
  )
  emm <- if (!is.null(emm_obj)) as.data.frame(emm_obj) else NULL

  pairwise <- if (!is.null(emm_obj)) {
    bw_pairwise_table(emm_obj, comparison_spec, ref_level, custom_contrasts)
  } else NULL

  # --- Summary text ---
  lines <- c(
    "=== BODY WEIGHT MIXED MODEL ===",
    "",
    sprintf("Response: %s", if (adjust_tumor_weight && has_volume) "Net Weight (body - tumor)" else "Body Weight"),
    sprintf("Estimation: %s", estimation),
    sprintf("Random effects: %s", if (model_simplified) "(1 | ID) [simplified]" else "(1 + Day | ID)"),
    sprintf("Fixed effects: %s", fixed_terms),
    sprintf("Observations: %d  |  Subjects: %d  |  Groups: %d",
            nrow(wd), length(unique(wd$ID)), length(levels(wd$Treatment))),
    ""
  )

  if (model_simplified) {
    lines <- c(lines,
      "NOTE: Model simplified due to sample size. Random slope removed; using random intercept only.",
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
  rd <- if (model_type == "lmm")
    build_residual_diagnostic_plots(model, title_prefix = "Body-weight LMM")
  else
    list(diag_qq_plot = NULL,
         diag_resid_fitted_plot = NULL,
         diag_scale_location_plot = NULL)
  diag_qq_plot             <- rd$diag_qq_plot
  diag_resid_fitted_plot   <- rd$diag_resid_fitted_plot
  diag_scale_location_plot <- rd$diag_scale_location_plot
  diag_re_qq_plot          <- if (model_type == "lmm")
    build_random_effects_qq_plot(model, title_prefix = "Body-weight LMM")
  else NULL
  # CODE_REVIEW.md DIAGNOSTICS gap (13) — LMM influence diagnostics.
  lmm_infl <- if (model_type == "lmm") build_lmm_influence(model) else NULL

  list(
    model          = model,
    fixed_effects  = fe,
    random_effects = re,
    emmeans_table  = emm,
    pairwise_comparisons = pairwise,
    comparison_family    = comparison_spec$family,
    p_adjust_method_used = comparison_spec$p_adjust_method,
    model_info     = list(
      estimation       = estimation,
      response         = response_col,
      fixed_formula    = fixed_terms,
      model_simplified = model_simplified,
      adjust_tumor     = adjust_tumor_weight && has_volume,
      tumor_density    = tumor_density,
      cage_in_model    = has_cage,
      n_obs            = nrow(wd),
      n_subjects       = length(unique(wd$ID)),
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
