# Copyright (c) 2026 mouseExperiment Contributors
# Licensed under the MIT License - see LICENSE file

#' Extrapolate missing data points for subjects
#'
#' @param df Data frame with tumor growth data
#' @param id_column Column name for subject IDs
#' @param treatment_column Column name for treatment groups
#' @param cage_column Column name for cage identifiers
#' @param time_column Column name for time points
#' @param volume_column Column name for tumor volume
#' @param verbose Print progress messages
#' @return Modified data frame with extrapolated points
#' @noRd
#' @keywords internal
tgs_extrapolate <- function(df, id_column, treatment_column, cage_column,
                            time_column, volume_column, verbose) {
  if (isTRUE(verbose)) message("Extrapolating points for subjects with missing data at the last timepoint")
  
  # Find the true maximum day across all subjects (global maximum day of the study)
  true_max_day <- max(df[[time_column]], na.rm = TRUE)
  if (isTRUE(verbose)) message("True maximum day of the study: ", true_max_day)
  
  # Make sure all data has the Extrapolated column
  if (!"Extrapolated" %in% colnames(df)) {
    df$Extrapolated <- FALSE
  }
  
  # Create a list to store results for each subject
  subjects_with_extrapolation <- list()
  
  # Process each unique subject
  unique_subjects <- unique(
    make_mouse_key(df[[id_column]], df[[treatment_column]], df[[cage_column]])
  )

  for (subject_id in unique_subjects) {
    # Parse the composite ID
    id_parts <- split_mouse_key(subject_id)
    id <- id_parts[1]
    treatment <- id_parts[2]
    cage <- id_parts[3]
    
    # Get data for this subject
    subject_data <- df[df[[id_column]] == id & 
                      df[[treatment_column]] == treatment & 
                      df[[cage_column]] == cage, ]
    
    # Get the max day for this subject
    max_subj_day <- max(subject_data[[time_column]])
    
    # Only extrapolate if the subject doesn't have data on the true max day
    if (max_subj_day < true_max_day) {
      # Need at least 2 points for extrapolation
      if (nrow(subject_data) >= 2) {
        # Use the last 3 points (or all if less than 3) to fit a linear model
        n_points <- min(3, nrow(subject_data))
        subject_data <- subject_data[order(subject_data[[time_column]]), ]
        last_points <- tail(subject_data, n_points)
        
        # Try to fit model
        tryCatch({
          lm_fit <- stats::lm(paste(volume_column, "~", time_column), data = last_points)
          
          # Predict at true_max_day
          new_data <- data.frame(time = true_max_day)
          names(new_data) <- time_column
          predicted_volume <- max(0, as.numeric(predict(lm_fit, newdata = new_data)))
          
          # Create a new row for the extrapolated point
          new_row <- subject_data[1, ]
          new_row[[time_column]] <- true_max_day
          new_row[[volume_column]] <- predicted_volume
          new_row$Extrapolated <- TRUE
          
          # Add the new extrapolated point to the subject data
          subject_data <- rbind(subject_data, new_row)
          
          if (isTRUE(verbose)) {
            message("Extrapolated subject ", id, " from day ", max_subj_day, " to day ", true_max_day)
          }
        }, error = function(e) {
          if (isTRUE(verbose)) {
            message("Failed to extrapolate subject ", id, ": ", conditionMessage(e))
          }
        })
      }
    }
    
    # Store the processed subject data
    subjects_with_extrapolation[[subject_id]] <- subject_data
  }
  
  # Combine all subject data
  df <- do.call(rbind, subjects_with_extrapolation)
  
  # Count extrapolated subjects for verbose output.
  # Previously: unique(df$Extrapolated[df$Extrapolated]) collapses to c(TRUE),
  # so the count was always 0 or 1. Use composite mouse keys for the actual
  # subject count.
  if (isTRUE(verbose)) {
    ex_rows <- df[df$Extrapolated, , drop = FALSE]
    n_extrapolated <- if (nrow(ex_rows) > 0L) {
      length(unique(make_mouse_key(
        ex_rows[[id_column]],
        ex_rows[[treatment_column]],
        ex_rows[[cage_column]]
      )))
    } else 0L
    if (n_extrapolated > 0) {
      message("Successfully extrapolated ", n_extrapolated, " subjects to day ", true_max_day)
    } else {
      message("No subjects needed or qualified for extrapolation to day ", true_max_day)
    }
  }
  df
}

#' Compute per-subject growth rates
#'
#' @param auc_df Untransformed data frame
#' @param treatment_column Column name for treatment groups
#' @param id_column Column name for subject IDs
#' @param cage_column Column name for cage identifiers
#' @param time_column Column name for time points
#' @param volume_column Column name for tumor volume
#' @return Data frame of growth rates or NULL on error
#' @noRd
#' @keywords internal
tgs_compute_growth_rates <- function(auc_df, treatment_column, id_column,
                                     cage_column, time_column, volume_column) {
  tryCatch({
    # Calculate growth rates for each subject
    # Use auc_df (untransformed volumes) so growth rates are not double-transformed
    # CODE_REVIEW.md R3.24 — split() on a list of factors allocates one element
    # per *combination of levels*, not per observed combination: 6 treatments x
    # 60 IDs x 12 cages is 4,320 elements of which ~60 are non-empty. Split on
    # the composite key instead.
    split_auc <- split(
      auc_df,
      make_mouse_key(auc_df[[treatment_column]],
                     auc_df[[id_column]],
                     auc_df[[cage_column]]),
      drop = TRUE
    )
    growth_rates_list <- lapply(split_auc, function(subject_data) {
      if (nrow(subject_data) >= 3) {
        # Sort by time
        subject_data <- subject_data[order(subject_data[[time_column]]), ]
        
        # Calculate log volume from original (untransformed) data.
        # CODE_REVIEW.md R3.20 — same all-non-positive guard as the main path;
        # a subject with no positive volumes is skipped rather than producing
        # log(Inf).
        raw_vol       <- subject_data[[volume_column]]
        positive_vals <- raw_vol[is.finite(raw_vol) & raw_vol > 0]
        if (length(positive_vals) == 0L) return(NULL)
        raw_vol[raw_vol <= 0] <- min(positive_vals, na.rm = TRUE) / 2
        log_volume <- log(raw_vol)
        
        # Fit linear model
        model <- stats::lm(log_volume ~ subject_data[[time_column]])
        growth_rate <- stats::coef(model)[2]
        ci <- tryCatch(stats::confint(model, level = 0.95)[2, ], error = function(e) c(NA_real_, NA_real_))

        # Return data frame with results including cage information
        data.frame(
          Treatment          = unique(subject_data[[treatment_column]]),
          ID                 = unique(subject_data[[id_column]]),
          Cage               = unique(subject_data[[cage_column]]),
          growth_rate        = growth_rate,
          growth_rate_lower  = ci[1],
          growth_rate_upper  = ci[2],
          R_squared          = round(summary(model)$r.squared, 4)
        )
      } else {
        NULL
      }
    })
    
    # Combine results
    do.call(rbind, growth_rates_list)
  }, error = function(e) {
    warning("Error calculating growth rates: ", e$message)
    NULL
  })
}

#' Analyse cage effects and collinearity with treatment
#'
#' @param analysis_df Transformed analysis data frame
#' @param cage_column Column name for cage identifiers
#' @param treatment_column Column name for treatment groups
#' @param volume_column Column name for tumor volume
#' @return List with collinearity_test and effects components
#' @noRd
#' @keywords internal
tgs_compute_cage_effects <- function(analysis_df, cage_column, treatment_column,
                                     volume_column) {
  # Cage effect analysis
  cage_analysis <- list()
  
  # Test for collinearity between cage and treatment
  cage_treatment_table <- table(analysis_df[[cage_column]], analysis_df[[treatment_column]])
  cage_analysis$collinearity_test <- tryCatch({
    stats::chisq.test(cage_treatment_table)
  }, error = function(e) {
    warning("Error in chi-square test: ", e$message)
    NULL
  })
  
  # Calculate cage-level effects
  cage_effects <- tryCatch({
    # Split data by cage and treatment
    split_data <- split(analysis_df, list(analysis_df[[cage_column]], analysis_df[[treatment_column]]))
    
    # Calculate statistics for each group
    cage_effects_list <- lapply(split_data, function(group_data) {
      if (nrow(group_data) > 0) {
        # Only create entries for non-empty groups
        data.frame(
          Cage = unique(group_data[[cage_column]]),
          Treatment = unique(group_data[[treatment_column]]),
          mean_volume = mean(group_data[[volume_column]], na.rm = TRUE),
          sd_volume = stats::sd(group_data[[volume_column]], na.rm = TRUE),
          n = nrow(group_data)
        )
      } else {
        # Skip empty groups
        NULL
      }
    })
    
    # Filter out NULL results and combine
    cage_effects_list <- cage_effects_list[!sapply(cage_effects_list, is.null)]
    if (length(cage_effects_list) > 0) {
      do.call(rbind, cage_effects_list)
    } else {
      NULL
    }
  }, error = function(e) {
    warning("Error calculating cage effects: ", e$message)
    NULL
  })
  
  cage_analysis$effects <- cage_effects
  cage_analysis
}

#' Fit lme4 mixed-effects model using caller-specified structure
#'
#' @param analysis_df Transformed analysis data frame
#' @param volume_column Column name for tumor volume
#' @param time_column Column name for time points
#' @param treatment_column Column name for treatment groups
#' @param id_column Column name for subject IDs
#' @param cage_column Column name for cage identifiers
#' @param random_effects_specification One of "intercept_only", "slope", "none"
#' @param cage_handling Output of \code{resolve_cage_handling()}: a list with
#'   \code{fixed} / \code{random} logicals and a \code{reason} string.
#' @param verbose Print progress messages
#' @return List with model, model_selection, and best_model
#' @noRd
#' @keywords internal
tgs_fit_lme4_models <- function(analysis_df, volume_column, time_column,
                                treatment_column, id_column, cage_column,
                                random_effects_specification,
                                cage_handling,
                                verbose,
                                include_necrotic_covariate = FALSE) {
  bt <- function(x) paste0("`", x, "`")

  # Linear time term only — higher-order polynomial time has been removed in
  # favour of model_type = "gam" (smoother-based non-linear trajectories).
  time_term <- bt(time_column)

  # CODE_REVIEW.md R3.17 — cage placement is decided by resolve_cage_handling()
  # from the design structure, not by a chi-square p-value here.
  include_cage_fixed  <- isTRUE(cage_handling$fixed)
  include_cage_random <- isTRUE(cage_handling$random)

  # Build fixed-effects string
  fixed_part <- paste(
    bt(volume_column), "~", time_term, "*", bt(treatment_column)
  )
  if (include_cage_fixed) {
    fixed_part <- paste(fixed_part, "+", bt(cage_column))
  }
  if (isTRUE(include_necrotic_covariate)) {
    fixed_part <- paste(fixed_part, "+ necrotic_cov_flag")
  }

  if (isTRUE(verbose)) {
    message("Random effects specification: ", random_effects_specification)
    message("Cage in fixed effects: ", include_cage_fixed)
    message("Cage as random effect: ", include_cage_random)
  }

  # No random effects — use lm()
  if (random_effects_specification == "none") {
    if (include_cage_random) {
      warning(
        "A cage random intercept was selected (", cage_handling$reason,
        ") but random_effects_specification = 'none' removes all random ",
        "effects; the cage term is dropped. Standard errors will not account ",
        "for cage-level clustering.", call. = FALSE
      )
    }
    model <- tryCatch(
      stats::lm(stats::as.formula(fixed_part), data = analysis_df),
      error = function(e) {
        warning("lm() failed (", e$message, "). Falling back to intercept-only lmer.")
        lme4::lmer(
          stats::as.formula(
            paste(fixed_part, "+ (1|", bt(id_column), ")")
          ),
          data = analysis_df,
          control = lme4::lmerControl(
            check.nobs.vs.nlev = "ignore",
            check.nobs.vs.nRE  = "ignore"
          )
        )
      }
    )
    model_selection <- list(
      aic = stats::AIC(model),
      bic = stats::BIC(model),
      selected_model = "none"
    )
    return(list(model = model, model_selection = model_selection,
                best_model = "none"))
  }

  # Build random-effects string
  # Use withCallingHandlers so boundary (singular) fit warnings don't discard
  # valid models — the model is still usable when a random-effect variance is
  # estimated near zero.
  re_part <- switch(random_effects_specification,
    intercept_only = paste0("(1|", bt(id_column), ")"),
    slope          = paste0("(", bt(time_column), "|", bt(id_column), ")")
  )
  if (include_cage_random) {
    re_part <- paste(re_part, "+", paste0("(1|", bt(cage_column), ")"))
  }

  formula_str <- paste(fixed_part, "+", re_part)

  model <- tryCatch({
    withCallingHandlers(
      lme4::lmer(
        stats::as.formula(formula_str),
        data = analysis_df,
        control = lme4::lmerControl(
          check.nobs.vs.nlev = "ignore",
          check.nobs.vs.nRE  = "ignore"
        )
      ),
      warning = function(w) {
        if (grepl("boundary|singular", w$message))
          message("Note: model has singular fit (random-effect variance near zero)")
        invokeRestart("muffleWarning")
      }
    )
  }, error = function(e) {
    warning("lmer() failed (", e$message, "). Falling back to intercept-only model.")
    fallback_fixed <- paste(
      bt(volume_column), "~", bt(time_column), "*", bt(treatment_column)
    )
    if (isTRUE(include_necrotic_covariate)) {
      fallback_fixed <- paste(fallback_fixed, "+ necrotic_cov_flag")
    }
    lme4::lmer(
      stats::as.formula(paste(fallback_fixed, "+ (1|", bt(id_column), ")")),
      data = analysis_df,
      control = lme4::lmerControl(
        check.nobs.vs.nlev = "ignore",
        check.nobs.vs.nRE  = "ignore"
      )
    )
  })

  model_selection <- list(
    aic = stats::AIC(model),
    bic = stats::BIC(model),
    selected_model = random_effects_specification
  )

  list(model = model, model_selection = model_selection,
       best_model = random_effects_specification)
}

#' Compute descriptive summary statistics by treatment and time
#'
#' @param analysis_df Transformed analysis data frame
#' @param treatment_column Column name for treatment groups
#' @param time_column Column name for time points
#' @param volume_column Column name for tumor volume
#' @return Data frame of summary statistics or NULL on error
#' @noRd
#' @keywords internal
tgs_compute_summary <- function(analysis_df, treatment_column, time_column,
                                volume_column) {
  # Create a basic summary of the data
  tryCatch({
    # Split data by treatment and time
    split_data <- split(analysis_df, list(analysis_df[[treatment_column]], analysis_df[[time_column]]))
    
    # Calculate summary statistics for each group
    summary_list <- lapply(split_data, function(group_data) {
      if (nrow(group_data) > 0) {
        # Only create entries for non-empty groups
        data.frame(
          Treatment = unique(group_data[[treatment_column]]),
          Day = unique(group_data[[time_column]]),
          mean_volume = mean(group_data[[volume_column]], na.rm = TRUE),
          sd_volume = stats::sd(group_data[[volume_column]], na.rm = TRUE),
          n = nrow(group_data)
        )
      } else {
        # Skip empty groups
        NULL
      }
    })
    
    # Filter out NULL results and combine
    summary_list <- summary_list[!sapply(summary_list, is.null)]
    if (length(summary_list) > 0) {
      # Combine results and calculate standard error
      summary_df <- do.call(rbind, summary_list)
      summary_df$sem_volume <- summary_df$sd_volume / sqrt(summary_df$n)
      summary_df
    } else {
      NULL
    }
  }, error = function(e) {
    warning("Error calculating data summary: ", e$message)
    NULL
  })
}

#' Calculate area under the curve for each subject
#'
#' @param auc_df Untransformed data frame
#' @param id_column Column name for subject IDs
#' @param treatment_column Column name for treatment groups
#' @param cage_column Column name for cage identifiers
#' @param time_column Column name for time points
#' @param volume_column Column name for tumor volume
#' @param verbose Print progress messages
#' @return List with individual and summary AUC results
#' @noRd
#' @keywords internal
#' Percentile bootstrap CI for the mean difference of two AUC vectors
#'
#' Per-mouse AUC is a derived statistic whose distribution can be skewed
#' (especially when extrapolation is in play), so the Welch's t-CI used
#' downstream is only approximate at small N. The bootstrap CI is
#' non-parametric and well-suited to this setting.
#'
#' Resamples each group with replacement \code{n_boot} times, computes the
#' mean difference per resample, and returns the empirical 2.5 / 97.5 %
#' percentiles. Returns \code{c(NA, NA)} when either group has fewer than
#' two observations or n_boot < 2.
#'
#' @noRd
tgs_boot_diff_ci <- function(x, y, n_boot = 999L, seed = NULL) {
  if (length(x) < 2L || length(y) < 2L || n_boot < 2L) {
    return(c(lower = NA_real_, upper = NA_real_))
  }
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
      get(".Random.seed", envir = .GlobalEnv)
    } else NULL
    on.exit({
      if (!is.null(old_seed)) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    })
    set.seed(seed)
  }
  diffs <- vapply(seq_len(n_boot), function(b) {
    mean(sample(x, length(x), replace = TRUE)) -
      mean(sample(y, length(y), replace = TRUE))
  }, numeric(1L))
  q <- stats::quantile(diffs, c(0.025, 0.975), names = FALSE, na.rm = TRUE)
  c(lower = q[1L], upper = q[2L])
}


tgs_compute_auc <- function(auc_df, id_column, treatment_column, cage_column,
                            time_column, volume_column, verbose) {
  # Calculate AUC for each subject
  if (isTRUE(verbose)) message("Calculating AUC for each subject")
  
  # Uses exported calculate_auc() utility from utils_auc.R
  
  # Distinguish mice that share an ID across cages / treatments via the
  # composite key.
  # CODE_REVIEW.md R3.23 — a sequential `unique_id` column was previously built
  # and merged on, then never used (the loop keys on `composite_id`). The merge
  # also silently reordered rows and would have duplicated them had the
  # combination table not been unique. Removed.
  auc_df_with_id <- auc_df
  composite_id <- make_mouse_key(
    auc_df_with_id[[id_column]],
    auc_df_with_id[[treatment_column]],
    auc_df_with_id[[cage_column]]
  )
  auc_rows <- vector("list", length(unique(composite_id)))
  auc_row_idx <- 0L
  
  # Get max experiment time to determine if extrapolation is needed
  max_experiment_time <- max(auc_df[[time_column]])
  
  for (unique_id in unique(composite_id)) {
    # Extract data for this unique ID
    id_parts <- split_mouse_key(unique_id)
    actual_id <- id_parts[1]
    treatment <- id_parts[2]
    cage <- id_parts[3]
    
    subject_data <- auc_df_with_id[composite_id == unique_id, ]
    subject_data <- subject_data[order(subject_data[[time_column]]), ]
    
    # Calculate AUC using trapezoidal method
    auc_value <- calculate_auc(subject_data[[time_column]], subject_data[[volume_column]])
    
    # Check if this subject's data contains any extrapolated points
    has_extrapolated <- FALSE
    if ("Extrapolated" %in% colnames(subject_data)) {
      has_extrapolated <- any(subject_data$Extrapolated, na.rm = TRUE)
    }
    
    # Get the true last observation time (excluding extrapolated points)
    if (has_extrapolated && any(subject_data$Extrapolated, na.rm = TRUE)) {
      # If there are extrapolated points, get the max day from non-extrapolated points
      non_extrapolated_data <- subject_data[!subject_data$Extrapolated, ]
      true_last_observation <- max(non_extrapolated_data[[time_column]])
    } else {
      # If no extrapolated points, just use the max day
      true_last_observation <- max(subject_data[[time_column]])
    }
    
    # Make sure we have NumPoints data
    if (!"Extrapolated" %in% colnames(subject_data)) {
      n_points <- nrow(subject_data)
    } else {
      n_points <- nrow(subject_data[!subject_data$Extrapolated, ])
    }
    
    # Add to results
    auc_row_idx <- auc_row_idx + 1L
    auc_rows[[auc_row_idx]] <- data.frame(
      ID = actual_id,
      Treatment = treatment,
      Cage = cage,
      Group = treatment, # Added Group column for compatibility with plot_auc
      AUC = auc_value,
      Last_Day = true_last_observation,
      First_Day = min(subject_data[[time_column]]),
      Extrapolated = has_extrapolated,
      NumPoints = n_points # Count only non-extrapolated points
    )
  }
  
  # Bind all rows at once
  auc_data <- do.call(rbind, auc_rows[seq_len(auc_row_idx)])
  
  # Calculate summary statistics
  auc_summary <- stats::aggregate(AUC ~ Treatment, data = auc_data, 
                                FUN = function(x) c(Mean = mean(x), 
                                                  SD = stats::sd(x), 
                                                  N = length(x),
                                                  SEM = stats::sd(x)/sqrt(length(x))))
  auc_summary <- do.call(data.frame, auc_summary)
  
  # Create the AUC analysis list
  auc_analysis <- list(
    individual = auc_data,
    summary = auc_summary
  )
}


#' Analyze Tumor Growth Using Various Statistical Methods
#'
#' @importFrom utils head tail
#' @importFrom stats coef
#' @importFrom dplyr group_by summarize mutate arrange filter
#'
#' @description
#' This function provides a comprehensive statistical analysis for tumor growth data using
#' different statistical approaches. It supports multiple methods including linear mixed-effects
#' models, generalized additive mixed models, and area under the curve analysis.
#'
#' @param df A data frame containing tumor growth data. Must include columns for time, volume,
#'        treatment group, cage identifier, and individual subject ID. Additional columns can be included.
#' @param time_column A character string specifying the column name for time points (e.g., "Day"). 
#'        This column should contain numeric values representing time since treatment started. Default is "Day".
#' @param volume_column A character string specifying the column name for tumor volume measurements (e.g., "Volume").
#'        This column should contain numeric tumor volume measurements. Default is "Volume".
#' @param treatment_column A character string specifying the column name for treatment groups (e.g., "Treatment").
#'        This column should contain categorical treatment identifiers. Default is "Treatment".
#' @param cage_column A character string specifying the column name for the cage identifier (e.g., "Cage").
#'        This column should contain cage identifiers. The function automatically tests for collinearity with treatment.
#'        Default is "Cage".
#' @param id_column A character string specifying the column name for individual subject identifiers (e.g., "ID").
#'        This column should contain unique identifiers for each animal/subject. Default is "ID".
#' @param dose_column Optional. A character string specifying the column name for dose levels, if available.
#'        This allows for dose-response analysis. Default is NULL.
#' @param transform A character string specifying the transformation to apply to volume data.
#'        Options are "log", "sqrt", "none". Default is "log", which is recommended for exponential growth.
#' @param model_type A character string specifying the type of model to fit. Options are:
#'        "lme4" (standard linear mixed effects model using lme4 package — linear time × treatment),
#'        "gam"  (generalized additive mixed model via gamm4 — group-specific smooth time terms;
#'                use this in place of fitting a high-order polynomial in time),
#'        "auc"  (area under the curve analysis).
#'        Default is "lme4".
#' @param random_effects_specification A character string specifying the random effects structure.
#'        "intercept_only" (default): (1|ID) - random intercepts by subject
#'        "slope": (Day|ID) - random intercepts and slopes by subject
#'        "none": No random effects
#' @param handle_cage_effects How cage enters the model. Default \code{"auto"}
#'   picks the correct treatment for the detected design structure:
#'   \itemize{
#'     \item \strong{crossed} (cages hold more than one treatment) — cage as a
#'       fixed effect;
#'     \item \strong{nested with replication} (one treatment per cage, >= 2
#'       cages per arm — the usual preclinical design) — cage as a random
#'       intercept \code{(1|cage)}. A fixed effect is not estimable here
#'       because cage is perfectly aliased with treatment, but the random
#'       intercept is, and omitting it understates every standard error by the
#'       design effect \eqn{1 + (m - 1)\rho} for \eqn{m} mice per cage and
#'       cage-level ICC \eqn{\rho};
#'     \item \strong{completely confounded} (one cage per arm) — cage omitted,
#'       with a warning that any cage effect is inseparable from the treatment
#'       effect.
#'   }
#'   The explicit options \code{"always_include"} (fixed), \code{"never_include"},
#'   and \code{"as_random_effect"} are honoured where estimable and rejected
#'   where not. \code{"include_if_not_collinear"} is retained for backward
#'   compatibility: it reproduces the pre-0.6.0 behaviour of dropping cage
#'   whenever it is nested, which discards estimable cage variance.
#'   The detected structure, the placement chosen, and the reason are returned
#'   in \code{summary$model_specification}; the cage ICC is in
#'   \code{cage_analysis$icc}.
#' @param auc_method Removed in v0.4.5. The AUC path now always uses the
#'   trapezoidal rule. For last-observation-carried-forward (LOCF) AUC, use
#'   \code{\link{tumor_auc_analysis}(method = "last_observation")} directly.
#' @param comparison_family Which set of comparisons to report and adjust over:
#'   "vs_reference" (default; each treatment against \code{reference_group}),
#'   "all_pairs" (every pairwise comparison), or "custom" (supply
#'   \code{custom_contrasts}). The multiplicity adjustment is applied over
#'   exactly this family, so the family and the correction always agree.
#' @param custom_contrasts Named list of contrast coefficient vectors over the
#'   treatment levels, required when \code{comparison_family = "custom"}. The
#'   names label the comparisons in the output. Not supported for
#'   \code{model_type = "auc"}, which runs independent Welch t-tests rather
#'   than fitting a joint model.
#' @param p_adjust_method Multiplicity adjustment applied across the comparison
#'   family: "bonferroni" (default, conservative), "holm" (step-down, less
#'   conservative), "fdr" (Benjamini-Hochberg), "dunnett" (exact many-to-one
#'   via the multivariate-t; requires \code{comparison_family = "vs_reference"}),
#'   "tukey" (HSD; requires \code{comparison_family = "all_pairs"}), or "none".
#'   "dunnett" and "tukey" need comparisons from a single fitted model and are
#'   therefore rejected for \code{model_type = "auc"}. The family and method
#'   actually applied are returned in \code{comparison_family} and
#'   \code{p_adjust_method_used}.
#' @param reference_group Optional. A character string specifying which treatment group should be used as the reference
#'        for statistical comparisons. If NULL, the first treatment group alphabetically will be used.
#' @param return_model Boolean. Should the full fitted model be returned? Default is TRUE.
#' @param include_diagnostics Boolean. Should model diagnostic information be included? Default is TRUE.
#' @param plots Boolean. Should standard plots be generated and returned? Default is TRUE.
#' @param verbose Boolean. Should detailed information be printed during analysis? Default is FALSE.
#' @param extrapolation_points Number of points to extrapolate for each subject (default: 0)
#' @param auc_bootstrap_n Integer >= 0. When > 0 and \code{model_type = "auc"},
#'   each pairwise Welch's t-test gains a non-parametric percentile-bootstrap
#'   95\% CI for the mean difference (\code{boot_ci_lower}, \code{boot_ci_upper}
#'   columns in \code{pairwise_comparisons}). Honest at small N when per-mouse
#'   AUC distributions are skewed (especially under extrapolation), where the
#'   parametric t-CI is only approximate. A typical value is \code{999};
#'   higher (e.g. \code{4999}) gives smoother percentile estimates. Default
#'   \code{0} (skip; \code{boot_ci_*} columns are \code{NA}).
#' @param auc_bootstrap_seed Optional integer seed for reproducible bootstrap
#'   resampling. \code{NULL} (default) uses the caller's RNG state.
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{model}{The fitted statistical model (lme4 or auc object)}
#'   \item{anova}{ANOVA table for fixed effects with Type III tests}
#'   \item{summary}{Detailed model summary with parameter estimates}
#'   \item{pairwise_comparisons}{Results of pairwise comparisons between treatment groups, showing treatment differences, p-values, and confidence intervals}
#'   \item{treatment_effects}{Estimated treatment effects on tumor growth, showing adjusted means for each treatment group}
#'   \item{growth_rates}{Data frame containing growth rates for each subject, calculated as the slope of log-volume over time. Higher values indicate faster tumor growth.}
#'   \item{cage_analysis}{Analysis of cage effects, including:
#'     \itemize{
#'       \item{collinearity_test}{Chi-squared test result for collinearity between cage and treatment}
#'       \item{effects}{Data frame of cage-level statistics (mean and SD of volume by cage and treatment)}
#'     }
#'   }
#'   \item{model_selection}{Results of model selection process, including:
#'     \itemize{
#'       \item{aic}{AIC values for each model specification}
#'       \item{bic}{BIC values for each model specification}
#'       \item{selected_model}{Name of the model with lowest BIC (most parsimonious)}
#'     }
#'   }
#'   \item{diagnostics}{Model diagnostic information including residuals, random effects, and variance components}
#'   \item{auc_data}{Area under the curve data and statistics (when model_type="auc")}
#'   \item{plots}{List of standard plots for visualizing results}
#'   \item{data_summary}{Descriptive statistics of the processed data by treatment and time point, including mean, SD, and SEM of volumes}
#' }
#'
#' @details
#' This function offers several statistical approaches for analyzing tumor growth data. The choice of model
#' depends on the experimental design, growth patterns, and research questions:
#'
#' 1. Linear mixed-effects models (model_type="lme4") are suitable for most tumor growth experiments and
#'    account for the correlation structure of repeated measurements on the same subjects. 
#'    The default formula is log(Volume) ~ Day * Treatment + (1|ID).
#'
#' 2. Area under the curve analysis (model_type="auc") reduces each subject's growth curve to a single
#'    summary metric (AUC), which is then compared between treatment groups.
#'
#' The function automatically handles:
#' - Missing data patterns using appropriate methods
#' - Checking for and addressing cage effects
#' - Transformations to normalize growth curves
#' - Polynomial terms for non-linear growth patterns
#' - Diagnostic checks of model assumptions
#' - Growth rate analysis for each treatment group
#' - Cage effect analysis and collinearity testing
#' - Model selection based on AIC/BIC criteria
#' - Treatment effect estimation with proper reference group handling
#'
#' @examples
#' # Load example data
#' data(combo_treatment_synthetic_data)
#' tumor_data <- calculate_volume(combo_treatment_synthetic_data)
#' tumor_data <- calculate_dates(tumor_data, start_date = "03/24/2025")
#' 
#' # Basic analysis with default settings
#' results <- tumor_growth_statistics(tumor_data)
#' 
#' # View ANOVA table to assess treatment effects
#' print(results$anova)
#' 
#' # View pairwise comparisons
#' print(results$pairwise_comparisons)
#' 
#' # Plot adjusted means for each treatment group
#' print(results$plots$adjusted_means)
#' 
#' # Analyze with a generalized additive mixed model (group-specific smooth
#' # time terms via gamm4) — preferred over a high-order polynomial in time.
#' results_gam <- tumor_growth_statistics(
#'   tumor_data,
#'   model_type = "gam",
#'   transform = "log"
#' )
#' print(results_gam$treatment_effects_over_time)
#'
#' @export
tumor_growth_statistics <- function(df,
                                  time_column = "Day",
                                  volume_column = "Volume",
                                  treatment_column = "Treatment",
                                  cage_column = "Cage",
                                  id_column = "ID",
                                  dose_column = NULL,
                                  transform = c("log", "sqrt", "none"),
                                  model_type = c("lme4", "gam", "auc"),
                                  random_effects_specification = c("intercept_only", "slope", "none"),
                                  handle_cage_effects = c("auto", "include_if_not_collinear",
                                                          "always_include", "never_include",
                                                          "as_random_effect"),
                                  comparison_family = c("vs_reference", "all_pairs", "custom"),
                                  custom_contrasts = NULL,
                                  p_adjust_method = c("bonferroni", "holm", "fdr",
                                                      "dunnett", "tukey", "none"),
                                  reference_group = NULL,
                                  return_model = TRUE,
                                  include_diagnostics = TRUE,
                                  plots = TRUE,
                                  verbose = FALSE,
                                  extrapolation_points = 0,
                                  necrotic_column  = NULL,
                                  necrotic_handling = c("exclude", "covariate", "none"),
                                  auc_bootstrap_n = 0L,
                                  auc_bootstrap_seed = NULL) {
  # Check for required packages
  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("Please install the lme4 package: install.packages('lme4')")
  }
  
  # Match arguments
  transform <- match.arg(transform)
  model_type <- match.arg(model_type, c("lme4", "gam", "auc"))
  random_effects_specification <- match.arg(random_effects_specification)
  handle_cage_effects <- match.arg(handle_cage_effects)
  comparison_family <- match.arg(comparison_family)
  p_adjust_method <- match.arg(p_adjust_method)
  necrotic_handling <- match.arg(necrotic_handling)

  # CODE_REVIEW.md R3.1 / G.1 — resolve the comparison family and adjustment
  # once, up front, so all three model paths use the same validated spec and
  # the same set of comparisons that the adjustment is calibrated for.
  # The AUC path's Welch t-tests are not from a joint model, so Tukey/Dunnett
  # are rejected there.
  comparison_spec <- resolve_comparison_spec(
    comparison_family, p_adjust_method, custom_contrasts,
    supports_joint = model_type %in% c("lme4", "gam")
  )
  
  if (isTRUE(verbose)) {
    message("Analyzing tumor growth data...")
    message("Model type: ", model_type)
    message("Transform: ", transform)
  }
  
  # Check if required columns exist (cage_column is optional — NULL means no cage data)
  required_cols <- c(time_column, volume_column, treatment_column, id_column)
  missing_cols <- required_cols[!required_cols %in% colnames(df)]
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # When cage_column is NULL or not present in the data, create a synthetic
  # single-value placeholder so all downstream [[cage_column]] references work
  # without modification.  This also suppresses spurious cage-effect analysis.
  no_cage_mode <- is.null(cage_column) || !cage_column %in% colnames(df)
  if (no_cage_mode) {
    cage_column <- ".cage_placeholder"
    df[[cage_column]] <- "1"
    if (isTRUE(verbose)) message("No cage column provided 2014 cage effects will not be analysed.")
  }
  
  # Set reference group if not specified
  treatment_groups <- unique(df[[treatment_column]])
  if (is.null(reference_group)) {
    reference_group <- sort(treatment_groups)[1]
  } else if (!reference_group %in% treatment_groups) {
    stop("Reference group '", reference_group, "' is not present in the data.")
  }
  if (isTRUE(verbose)) {
    message("Using ", reference_group, " as reference group for statistical comparisons")
  }
  
  # Extrapolate data points if requested
  if (extrapolation_points > 0) {
    df <- tgs_extrapolate(df, id_column, treatment_column, cage_column,
                          time_column, volume_column, verbose)
  }
  
  # Create a copy of the data for analysis
  analysis_df <- df
  
  # Create a copy for AUC calculation before transformation
  auc_df <- df

  # Handle necrotic observations (before transformation, after extrapolation)
  necrosis_summary <- NULL
  include_necrotic_covariate <- FALSE
  if (!is.null(necrotic_column) && necrotic_column %in% colnames(analysis_df)) {
    nec <- tgs_handle_necrosis(analysis_df, necrotic_column, necrotic_handling,
                               treatment_column, id_column, time_column)
    analysis_df      <- nec$data
    necrosis_summary <- nec$necrosis_summary
    include_necrotic_covariate <- necrotic_handling == "covariate"
    if (isTRUE(verbose)) {
      n_excl <- attr(analysis_df, ".n_necrotic_excluded") %||% 0L
      message("Necrosis handling: '", necrotic_handling, "'. Excluded: ", n_excl, " observations.")
    }
  }

  # Apply transformations if needed
  if (transform == "log") {
    # Replace zeros/negatives with half the smallest positive value.
    # CODE_REVIEW.md R3.20 — guard the all-non-positive case. min(numeric(0))
    # is Inf, so without this the fill silently becomes Inf and log(Inf) = Inf
    # corrupts the fit. The Bayesian path was hardened against this in B3.3;
    # the frequentist entry point was not.
    vol           <- analysis_df[[volume_column]]
    positive_vals <- vol[is.finite(vol) & vol > 0]
    if (length(positive_vals) == 0L) {
      stop("No positive values in '", volume_column,
           "' after filtering. Cannot apply log transform.", call. = FALSE)
    }
    vol[vol <= 0] <- min(positive_vals, na.rm = TRUE) / 2
    analysis_df[[volume_column]] <- log(vol)
    if (isTRUE(verbose)) message("Applied log transformation to volume data")
  } else if (transform == "sqrt") {
    analysis_df[[volume_column]] <- sqrt(analysis_df[[volume_column]])
    if (isTRUE(verbose)) message("Applied square root transformation to volume data")
  }
  
  # Growth rate analysis
  growth_rates <- tgs_compute_growth_rates(auc_df, treatment_column, id_column,
                                           cage_column, time_column, volume_column)

  # Cage effect analysis
  cage_analysis <- tgs_compute_cage_effects(analysis_df, cage_column, treatment_column, volume_column)
  cage_effects <- cage_analysis$effects

  # Model fitting and selection.
  #
  # CODE_REVIEW.md R3.17 / G.2 — cage inclusion is now decided by the design's
  # STRUCTURE (does any cage hold more than one treatment; are there >= 2 cages
  # per arm), not by a chi-square p-value on observation counts. The old test
  # dropped cage entirely in the maintainer's standard design (one arm per
  # cage, >= 2 cages per arm) even though (1|cage) is estimable and needed
  # there — understating every standard error by the design effect.
  cage_structure <- classify_cage_structure(
    analysis_df, cage_column, id_column, treatment_column
  )
  cage_analysis$structure <- cage_structure
  cage_handling <- resolve_cage_handling(cage_structure, handle_cage_effects,
                                         verbose = verbose)
  cage_analysis$handling <- cage_handling

  # Retained for the descriptive chi-square in cage_analysis; no longer drives
  # model structure.
  cage_collinear <- !is.null(cage_analysis$collinearity_test) &&
    isTRUE(cage_analysis$collinearity_test$p.value < 0.05)

  # ── GAM path (early return) ──────────────────────────────────────────────
  # Use gamm4 (lme4 + mgcv) when model_type == "gam". Skips the LME4 fit
  # entirely and builds a result list that mirrors the LME4 shape so the
  # dashboard render pipeline is shared.
  if (model_type == "gam") {
    # data_summary is computed inside the GAM branch (the LME4 path computes
    # it later); same helper, same shape.
    data_summary <- tgs_compute_summary(
      analysis_df, treatment_column, time_column, volume_column
    )

    gam_result <- tgs_fit_gamm4_model(
      analysis_df, volume_column, time_column,
      treatment_column, id_column, cage_column,
      cage_handling, verbose,
      include_necrotic_covariate = include_necrotic_covariate
    )
    if (is.null(gam_result)) {
      gam_err <- attr(gam_result, "gamm4_error")
      stop("model_type = 'gam' failed to fit",
           if (!is.null(gam_err)) paste0(": ", gam_err) else ".",
           call. = FALSE)
    }

    gam_fit  <- gam_result$model
    gam_obj  <- gam_fit$gam
    mer_obj  <- gam_fit$mer

    day_range <- sort(unique(analysis_df[[time_column]]))
    mean_day  <- mean(analysis_df[[time_column]])

    treatment_effects   <- tgs_gam_treatment_effects(
      gam_obj, treatment_column, time_column, mean_day, reference_group
    )
    emm_time_df         <- tgs_gam_emm_time(
      gam_obj, treatment_column, time_column, day_range
    )
    pairwise_comp_df    <- tgs_gam_pairwise(
      gam_obj, treatment_column, time_column, day_range, reference_group,
      comparison_spec = comparison_spec, custom_contrasts = custom_contrasts
    )
    anova_table         <- tgs_gam_anova_table(gam_obj)
    gam_diagnostics     <- if (include_diagnostics) {
      tgs_gam_diagnostics(gam_fit, id_column)
    } else NULL

    analysis_summary <- list(
      analysis_type = "Generalized Additive Mixed Model (gamm4)",
      data_description = list(
        subjects         = length(unique(make_mouse_key(
          analysis_df[[id_column]],
          analysis_df[[treatment_column]],
          analysis_df[[cage_column]]
        ))),
        treatment_groups = length(unique(analysis_df[[treatment_column]])),
        time_points      = length(unique(analysis_df[[time_column]])),
        reference_group  = reference_group
      ),
      model_specification = list(
        fixed_effects  = paste(volume_column, "~", treatment_column,
                               "+ s(", time_column, ", by =",
                               treatment_column, ")"),
        random_effects = if (isTRUE(cage_handling$random)) {
          paste0("(1|", cage_column, ") + (1|", id_column, ")")
        } else {
          paste0("(1|", id_column, ")")
        },
        cage_effects   = handle_cage_effects,
        cage_structure = cage_structure$structure,
        cage_placement = if (isTRUE(cage_handling$fixed)) "fixed effect"
                         else if (isTRUE(cage_handling$random)) "random intercept"
                         else "omitted",
        cage_reason    = cage_handling$reason,
        smoother_k     = gam_result$model_selection$k_basis
      ),
      model_selection = gam_result$model_selection,
      methods = list(
        volume_transformation = transform,
        anova_method          = "Smooth-term significance (mgcv summary)",
        posthoc_method        = paste0(
          "Difference of group-specific smooths at five study-day quantiles ",
          "(min / Q1 / median / Q3 / max), pairwise vs ", reference_group
        )
      ),
      notes = c(
        if (transform != "none") paste("Volume data was", transform,
                                       "transformed prior to analysis")
        else "No transformation applied to volume data",
        paste0("Smoother basis dimension k = ",
               gam_result$model_selection$k_basis,
               "; check k_check in diagnostics if k.index << 1")
      )
    )

    return(list(
      model                       = if (return_model) gam_fit else NULL,
      model_type_used             = "gam",
      transform_used              = transform,   # CODE_REVIEW.md R3.29
      comparison_family     = comparison_spec$family,
      p_adjust_method_used  = comparison_spec$p_adjust_method,
      anova                       = anova_table,
      summary                     = analysis_summary,
      pairwise_comparisons        = pairwise_comp_df,
      posthoc = list(
        method   = "GAM smooth-difference contrasts at quantile days",
        pairwise = pairwise_comp_df
      ),
      treatment_effects           = treatment_effects,
      treatment_effects_over_time = emm_time_df,
      growth_rates                = growth_rates,
      cage_analysis               = cage_analysis,
      model_selection             = gam_result$model_selection,
      diagnostics                 = gam_diagnostics,
      # Promote GAM-specific diagnostics to top-level fields so the
      # dashboard can read them with the same shape as LME4 / AUC.
      diag_qq_plot                = if (!is.null(gam_diagnostics)) gam_diagnostics$diag_qq_plot,
      diag_resid_fitted_plot      = if (!is.null(gam_diagnostics)) gam_diagnostics$diag_resid_fitted_plot,
      diag_scale_location_plot    = if (!is.null(gam_diagnostics)) gam_diagnostics$diag_scale_location_plot,
      diag_re_qq_plot             = NULL,    # gamm4's random-effects are inside $mer; not extracted here
      gam_k_check                 = if (!is.null(gam_diagnostics)) gam_diagnostics$k_check,
      gam_concurvity              = if (!is.null(gam_diagnostics)) gam_diagnostics$concurvity,
      gam_deviance_explained      = if (!is.null(gam_diagnostics)) gam_diagnostics$deviance_explained,
      data_summary                = data_summary,
      plots                       = NULL,
      necrosis_summary            = necrosis_summary
    ))
  }

  lme4_result <- tgs_fit_lme4_models(
    analysis_df, volume_column, time_column,
    treatment_column, id_column, cage_column,
    random_effects_specification, cage_handling, verbose,
    include_necrotic_covariate = include_necrotic_covariate
  )
  model <- lme4_result$model
  model_selection <- lme4_result$model_selection
  best_model <- lme4_result$best_model

  # Diagnostic plots
  diagnostics <- list()
  
  if (include_diagnostics) {
    # Residual plots
    diagnostics$residuals <- list(
      fitted = stats::fitted(model),
      residuals = stats::residuals(model),
      qq_plot = stats::qqnorm(stats::residuals(model), plot.it = FALSE)
    )
    
    # Random effects plots (only available for lmer models)
    if (inherits(model, "lmerMod")) {
      diagnostics$random_effects <- list(
        intercepts = lme4::ranef(model)[[id_column]],
        slopes = if (best_model == "slope") {
          lme4::ranef(model)[[id_column]][, 2]
        } else NULL
      )
      diagnostics$variance_components <- lme4::VarCorr(model)
    } else {
      diagnostics$random_effects <- NULL
      diagnostics$variance_components <- NULL
    }
  }

  # Create a basic summary of the data
  data_summary <- tgs_compute_summary(analysis_df, treatment_column, time_column, volume_column)

  # Only for auc model type or when additional AUC analysis is requested
  if (model_type == "auc" || (model_type == "lme4" && include_diagnostics)) {
    auc_analysis <- tgs_compute_auc(auc_df, id_column, treatment_column, cage_column,
                                    time_column, volume_column, verbose)
  }

  # Return the results for AUC model
  if (model_type == "auc") {
    # CODE_REVIEW.md D.1 — the AUC path body lives in R/tgs_path_auc.R
    return(tgs_path_auc(
      auc_analysis        = auc_analysis,
      growth_rates        = growth_rates,
      auc_df              = auc_df,
      transform           = transform,
      comparison_spec     = comparison_spec,
      reference_group     = reference_group,
      return_model        = return_model,
      include_diagnostics = include_diagnostics,
      auc_bootstrap_n     = auc_bootstrap_n,
      auc_bootstrap_seed  = auc_bootstrap_seed,
      cage_analysis       = cage_analysis,
      data_summary        = data_summary,
      necrosis_summary    = necrosis_summary,
      id_column           = id_column,
      treatment_column    = treatment_column,
      cage_column         = cage_column,
      time_column         = time_column
    ))
  } else {
    # Create ANOVA table
    if (requireNamespace("car", quietly = TRUE)) {
      anova_table <- car::Anova(model, type = "III")
      anova_type <- "Type III ANOVA (car package)"
    } else {
      anova_table <- stats::anova(model)
      anova_type <- "ANOVA (stats package)"
      warning("Package 'car' not available. Using stats::anova instead of Type III tests.")
    }
    
    # Create pairwise comparisons
    if (requireNamespace("emmeans", quietly = TRUE)) {
      # Marginalise at the mean time point for the primary treatment_effects
      # summary (backward-compatible). Additionally compute EMMs at the five
      # study-day quantiles (min, Q1, median, Q3, max) to capture how the
      # treatment × time interaction unfolds — marginalising at a single point
      # discards the interaction structure that is the primary analysis target.
      mean_day <- mean(analysis_df[[time_column]])
      day_range <- unique(sort(analysis_df[[time_column]]))
      quant_days <- unique(round(stats::quantile(
        day_range,
        probs = c(0, 0.25, 0.5, 0.75, 1),
        type  = 1
      )))

      lsmeans_obj <- emmeans::emmeans(
        model, specs = treatment_column,
        at = stats::setNames(list(mean_day), time_column)
      )

      # EMMs over time: one row per (Treatment, Day) combination
      lsmeans_time <- emmeans::emmeans(
        model,
        specs = c(treatment_column, time_column),
        at    = stats::setNames(list(quant_days), time_column)
      )
      emm_time_df <- as.data.frame(summary(lsmeans_time))
      names(emm_time_df)[names(emm_time_df) == "emmean"]   <- "Adjusted_Mean"
      names(emm_time_df)[names(emm_time_df) == "lower.CL"] <- "Lower_CL"
      names(emm_time_df)[names(emm_time_df) == "upper.CL"] <- "Upper_CL"

      # Extract treatment effects
      emm_summary <- summary(lsmeans_obj)
      treatment_effects <- data.frame(
        Group = emm_summary[[treatment_column]],
        Adjusted_Mean = emm_summary$emmean,
        SE = emm_summary$SE,
        DF = emm_summary$df,
        Lower_CL = emm_summary$lower.CL,
        Upper_CL = emm_summary$upper.CL
      )

      # Add reference group note
      treatment_effects$Note <- ifelse(treatment_effects$Group == reference_group, "Reference group", "")

      # Reorder treatment effects to put reference group first
      ref_idx <- which(treatment_effects$Group == reference_group)
      if (ref_idx > 1) {
        treatment_effects <- rbind(
          treatment_effects[ref_idx, ],
          treatment_effects[-ref_idx, ]
        )
      }

      # Format numeric columns
      treatment_effects$Adjusted_Mean <- round(treatment_effects$Adjusted_Mean, 3)
      treatment_effects$SE <- round(treatment_effects$SE, 3)
      treatment_effects$Lower_CL <- round(treatment_effects$Lower_CL, 3)
      treatment_effects$Upper_CL <- round(treatment_effects$Upper_CL, 3)

      # Calculate the requested comparisons WITH the requested adjustment.
      #
      # CODE_REVIEW.md R3.1 — this previously built a custom contrast list and
      # passed it to emmeans::contrast() with no `adjust`. emmeans defaults to
      # adjust = "none" for a user-supplied coefficient list (it only defaults
      # to tukey/dunnettx for the named "pairwise"/"trt.vs.ctrl" methods), so
      # the default model path reported completely unadjusted p-values while
      # `p_adjust_method` documented "bonferroni" as the default. The argument
      # was never referenced anywhere in this path.
      pairwise_comp <- build_requested_contrasts(
        lsmeans_obj, comparison_spec,
        reference_group  = reference_group,
        custom_contrasts = custom_contrasts
      )

      posthoc_method <- paste0(
        "Estimated marginal means at the mean study day; ",
        switch(comparison_spec$family,
          vs_reference = paste0("each group vs ", reference_group),
          all_pairs    = "all pairwise comparisons",
          custom       = "user-supplied contrasts"
        ),
        " with ", comparison_spec$p_adjust_method, " adjustment"
      )

    } else {
      pairwise_comp <- NULL
      treatment_effects    <- NULL
      emm_time_df          <- NULL
      posthoc_method <- NA
      warning("Package 'emmeans' not available. Pairwise comparisons and treatment effects not calculated.")
    }
    
    # Create plots
    if (plots) {
      plot_data <- data_summary
      plots_list <- list(
        data_summary = plot_data,
        growth_rates = growth_rates,
        cage_effects = cage_effects,
        diagnostics = diagnostics
      )
      
      # Add adjusted means plot if available
      if (!is.null(treatment_effects)) {
        plots_list$adjusted_means <- treatment_effects
      }
    } else {
      plots_list <- NULL
    }
    
    # Create a descriptive summary
    analysis_summary <- list(
      analysis_type = "Linear Mixed Effects Model Analysis",
      data_description = list(
        subjects = length(unique(make_mouse_key(
          analysis_df[[id_column]], analysis_df[[treatment_column]], analysis_df[[cage_column]]
        ))),
        treatment_groups = length(unique(analysis_df[[treatment_column]])),
        time_points = length(unique(analysis_df[[time_column]])),
        reference_group = reference_group
      ),
      model_specification = list(
        fixed_effects = paste(volume_column, "~", time_column, "*", treatment_column),
        random_effects = switch(random_effects_specification,
                               "intercept_only" = paste("(1|", id_column, ")"),
                               "slope" = paste("(", time_column, "|", id_column, ")"),
                               "none" = "None"),
        cage_effects   = handle_cage_effects,
        cage_structure = cage_structure$structure,
        cage_placement = if (isTRUE(cage_handling$fixed)) "fixed effect"
                         else if (isTRUE(cage_handling$random)) "random intercept"
                         else "omitted",
        cage_reason    = cage_handling$reason
      ),
      model_selection = list(
        criteria = paste("Random effects:", random_effects_specification),
        selected_model = if (!is.null(model_selection$selected_model)) {
          model_selection$selected_model
        } else random_effects_specification
      ),
      methods = list(
        volume_transformation = transform,
        anova_method = anova_type,
        posthoc_method = posthoc_method,
        growth_rate_calculation = paste0(
          "Growth rates are calculated by fitting a linear regression model to log-transformed volume data over time for each subject ",
          "(non-positive volumes are replaced with half the smallest positive value before the log). ",
          "The slope coefficient from this model represents the exponential growth rate. ",
          "A value of 0.1 indicates approximately 10% tumor volume increase per day. ",
          "Only subjects with 3 or more time points are included in growth rate calculations."
        )
      ),
      notes = c(
        if(transform != "none") paste("Volume data was", transform, "transformed prior to analysis") else "No transformation applied to volume data",
        paste("Cage design:", cage_structure$description),
        paste("Cage handling:", cage_handling$reason)
      )
    )
    
    # CODE_REVIEW.md DIAGNOSTICS bug (1) — build the ggplots the dashboard's
    # TG Diagnostics tab expects (diag_qq_plot, diag_resid_fitted_plot,
    # diag_scale_location_plot, diag_re_qq_plot). Previously the function
    # only stored qqnorm() *data* under diagnostics$residuals; the dashboard
    # looked for rendered ggplots on the top-level result list and always
    # showed "not available". Mirrors the analyze_body_weight pattern.
    rd <- if (include_diagnostics)
      build_residual_diagnostic_plots(model,
                                      title_prefix = "Tumour growth LMM")
    else
      list(diag_qq_plot = NULL,
           diag_resid_fitted_plot = NULL,
           diag_scale_location_plot = NULL)
    diag_re_qq_plot <- if (include_diagnostics)
      build_random_effects_qq_plot(model, id_column,
                                   title_prefix = "Tumour growth LMM")
    else NULL
    # CODE_REVIEW.md DIAGNOSTICS gap (13) — LMM influence diagnostics
    # (Cook's distance + DFBETAS). Refits the model leaving each observation
    # out; can be slow on large designs, so gated on include_diagnostics.
    lmm_infl <- if (include_diagnostics) build_lmm_influence(model) else NULL

    # CODE_REVIEW.md R3.17 / G.2 — report the cage-level intraclass correlation
    # so the reader can see how much the clustering mattered. The design effect
    # on the effective sample size is 1 + (mice_per_cage - 1) * ICC.
    cage_analysis$icc <- if (isTRUE(cage_handling$random)) {
      cage_icc(model, cage_column)
    } else NULL

    # Return the results
    results <- list(
      model = if (return_model) model else NULL,
      # CODE_REVIEW.md R3.29 — report the modelling scale and model type so a
      # downstream consumer can back-transform correctly without having to know
      # what was requested.
      model_type_used = "lme4",
      transform_used  = transform,
      comparison_family     = comparison_spec$family,
      p_adjust_method_used  = comparison_spec$p_adjust_method,
      anova = anova_table,
      summary = analysis_summary,
      pairwise_comparisons = pairwise_comp,
      posthoc = list(
        method = posthoc_method,
        pairwise = if (!is.null(pairwise_comp)) as.data.frame(summary(pairwise_comp)) else NULL
      ),
      treatment_effects          = treatment_effects,
      treatment_effects_over_time = emm_time_df,
      growth_rates = growth_rates,
      cage_analysis = cage_analysis,
      model_selection = model_selection,
      diagnostics = if (include_diagnostics) diagnostics else NULL,
      diag_qq_plot             = rd$diag_qq_plot,
      diag_resid_fitted_plot   = rd$diag_resid_fitted_plot,
      diag_scale_location_plot = rd$diag_scale_location_plot,
      diag_re_qq_plot          = diag_re_qq_plot,
      diag_cooks_distance      = if (!is.null(lmm_infl)) lmm_infl$cooks_distance,
      diag_dfbetas             = if (!is.null(lmm_infl)) lmm_infl$dfbetas,
      data_summary = data_summary,
      plots = plots_list,
      necrosis_summary = necrosis_summary
    )

    return(results)
  }
}