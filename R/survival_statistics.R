# Copyright (c) 2026 mouseExperiment Contributors
# Licensed under the MIT License - see LICENSE file

#' Perform Survival Analysis for Mouse Tumor Experiments
#'
#' @description
#' Performs comprehensive survival analysis using appropriate statistical methods
#' based on data characteristics. Automatically selects between Cox Proportional 
#' Hazards Model, Firth's bias-reduced estimation, or Log-Rank Test depending on
#' the presence of complete or quasi-complete separation in the data.
#'
#' @param df Data frame containing survival data, with **one row per animal**.
#'   Each row must carry that animal's time-to-event and event indicator. A
#'   longitudinal frame (one row per measurement occasion) is rejected with an
#'   error, because fitting a Cox model to it would treat every measurement as
#'   an independent subject at risk. Reduce to one row per animal first — e.g.
#'   each animal's last observation, carrying its event indicator.
#' @param time_column Name of column containing time-to-event data. Default: "Day"
#' @param censor_column Name of column containing censoring indicator (1=event, 0=censored). Default: "Survival_Censor"
#' @param treatment_column Name of column containing treatment groups. Default: "Treatment"
#' @param cage_column Name of column containing cage identifiers. Default: "Cage"
#' @param id_column Name of column containing individual mouse identifiers. Default: "ID"
#' @param randomisation_unit What treatment was assigned to: \code{"mouse"}
#'   (default, individual animals) or \code{"cage"} (whole cages). CODE_REVIEW.md
#'   R20.6: the function used to add \code{cluster(cage)} whenever cages were
#'   replicated, and a sandwich variance from 4-10 clusters gave a 28%
#'   false-positive rate under the null.
#'   \itemize{
#'     \item \code{"mouse"}: Cox (or Firth) with ordinary model-based standard
#'       errors; cage is not modelled.
#'     \item \code{"cage"}: each arm is compared with the reference by a
#'       cage-level permutation log-rank, which moves whole cages between the
#'       two arms (exact when there are at most 5,000 assignments). Hazard
#'       ratios are reported as point estimates without an interval, because
#'       no valid cage-level interval is available. The smallest p-value such a
#'       comparison can reach is reported: with 2 cages per arm it is 1/3, and
#'       with 3 it is 0.1, so neither can reach p < 0.05.
#'   }
#' @param dose_column Optional name of column containing dose information. Default: NULL
#' @param reference_group Treatment group to use as reference. Default: NULL (uses first alphabetically)
#' @param firth_correction Whether to apply Firth's correction for separation issues. Default: TRUE
#' @param p_adjust_method Multiplicity adjustment applied across the k-1
#'   treatment-vs-reference comparisons: "bonferroni" (default), "holm", "fdr"
#'   (Benjamini-Hochberg), or "none". The unadjusted values are retained in
#'   `P_Value_Unadjusted` and the method used in `P_Adjust_Method`.
#' @param permutation_logrank Logical. When the log-rank fallback is used (a
#'   group has zero events, so the Cox partial likelihood is not estimable), use
#'   an exact-style permutation log-rank via the \pkg{coin} package rather than
#'   the asymptotic chi-square approximation. That regime — few events, small n —
#'   is where the asymptotic test is least trustworthy and where permutation is
#'   valid without any large-sample appeal. Falls back to \code{survdiff()} with
#'   a message when \pkg{coin} is unavailable. Default TRUE.
#'   CODE_REVIEW.md H.4.
#' @param verbose Whether to print analysis details to the console. Default: TRUE
#'
#' @return A list containing:
#' \describe{
#'   \item{model}{The fitted statistical model object}
#'   \item{results}{Data frame with hazard ratios, confidence intervals, p-values, and median survival times}
#'   \item{reference_group}{The treatment group used as reference}
#'   \item{method_used}{The statistical method used ("cox", "coxphf", or "logrank")}
#'   \item{ph_test}{A \code{cox.zph} object (Schoenfeld-residual proportional-hazards
#'     test) when the standard Cox path is used; \code{NULL} for the Firth and
#'     log-rank fallbacks. A small global p-value indicates the PH assumption is
#'     violated and hazard ratios should be interpreted with caution.}
#' }
#'
#' @details
#' The function adapts to data characteristics:
#' * For well-behaved data: Standard Cox proportional hazards model
#' * For groups with no events: Log-rank test
#' * For groups with few events: Firth's bias-reduced Cox model
#'
#' @importFrom survival Surv survfit coxph survdiff
#' @importFrom stats confint model.matrix as.formula pchisq
#'
#' @examples
#' # master_synthetic_data is longitudinal (one row per measurement occasion)
#' # and carries a Survival_Censor column. Reduce it to one row per animal
#' # before calling — the function requires that and will error otherwise.
#' data(master_synthetic_data)
#'
#' per_animal <- do.call(rbind, lapply(
#'   split(master_synthetic_data,
#'         list(master_synthetic_data$ID, master_synthetic_data$Treatment),
#'         drop = TRUE),
#'   function(animal) {
#'     animal <- animal[order(animal$Day), ]
#'     last <- animal[nrow(animal), ]
#'     # an event anywhere in the record is an event for that animal
#'     last$Survival_Censor <- max(animal$Survival_Censor, na.rm = TRUE)
#'     last
#'   }
#' ))
#'
#' # Run survival analysis
#' results <- survival_statistics(
#'   df = per_animal,
#'   reference_group = "Vehicle"
#' )
#'
#' # Access results
#' print(results$results)  # Hazard ratios, CIs, p-values, and median survival times
#' 
#' # Extract median survival times
#' median_surv <- results$results$Median_Survival
#' names(median_surv) <- results$results$Group
#' print(median_surv)
#'
#' @export
survival_statistics <- function(df,
                              time_column = "Day",
                              censor_column = "Survival_Censor",
                              treatment_column = "Treatment",
                              cage_column = "Cage",
                              id_column = "ID",
                              randomisation_unit = c("mouse", "cage"),
                              dose_column = NULL,
                              reference_group = NULL,
                              firth_correction = TRUE,
                              p_adjust_method = c("bonferroni", "holm", "fdr", "none"),
                              permutation_logrank = TRUE,
                              verbose = TRUE) {

  p_adjust_method <- match.arg(p_adjust_method)
  randomisation_unit <- match.arg(randomisation_unit)

  # Validate inputs
  validate_inputs(df, time_column, censor_column, treatment_column)
  # R15.1: cage_column arrives either as NULL (the dashboard passes NULL when no
  # cage column is mapped) or as a name that may not be present -- the default is
  # "Cage", which is absent from plenty of real uploads. Both cases used to fail
  # opaquely: NULL gave "argument is of length zero" from an unguarded
  # `if (cage_column %in% colnames(df))`, and a missing name gave "undefined
  # columns selected" from a bare data-frame subset. There was no way to run this
  # function on data without a cage column, and the dashboard's Survival tab hit
  # the NULL branch for every such upload.
  #
  # Normalise here so everything below can assume NULL means "no cage information".
  if (!is.null(cage_column) && !cage_column %in% colnames(df)) cage_column <- NULL

  # CODE_REVIEW.md R20.15 -- the ID column is required. The dashboard never
  # passed id_column, so any upload whose ID column was not literally "ID" failed
  # with "undefined columns selected".
  if (is.null(id_column) || !id_column %in% colnames(df)) {
    stop("ID column '", if (is.null(id_column)) "NULL" else id_column,
         "' not found. Pass the column ",
         "that identifies each animal as id_column.", call. = FALSE)
  }

  # CODE_REVIEW.md R20.43 -- formulas below are built by pasting column names,
  # and some sites quoted them while others did not, so a column called
  # "Study Day" failed with "unexpected symbol". Copy the columns into fixed
  # internal names once, here.
  work <- data.frame(
    Time      = as.numeric(df[[time_column]]),
    Event     = df[[censor_column]],
    Treatment = df[[treatment_column]],
    ID        = as.character(df[[id_column]]),
    stringsAsFactors = FALSE
  )
  if (!is.null(cage_column)) work$Cage <- as.character(df[[cage_column]])
  df               <- work
  time_column      <- "Time"
  censor_column    <- "Event"
  treatment_column <- "Treatment"
  id_column        <- "ID"
  if (!is.null(cage_column)) cage_column <- "Cage"

  validate_one_row_per_subject(df, id_column, treatment_column, cage_column)
  
  # Setup parameters
  treatment_groups <- unique(df[[treatment_column]])
  if (is.null(reference_group)) {
    reference_group <- sort(treatment_groups)[1]
  } else if (!reference_group %in% treatment_groups) {
    stop("Reference group '", reference_group, "' is not present in the data.")
  }
  message("Using ", reference_group, " as the reference group for hazard ratios.")
  
  # Check cage distribution
  check_cage_distribution(df, treatment_column, cage_column)

  # CODE_REVIEW.md R3.13 / G.2 — cage was accepted, printed about, and then
  # never entered the model: the Cox fit was `Surv(...) ~ Treatment` with no
  # cluster() or frailty() term. Co-housed animals share an environment, so
  # ignoring that correlation makes the standard errors and p-values
  # anti-conservative wherever cage is not perfectly confounded with treatment.
  # Use the same structural classification as tumor_growth_statistics(): a
  # robust sandwich variance via cluster() when cage is nested with replication
  # or crossed, and nothing (with a warning) when it is completely confounded,
  # since there is no replication to estimate from.
  cage_structure <- classify_cage_structure(df, cage_column, id_column,
                                            treatment_column)
  # CODE_REVIEW.md R20.6 -- this used cluster(cage) whenever cages were
  # replicated or crossed. A sandwich variance from 4-10 clusters is badly
  # downward-biased: 28 % false positives with 2 cages x 5 mice per arm. The
  # declared unit of randomisation decides instead (see randomisation_unit).
  use_cage_cluster <- FALSE
  if (cage_structure$structure == "nested_confounded") {
    warning("Cage and treatment are completely confounded: ",
            cage_structure$description,
            " Hazard ratios therefore include any cage effect.",
            call. = FALSE)
  }
  if (randomisation_unit == "cage") {
    if (is.null(cage_column)) {
      stop("randomisation_unit = 'cage' needs a cage column.", call. = FALSE)
    }
    if (cage_structure$structure == "crossed") {
      stop("randomisation_unit = 'cage', but some cages hold more than one ",
           "treatment, so treatment cannot have been assigned to whole cages. ",
           "Check the cage column, or declare randomisation_unit = 'mouse'.",
           call. = FALSE)
    }
  }
  if (isTRUE(verbose)) {
    message("Cage structure: ", cage_structure$structure, ". Randomisation unit: ",
            randomisation_unit, ".")
  }
  
  # Check for separation issues
  surv_obj <- survival::Surv(df[[time_column]], df[[censor_column]])
  surv_formula_str <- paste("surv_obj ~", treatment_column)
  cox_formula <- stats::as.formula(surv_formula_str)
  separation_info <- check_separation(df, treatment_column, censor_column)
  
  # Choose and fit appropriate model
  model_results <- fit_survival_model(
    df,
    surv_obj,
    cox_formula,
    treatment_column,
    treatment_groups,
    reference_group,
    time_column,
    censor_column,
    separation_info,
    firth_correction,
    verbose = verbose,
    p_adjust_method = p_adjust_method,
    cage_column = NULL,   # no cluster() term (R20.6)
    permutation_logrank = permutation_logrank
  )
  
  # Extract model results
  model <- model_results$model
  results <- model_results$results
  method_used <- model_results$method_used

  # R20.6 -- under cage randomisation the p-values come from permuting whole
  # cages, which is the design-faithful test. The Cox standard errors assume
  # independent animals, so their intervals are withheld.
  cage_permutation <- NULL
  if (randomisation_unit == "cage" && nrow(results) > 0L) {
    perm_rows <- list()
    for (g in setdiff(as.character(treatment_groups), reference_group)) {
      cp  <- cage_permutation_logrank(df, g, reference_group)
      idx <- which(results$Group == g)
      results$P_Value[idx]  <- cp$p_value
      results$CI_Lower[idx] <- NA_real_
      results$CI_Upper[idx] <- NA_real_
      perm_rows[[g]] <- data.frame(
        Group = g, Cages_Treated = cp$cages_trt, Cages_Reference = cp$cages_ref,
        Assignments = cp$n_assignments, Exact = cp$exact, P_Value = cp$p_value,
        Min_Attainable_P = cp$min_attainable_p, stringsAsFactors = FALSE)
    }
    cage_permutation <- do.call(rbind, perm_rows)
    rownames(cage_permutation) <- NULL
    results$P_Method <- ifelse(results$Group == reference_group, NA_character_,
                               "cage-level permutation log-rank")
    floor_hit <- cage_permutation$Min_Attainable_P > 0.05
    if (any(floor_hit)) {
      warning("With so few cages, ",
              paste(sprintf("%s vs %s cannot reach p < %.3g", cage_permutation$Group[floor_hit],
                            reference_group, cage_permutation$Min_Attainable_P[floor_hit]),
                    collapse = "; "),
              " at any effect size: a cage-randomised comparison has only as many ",
              "distinct outcomes as ways to assign the cages.", call. = FALSE)
    }
  }

  # CODE_REVIEW.md R3.1 / G.1 — apply the multiplicity adjustment once, here,
  # rather than in each of the three model branches. The comparison family is
  # the k-1 treatment-vs-reference contrasts actually reported, so the
  # adjustment is over exactly the set of tests in the returned table. Keep the
  # unadjusted values and the method used so the result is self-describing.
  if (nrow(results) > 0L && "P_Value" %in% names(results)) {
    non_ref <- results$Group != reference_group
    results$P_Value_Unadjusted <- results$P_Value
    if (any(non_ref)) {
      results$P_Value[non_ref] <- stats::p.adjust(results$P_Value[non_ref],
                                                  method = p_adjust_method)
    }
    results$P_Adjust_Method <- p_adjust_method
    results$Comparison_Family <- "vs_reference"
  }
  
  # Create a separate survival fit for median survival calculation
  message("\nCalculating median survival times...")
  surv_formula_str <- paste("Surv(", time_column, ",", censor_column, ") ~ ", 
                            treatment_column)
  surv_formula <- stats::as.formula(surv_formula_str)
  km_fit <- survival::survfit(surv_formula, data = df)
  
  # Display median survival information
  if (isTRUE(verbose)) message(paste(utils::capture.output(print(km_fit)), collapse = "\n"))
  
  # Calculate and add median survival times
  median_survival <- NULL
  tryCatch({
    fit_summary <- summary(km_fit)
    if ("table" %in% names(fit_summary) && is.matrix(fit_summary$table) && 
        "median" %in% colnames(fit_summary$table)) {
      median_survival <- fit_summary$table[, "median"]
      names(median_survival) <- rownames(fit_summary$table)
      
    } else if (!is.null(km_fit$median)) {
      median_survival <- km_fit$median
      names(median_survival) <- names(km_fit$strata)
    }
    
    if (!is.null(median_survival)) {
      # Clean up strata names
      if (!is.null(names(median_survival))) {
        names(median_survival) <- gsub(paste0(treatment_column, "="), "", names(median_survival))
      }
      
      # Add to results data frame
      results$Median_Survival <- median_survival[match(results$Group, names(median_survival))]
      
      # Display the median survival times
      message("\nMedian Survival Times:")
      for (i in seq_along(median_survival)) {
        group_name <- names(median_survival)[i]
        med_surv_val <- median_survival[i]
        if (!is.na(med_surv_val)) {
          message(sprintf("%s: %.1f days", group_name, med_surv_val))
        } else {
          message(sprintf("%s: NA days", group_name))
        }
      }
    }
  }, error = function(e) {
    warning("Error calculating median survival: ", e$message)
  })
  
  # Add Events and Total columns
  # Improved approach to count unique subjects and their events per treatment group
  # First, find the unique subjects (including cage information) in each treatment group
  subject_treatment <- unique(df[, c(id_column, treatment_column, cage_column)])
  total_counts <- table(subject_treatment[[treatment_column]])
  
  # Next, find subjects with events
  # We need to handle possible duplicates in the data (multiple rows per subject)
  # For each subject, if any row has an event, count it as an event
  event_data <- df[, c(id_column, treatment_column, censor_column, cage_column)]
  # Aggregate to get maximum event per subject (1 if any event occurred, 0 otherwise)
  event_by_subject <- stats::aggregate(
    event_data[[censor_column]], 
    by = c(
      list(
        ID = event_data[[id_column]],
        Treatment = event_data[[treatment_column]]
      ),
      # Grouping by cage only when there is one; event_data[[NULL]] errors.
      if (!is.null(cage_column)) list(Cage = event_data[[cage_column]]) else NULL
    ), 
    FUN = max
  )
  
  # Count events per treatment group
  event_counts <- tapply(event_by_subject$x, event_by_subject$Treatment, sum)
  
  # Assign to results data frame
  results$Events <- event_counts[match(results$Group, names(event_counts))]
  results$Total <- total_counts[match(results$Group, names(total_counts))]
  
  # Calculate event rates for each group
  results$Event_Rate <- results$Events / results$Total
  
  # Verify median survival - if event rate > 0.5 but median is NA, there's likely an issue
  if ("Median_Survival" %in% colnames(results)) {
    # For each group, check if we have > 50% events but NA median
    for (i in 1:nrow(results)) {
      if (is.na(results$Median_Survival[i]) && results$Event_Rate[i] > 0.5) {
        # We should be able to calculate median survival when >50% of subjects have events
        message(sprintf("Group %s has > 50%% events (%.1f%%) but no median survival calculated. Attempting to calculate it now.", 
                        results$Group[i], results$Event_Rate[i] * 100))
        
        # Try to calculate the median for this group
        group_data <- df[df[[treatment_column]] == results$Group[i], ]
        if (nrow(group_data) > 0) {
          # Create a separate survfit object for just this group
          group_surv_formula <- stats::as.formula(paste("Surv(", time_column, ",", censor_column, ") ~ 1"))
          group_km_fit <- survival::survfit(group_surv_formula, data = group_data)
          
          # Extract median (at 0.5)
          if (!is.null(group_km_fit$median)) {
            med_surv <- group_km_fit$median
            if (!is.na(med_surv) && med_surv > 0) {
              results$Median_Survival[i] <- med_surv
              message(sprintf("Successfully calculated median survival for group %s: %.1f days", 
                              results$Group[i], results$Median_Survival[i]))
            } else {
              message(sprintf("Could not calculate valid median survival for group %s despite >50%% events.", 
                              results$Group[i]))
            }
          } else {
            # Try alternate approach using quantiles
            group_quantiles <- summary(group_km_fit)$quantile
            if (!is.null(group_quantiles) && "50%" %in% colnames(group_quantiles)) {
              results$Median_Survival[i] <- group_quantiles["50%"]
              message(sprintf("Calculated median survival using quantiles for group %s: %.1f days", 
                              results$Group[i], results$Median_Survival[i]))
            } else {
              message(sprintf("Could not extract median or quantiles for group %s.", results$Group[i]))
            }
          }
        } else {
          message(sprintf("No data available for group %s to calculate median survival.", results$Group[i]))
        }
      }
    }
  }
  
  # Add reference group note
  results$Note <- ifelse(results$Group == reference_group, "Reference group", "")
  
  # Print formatted results
  if (isTRUE(verbose)) print_results(results, df, treatment_column, time_column, censor_column)
  
  # Build our result list
  result_list <- list(
    results = results,
    reference_group = reference_group,
    method_used = method_used,
    survival_data = data.frame(
      Time = df[[time_column]],
      Event = df[[censor_column]],
      Treatment = df[[treatment_column]]
    )
  )
  
  # Add model if it exists
  if (!is.null(model)) {
    result_list$model <- model
  }

  # Add cox.zph proportional-hazards check when available (cox path only)
  if (!is.null(model_results$ph_test)) {
    result_list$ph_test <- model_results$ph_test
  }

  # CODE_REVIEW.md B7.1 / B7.2 -- survival reports hazard ratios, not marginal
  # means, so its table legitimately differs from the treatment-effects schema.
  # `meta` says so explicitly rather than leaving a consumer to discover it.
  result_list$meta <- me_result_meta(
    analysis_type     = "Survival analysis (Cox PH / Firth / log-rank)",
    model_type_used   = method_used,
    inference         = "frequentist",
    interval_type     = "confidence",
    transform_used    = "none",
    estimate_scale    = "hazard ratio",
    comparison_family = "vs_reference",
    p_adjust_method   = p_adjust_method,
    extra = list(interval_columns_override = c(lower = "CI_Lower",
                                               upper = "CI_Upper"))
  )

  # Cage structure, the declared unit and, under cage randomisation, the
  # permutation record (R20.6). cage_cluster_used stays for older readers; it is
  # always FALSE now.
  result_list$cage_structure     <- cage_structure
  result_list$randomisation_unit <- randomisation_unit
  result_list$cage_permutation   <- cage_permutation
  result_list$cage_cluster_used  <- FALSE
  # Animals randomised individually but housed by arm: the Cox standard errors
  # treat cage-mates as independent. In simulation (2-4 cages per arm, cage
  # frailty SD 0.8, HR = 1) that gave 16-19 % false positives at alpha 0.05, and
  # 2-4 % with no cage effect. Say so; do not guess which applies.
  result_list$cage_caveat <- if (randomisation_unit == "mouse" &&
                                 cage_structure$structure == "nested_replicated") {
    paste0("Each cage holds one treatment. These p-values treat cage-mates as ",
           "independent animals; if cage-mates are more alike than other animals, ",
           "the p-values are too small. If whole cages were assigned to ",
           "treatments, declare the unit of randomisation as cage.")
  } else NULL

  # Add concordance / C-index when available (cox path only)
  if (!is.null(model_results$c_index)) {
    result_list$c_index <- model_results$c_index
  }

  return(result_list)
}

#' Cage-level permutation log-rank: one treated arm against the reference
#'
#' CODE_REVIEW.md R20.6. Whole cages are reassigned between the two arms and the
#' log-rank chi-square recomputed for every assignment (exactly when there are at
#' most \code{max_exact}, otherwise from \code{n_perm} random assignments). The
#' p-value is the share of assignments at least as extreme as the one observed.
#' With c cages per arm there are only choose(2c, c) assignments, so the smallest
#' attainable p is reported alongside: 1/3 for 2 cages per arm, 0.1 for 3.
#'
#' @param df Internal survival frame with Time, Event, Treatment and Cage.
#' @param group,reference_group Arms to compare.
#' @return A list: p_value, min_attainable_p, n_assignments, exact, cages_trt,
#'   cages_ref.
#' @noRd
#' @keywords internal
cage_permutation_logrank <- function(df, group, reference_group,
                                     max_exact = 5000L, n_perm = 4999L,
                                     seed = 20260930L) {
  pair  <- df[df$Treatment %in% c(reference_group, group), , drop = FALSE]
  cages <- unique(pair[, c("Cage", "Treatment")])
  trt_cages <- cages$Cage[cages$Treatment == group]
  n_trt <- length(trt_cages)
  n_ref <- sum(cages$Treatment == reference_group)
  cage_ids <- cages$Cage

  stat_for <- function(treated) {
    arm <- ifelse(pair$Cage %in% treated, "T", "R")
    if (length(unique(arm)) < 2L) return(NA_real_)
    fit <- tryCatch(
      survival::survdiff(survival::Surv(Time, Event) ~ arm,
                         data = data.frame(Time = pair$Time, Event = pair$Event,
                                           arm = arm)),
      error = function(e) NULL)
    if (is.null(fit)) NA_real_ else fit$chisq
  }

  observed <- stat_for(trt_cages)
  n_assign <- choose(n_ref + n_trt, n_trt)
  exact <- n_assign <= max_exact
  stats <- if (exact) {
    vapply(utils::combn(cage_ids, n_trt, simplify = FALSE), stat_for, numeric(1))
  } else {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
      get(".Random.seed", envir = .GlobalEnv)
    } else NULL
    on.exit(if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = .GlobalEnv),
            add = TRUE)
    set.seed(seed)
    vapply(seq_len(n_perm), function(i) stat_for(sample(cage_ids, n_trt)), numeric(1))
  }
  stats <- stats[is.finite(stats)]
  tol <- 1e-9
  if (!is.finite(observed) || !length(stats)) {
    return(list(p_value = NA_real_, min_attainable_p = NA_real_,
                n_assignments = n_assign, exact = exact,
                cages_trt = n_trt, cages_ref = n_ref))
  }
  if (exact) {
    p_value <- mean(stats >= observed - tol)
    min_p   <- mean(stats >= max(stats) - tol)
  } else {
    p_value <- (1 + sum(stats >= observed - tol)) / (length(stats) + 1)
    min_p   <- 1 / (length(stats) + 1)
  }
  list(p_value = p_value, min_attainable_p = min_p, n_assignments = n_assign,
       exact = exact, cages_trt = n_trt, cages_ref = n_ref)
}

#' Validate the one-row-per-subject precondition
#'
#' CODE_REVIEW.md R3.14 — `survival_statistics()` builds its `Surv` object from
#' `df` directly, so every row is treated as an independent subject at risk.
#' Passed a longitudinal frame (one row per measurement occasion), a mouse
#' measured ten times contributes nine censored pseudo-subjects plus one event:
#' risk sets, hazard ratios, log-rank statistics and the KM curve are all wrong
#' and n is inflated by the number of timepoints. The contract was previously
#' neither documented nor checked. Fail loudly instead.
#'
#' @noRd
#' @keywords internal
validate_one_row_per_subject <- function(df, id_column, treatment_column,
                                         cage_column,
                                         caller = "survival_statistics()") {
  if (is.null(id_column) || !id_column %in% colnames(df)) return(invisible(NULL))

  key_parts <- list(as.character(df[[id_column]]))
  if (!is.null(treatment_column) && treatment_column %in% colnames(df)) {
    key_parts <- c(key_parts, list(as.character(df[[treatment_column]])))
  }
  if (!is.null(cage_column) && cage_column %in% colnames(df)) {
    key_parts <- c(key_parts, list(as.character(df[[cage_column]])))
  }
  keys <- do.call(make_mouse_key, key_parts)

  dup_n <- sum(duplicated(keys))
  if (dup_n > 0L) {
    stop(
      caller, " requires one row per animal, but ", dup_n,
      " duplicate subject key(s) were found (", length(unique(keys)),
      " unique animals across ", nrow(df), " rows).\n",
      "This looks like a longitudinal data frame. Fitting a survival model to it ",
      "would treat every measurement occasion as an independent subject.\n",
      "Reduce to one row per animal first — e.g. each animal's last ",
      "observation, carrying its event indicator — then call this function.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' Validate Required Inputs
#' @noRd
validate_inputs <- function(df, time_column, censor_column, treatment_column) {
  required_cols <- c(time_column, censor_column, treatment_column)
  missing_cols <- required_cols[!required_cols %in% colnames(df)]
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
}

#' Check Cage Distribution
#' @noRd
check_cage_distribution <- function(df, treatment_column, cage_column) {
  # Only check if a cage column was supplied AND exists. `NULL %in% colnames(df)`
  # is logical(0), which makes a bare `if` error rather than skip (R15.1).
  if (!is.null(cage_column) && cage_column %in% colnames(df)) {
    cage_treatment_table <- table(df[[cage_column]], df[[treatment_column]])
    message("Cage distribution across treatment groups:")
    message(paste(utils::capture.output(print(cage_treatment_table)), collapse = "\n"))
    
    # Check for collinearity between cage and treatment
    cage_treatment_df <- data.frame(
      Cage = df[[cage_column]],
      Treatment = df[[treatment_column]]
    )
    cage_treatment_counts <- table(cage_treatment_df$Cage, cage_treatment_df$Treatment)
    cage_has_multiple_treatments <- rowSums(cage_treatment_counts > 0) > 1
    if (!any(cage_has_multiple_treatments)) {
      message("Detected collinearity between Cage and Treatment. Using Treatment only in the model.")
    }
  }
}

#' Check for Separation Issues in Survival Data
#' @noRd
check_separation <- function(df, treatment_column, censor_column) {
  # Check for complete separation (groups with all or no events)
  event_by_treatment <- tapply(df[[censor_column]], df[[treatment_column]], function(x) {
    c(sum(x), length(x), sum(x) / length(x))
  })
  
  groups_no_events <- names(event_by_treatment)[sapply(event_by_treatment, function(x) x[1] == 0)]
  groups_all_events <- names(event_by_treatment)[sapply(event_by_treatment, function(x) x[1] == x[2])]
  
  # CODE_REVIEW.md R3.27 — only a group with *zero* events causes separation in
  # a Cox model. A group in which every animal has an event is perfectly
  # estimable; routing it to Firth was contradicted by this function's own
  # message below ("this is not a problem for Cox models") and needlessly
  # replaced an exact partial likelihood with a penalised approximation.
  has_separation <- length(groups_no_events) > 0
  
  if (has_separation) {
    message("Warning: Some groups have perfect separation (no events). This may affect hazard ratio estimates.")
    if (length(groups_no_events) > 0) {
      message("Groups with no events: ", paste(groups_no_events, collapse = ", "))
    }
    if (length(groups_all_events) > 0) {
      message("Note: Groups with all events: ", paste(groups_all_events, collapse = ", "), 
              " (this is not a problem for Cox models)")
    }
  }
  
  return(list(
    has_separation = has_separation,
    groups_no_events = groups_no_events,
    groups_all_events = groups_all_events
  ))
}

#' Fit Appropriate Survival Model
#' @noRd
fit_survival_model <- function(df, surv_obj, cox_formula, treatment_column, treatment_groups,
                              reference_group, time_column, censor_column, separation_info,
                              firth_correction, verbose = TRUE,
                              p_adjust_method = "bonferroni",
                              cage_column = NULL,
                              permutation_logrank = TRUE) {
  
  # Try standard Cox model first
  cox_model <- tryCatch({
    # Create a factor version of the treatment column with the reference level set explicitly
    df$treatment_factor <- factor(df[[treatment_column]], levels = c(reference_group, setdiff(treatment_groups, reference_group)))

    # CODE_REVIEW.md R3.13 — add cluster(cage) when the caller resolved the
    # design as one where cage-level correlation is estimable. cluster() leaves
    # the point estimates unchanged and replaces the model-based variance with a
    # robust sandwich estimate that accounts for within-cage dependence.
    rhs <- "treatment_factor"
    if (!is.null(cage_column) && cage_column %in% colnames(df)) {
      df$.cage_cluster <- factor(df[[cage_column]])
      rhs <- paste(rhs, "+ cluster(.cage_cluster)")
    }
    new_formula <- stats::as.formula(paste("surv_obj ~", rhs))

    # Fit model with explicit reference level
    survival::coxph(new_formula, data = df)
  }, error = function(e) {
    message("Standard Cox model failed: ", e$message)
    NULL
  })
  
  # Check for potential issues
  has_issues <- is.null(cox_model) || separation_info$has_separation
  
  if (has_issues && firth_correction) {
    # Use Firth's bias-reduced Cox model
    method_used <- "coxphf"
    message("Using Firth's bias-reduced Cox model: Surv(time, status) ~ group")
    
    results <- tryCatch({
      # Create analysis data frame with safer group naming
      analysis_df <- data.frame(
        time = df[[time_column]],
        status = df[[censor_column]],
        group = factor(df[[treatment_column]], levels = c(reference_group, setdiff(treatment_groups, reference_group)))
      )
      
      # Fit model using the simplified data frame
      if (requireNamespace("coxphf", quietly = TRUE)) {
        model <- coxphf::coxphf(
          survival::Surv(time, status) ~ group, 
          data = analysis_df
        )
        
        # Extract results
      coefs <- model$coefficients
      hazard_ratios <- exp(coefs)
      confidence_intervals <- exp(confint(model))
      p_values <- model$prob
      
        # Create results data frame with all treatment groups
        results <- data.frame(
          Group = treatment_groups,
          HR = NA,
          CI_Lower = NA,
          CI_Upper = NA,
          P_Value = NA,
          stringsAsFactors = FALSE
        )
        
        # Set reference group values
        ref_idx <- which(results$Group == reference_group)
        results$HR[ref_idx] <- 1
        results$CI_Lower[ref_idx] <- 1
        results$CI_Upper[ref_idx] <- 1
        
        # Fill in values for non-reference groups
        for (i in seq_along(coefs)) {
          group_name <- levels(analysis_df$group)[i + 1]  # +1 because first level is reference
          if(group_name %in% results$Group) {
            idx <- which(results$Group == group_name)
            results$HR[idx] <- hazard_ratios[i]
            results$CI_Lower[idx] <- confidence_intervals[i, 1]
            results$CI_Upper[idx] <- confidence_intervals[i, 2]
            results$P_Value[idx] <- p_values[i]
          }
        }

        # Fallback: for non-reference groups where coxphf failed to converge
        # (p_value is NA), compute a pairwise log-rank p-value instead.
        non_ref_na <- which(results$Group != reference_group & is.na(results$P_Value))
        for (idx in non_ref_na) {
          grp <- results$Group[idx]
          pair_df <- df[df[[treatment_column]] %in% c(reference_group, grp), , drop = FALSE]
          pair_df[[treatment_column]] <- factor(pair_df[[treatment_column]])
          tryCatch({
            pair_formula <- stats::as.formula(
              paste0("survival::Surv(", time_column, ", ", censor_column, ") ~ ", treatment_column)
            )
            lr <- survival::survdiff(pair_formula, data = pair_df)
            results$P_Value[idx] <- round(1 - stats::pchisq(lr$chisq, df = 1), 4)
            message(sprintf("  Pairwise log-rank p-value used for %s (coxphf did not converge): %.4f",
                            grp, results$P_Value[idx]))
          }, error = function(e) {
            message(sprintf("  Could not compute fallback p-value for %s: %s", grp, e$message))
          })
        }
        
        return(list(
          model = model,
          results = results,
          method_used = "coxphf"
        ))
    } else {
        stop("Package 'coxphf' is required but not available")
      }
    }, error = function(e) {
      message("Firth model failed: ", e$message)
      NULL
    })
    
    if (is.null(results)) {
      method_used <- "logrank"
      model <- NULL
      results <- data.frame()
    } else {
      method_used <- results$method_used
      model <- results$model
      results <- results$results
    }
    
  } else if (has_issues) {
    # Use pairwise log-rank tests as fallback (one per non-reference group).
    # An omnibus logrank p-value must not be assigned to all groups — it is a
    # single test for any difference and is not a valid per-comparison p-value.
    method_used <- "logrank"
    message("One or more groups have zero events. Using pairwise log-rank tests.")

    # Omnibus test (stored separately for reference; not used as per-group p-values)
    surv_diff <- survival::survdiff(cox_formula, data = df)
    if (isTRUE(verbose)) message(paste(utils::capture.output(print(surv_diff)), collapse = "\n"))
    omnibus_p <- 1 - stats::pchisq(surv_diff$chisq, df = length(treatment_groups) - 1L)

    # Create basic results without HRs (not estimable when a group has zero events)
    results <- data.frame(
      Group    = treatment_groups,
      HR       = NA,
      CI_Lower = NA,
      CI_Upper = NA,
      P_Value  = NA,
      stringsAsFactors = FALSE
    )

    ref_idx <- which(results$Group == reference_group)
    results$HR[ref_idx]       <- 1
    results$CI_Lower[ref_idx] <- 1
    results$CI_Upper[ref_idx] <- 1

    # Pairwise log-rank: compare each treatment group vs. reference individually.
    #
    # CODE_REVIEW.md R3.2 — the formula MUST be built from column names here.
    # `cox_formula` is `surv_obj ~ Treatment`, where `surv_obj` is a Surv object
    # of length nrow(df) living in the caller's environment. Passing it with a
    # subsetted `data =` makes model.frame() resolve `surv_obj` from that
    # environment (full length) while `Treatment` comes from the subset, so
    # every call failed with "variable lengths differ" and the old
    # `error = function(e) omnibus_p` handler silently substituted the omnibus
    # p-value for every group — reproducing the exact defect Round 1 2.8 was
    # written to remove. Build a self-contained formula instead, and surface
    # failures as NA + a warning rather than a plausible-looking wrong number.
    pair_formula <- stats::as.formula(
      paste0("survival::Surv(`", time_column, "`, `", censor_column, "`) ~ `",
             treatment_column, "`")
    )
    # CODE_REVIEW.md H.4 — this branch runs precisely when a group has zero
    # events, i.e. when the Cox partial likelihood is not estimable and the
    # asymptotic chi-square approximation to the log-rank statistic is at its
    # worst (few events, small n). A permutation log-rank is exactly valid in
    # that regime: it makes no large-sample appeal and needs no bias
    # correction, because it is not estimating a hazard ratio at all. Prefer it
    # when `coin` is available and fall back to the asymptotic test otherwise.
    #
    # Note the division of labour with Firth: Firth corrects the small-sample
    # bias of the *estimate*; permutation gives a valid *p-value*.
    # coin is a hard Import as of v0.10.0, so this is purely the caller's
    # choice now -- there is no "coin is missing" degradation path left.
    use_perm <- isTRUE(permutation_logrank)
    logrank_flavour <- if (use_perm) {
      "permutation log-rank (coin)"
    } else "asymptotic log-rank (survdiff)"

    for (grp in treatment_groups[treatment_groups != reference_group]) {
      pair_data <- df[df[[treatment_column]] %in% c(reference_group, grp), , drop = FALSE]
      pair_data[[treatment_column]] <- factor(
        pair_data[[treatment_column]], levels = c(reference_group, grp)
      )
      pair_p <- tryCatch({
        if (use_perm) {
          # Prefer the EXACT permutation distribution: it is deterministic, so
          # re-running gives the same p-value. Fall back to a Monte Carlo
          # approximation only when the sample is too large to enumerate, and
          # seed it so that path is reproducible too — an analysis function must
          # not return a different number each time it is called.
          n_pair <- nrow(pair_data)
          lt <- if (n_pair <= 30L) {
            tryCatch(
              coin::logrank_test(pair_formula, data = pair_data,
                                 distribution = "exact"),
              error = function(e) NULL
            )
          } else NULL
          if (is.null(lt)) {
            withr_seed <- 20260729L
            old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
              get(".Random.seed", envir = .GlobalEnv)
            } else NULL
            set.seed(withr_seed)
            lt <- coin::logrank_test(
              pair_formula, data = pair_data,
              distribution = coin::approximate(nresample = 20000L)
            )
            if (!is.null(old_seed)) {
              assign(".Random.seed", old_seed, envir = .GlobalEnv)
            }
          }
          as.numeric(coin::pvalue(lt))
        } else {
          pair_diff <- survival::survdiff(pair_formula, data = pair_data)
          1 - stats::pchisq(pair_diff$chisq, df = 1L)
        }
      }, error = function(e) {
        warning("Pairwise log-rank test failed for group '", grp, "': ",
                conditionMessage(e), ". P-value reported as NA.", call. = FALSE)
        NA_real_
      })
      results$P_Value[results$Group == grp] <- pair_p
    }
    attr(results, "logrank_flavour") <- logrank_flavour

    # Retain the omnibus test as an explicitly-labelled separate quantity so it
    # is available without ever being mistaken for a per-group comparison.
    attr(results, "omnibus_logrank_p") <- omnibus_p

    model <- surv_diff
    
  } else {
    # Use standard Cox model
    method_used <- "cox"
    model <- cox_model

    # Proportional-hazards check via Schoenfeld residuals (cox.zph).
    # Surface global + per-covariate test for downstream display.
    ph_test <- tryCatch(survival::cox.zph(model), error = function(e) NULL)
    if (!is.null(ph_test) && isTRUE(verbose)) {
      global_p <- tryCatch(ph_test$table["GLOBAL", "p"], error = function(e) NA_real_)
      if (is.finite(global_p) && global_p < 0.05) {
        message(sprintf("Proportional-hazards check: global p = %.4f (< 0.05) — PH assumption may be violated.", global_p))
      }
    }

    # Concordance (C-index) — survival analogue of AUC-ROC; values around
    # 0.5 indicate no discrimination, 1.0 perfect discrimination.
    c_index <- tryCatch(
      survival::concordance(model),
      error = function(e) NULL
    )

    # Extract the hazard ratios, CIs, and p-values
    model_summary <- summary(model)
    
    # Create results data frame with all treatment groups
    results <- data.frame(
      Group = treatment_groups,
      HR = NA,
      CI_Lower = NA,
      CI_Upper = NA,
      P_Value = NA,
      stringsAsFactors = FALSE
    )
    
    # Set reference group values
    ref_idx <- which(results$Group == reference_group)
    results$HR[ref_idx] <- 1
    results$CI_Lower[ref_idx] <- 1
    results$CI_Upper[ref_idx] <- 1
    
    # Extract coefficient names which should match "treatment_factorTreatmentName"
    coef_names <- rownames(model_summary$coefficients)
    
    # For non-reference groups, extract HR, CI, and p-value
    for (i in seq_along(coef_names)) {
      # Extract treatment group name from coefficient name
      group_name <- gsub("treatment_factor", "", coef_names[i])
      
      # Find corresponding row in results
      idx <- which(results$Group == group_name)
      
      if (length(idx) > 0) {
        # Use summary(coxph)$conf.int directly — already contains exp(coef),
        # lower .95, and upper .95 at the proper qnorm(0.975) ≈ 1.959964.
        ci_row   <- model_summary$conf.int[i, , drop = TRUE]
        hr       <- unname(ci_row["exp(coef)"])
        ci_lower <- unname(ci_row["lower .95"])
        ci_upper <- unname(ci_row["upper .95"])
        p_value  <- model_summary$coefficients[i, "Pr(>|z|)"]
        
        # Assign values
        results$HR[idx] <- hr
        results$CI_Lower[idx] <- ci_lower
        results$CI_Upper[idx] <- ci_upper
        results$P_Value[idx] <- p_value
      }
    }
  }
  
  return(list(
    model = model,
    results = results,
    method_used = method_used,
    ph_test = if (exists("ph_test", inherits = FALSE)) ph_test else NULL,
    c_index = if (exists("c_index", inherits = FALSE)) c_index else NULL
  ))
}

#' Print Formatted Results
#' @noRd
print_results <- function(results, df = NULL, treatment_column = NULL, time_column = NULL, censor_column = NULL) {
  message("\nSurvival Analysis Results:")
  message("=======================")
  
  # Debug output removed for cleaner presentation
  
  for(i in 1:nrow(results)) {
    message(sprintf("\nGroup: %s", results$Group[i]))
    
    # Safely handle HR values
    hr_na <- is.na(results$HR[i])
    ci_lower_na <- is.na(results$CI_Lower[i]) 
    ci_upper_na <- is.na(results$CI_Upper[i])
    
    hr_text <- if(hr_na || (!hr_na && results$HR[i] == 0)) {
      "Hazard Ratio: Not estimable"
    } else {
      sprintf("Hazard Ratio: %.3f (%.3f-%.3f)", 
              results$HR[i], 
              results$CI_Lower[i], 
              results$CI_Upper[i])
    }
    message(hr_text)
    
    if (!is.na(results$P_Value[i])) {
      message(sprintf("P-value: %.4f", results$P_Value[i]))
    }
    
    if ("Median_Survival" %in% colnames(results)) {
      if (!is.na(results$Median_Survival[i])) {
        message(sprintf("Median Survival: %.1f days", results$Median_Survival[i]))
      } else {
        message("Median Survival: Not reached")
      }
    }
    
    # Ensure event counts are properly displayed 
    if (!is.na(results$Events[i]) && !is.na(results$Total[i])) {
      message(sprintf("Events: %d/%d", results$Events[i], results$Total[i]))
    }
    
    if(!is.na(results$Note[i]) && results$Note[i] != "") {
      message(sprintf("Note: %s", results$Note[i]))
    }
  }
  message("\n")
  
  # Print summary table
  message("Summary Table:")
  formatted_table <- data.frame(
    Group = results$Group,
    "HR (95% CI)" = sapply(1:nrow(results), function(i) {
      if(is.na(results$HR[i]) || (!is.na(results$HR[i]) && results$HR[i] == 0)) {
        "Not estimable"
      } else {
        sprintf("%.2f (%.2f-%.2f)", 
                results$HR[i], 
                results$CI_Lower[i], 
                results$CI_Upper[i])
      }
    }),
    "P-value" = sapply(1:nrow(results), function(i) {
      is_ref <- !is.na(results$Note[i]) && results$Note[i] == "Reference group"
      if (is_ref) {
        "Ref"
      } else if (is.na(results$P_Value[i])) {
        "NC"  # Not Converged
      } else {
        sprintf("%.4f", results$P_Value[i])
      }
    }),
    "Events/Total" = sapply(1:nrow(results), function(i) {
      if (!is.na(results$Events[i]) && !is.na(results$Total[i])) {
        sprintf("%d/%d", results$Events[i], results$Total[i])
      } else {
        "NA/NA"
      }
    }),
    stringsAsFactors = FALSE
  )
  
  # Add median survival to table if available
  if ("Median_Survival" %in% colnames(results)) {
    formatted_table$"Median Survival" <- sapply(seq_len(nrow(results)), function(i) {
      if (is.na(results$Median_Survival[i])) "Not reached"
      else sprintf("%.1f days", results$Median_Survival[i])
    })
  }
  
  message(paste(utils::capture.output(print(formatted_table)), collapse = "\n"))
  
  # Return the formatted table invisibly for further use if needed
  invisible(formatted_table)
}