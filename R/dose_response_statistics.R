#' Test for Dose-Response Relationship in Tumor Growth Data
#' 
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#'
#' @description
#' Performs statistical tests to determine if there is a significant dose-response 
#' relationship between drug dose levels and tumor volume, testing both linear and 
#' non-linear relationships.
#'
#' @param df Data frame containing tumor growth and dose data.
#' @param dose_column Column containing dose concentrations. Default: "Dose".
#' @param treatment_column Column containing treatment names. Default: "Treatment".
#' @param volume_column Column storing tumor volume measurements. Default: "Volume".
#' @param day_column Column with number of days since experiment start. Default: "Day".
#' @param id_column Column with individual mouse identifiers. Default: "ID".
#' @param cage_column Optional cage column, part of the animal key (treatment +
#'   ID + cage), so an ID reused across cages or arms is not one animal
#'   (CODE_REVIEW.md T1). NULL or absent means no cage information.
#' @param time_point Day to analyse. Default \code{NULL}: the last day on which
#'   every dose group is evaluable (see \code{\link{evaluable_days}}). A day that
#'   is not evaluable is an error. Before v0.26.0 the default was each animal's
#'   own last observation (CODE_REVIEW.md R20.17).
#' @param control_group_name Name of the control arm in \code{treatment_column}.
#'   It must exist, and its dose must be 0 or missing (a missing dose is read
#'   as 0). \code{NULL} (default): the rows with dose 0 are the control. Before
#'   v0.27.0 a name that did not exist was ignored (CODE_REVIEW.md R20.16).
#' @param treatments The arms that form the dose series, besides the control:
#'   one agent at several doses. \code{NULL} (default) takes every other arm,
#'   and is an error when there is more than one, because the analysis would
#'   then put different agents on one dose axis (R20.16). A dose level held by
#'   more than one of the chosen arms is also an error.
#' @param verbose Logical; if TRUE, prints model summaries and statistics to the console. Default: TRUE.
#'
#' @return A list containing:
#'   \item{dose_effect_test}{Statistical test results for dose-dependency}
#'   \item{trend_test}{Results of trend tests. \code{jonckheere_test} is
#'     two-sided since v0.27.0; its \code{direction} field describes the
#'     trend of the dose-group means without having chosen the test
#'     (CODE_REVIEW.md R20.20).}
#'   \item{linear_model}{Linear regression model}
#'   \item{anova_model}{ANOVA model comparing dose groups}
#'   \item{plots}{List of data visualizations}
#'   \item{summary_table}{Data frame summarizing results for each dose level}
#'   \item{endpoint_day}{The day analysed.}
#'   \item{evaluability}{The evaluable-day record: rule, days, and the days and
#'     dose groups it excluded.}
#'   \item{tgi_table}{TGI per dose group at \code{endpoint_day}, from the
#'     endpoint model shared with \code{analyze_drug_synergy()} and
#'     \code{therapeutic_window_metric()}, with 95 % intervals from draws of
#'     its fixed effects. \code{NULL} when the model cannot be fitted.}
#'   \item{series}{The arms analysed: \code{control} and \code{treatments}.}
#'   \item{statistics}{Among others, the EC50 of the fitted log-logistic curve:
#'     \code{ec50} is the dose giving half of the fitted response range
#'     (\code{drc::ED(model, 50)}), with \code{ec50_ci}, a 95 % interval
#'     computed on log dose so it stays positive. \code{ec50_in_range} is
#'     FALSE when the EC50 lies outside the tested doses, and \code{ec50_note}
#'     gives any reason to read it as descriptive. The lower asymptote is
#'     constrained to be at least 0, and the 5-parameter curve is considered
#'     only with at least 6 dose levels (CODE_REVIEW.md R20.18).}
#'
#' @import drc ggplot2 dplyr stats
#' @importFrom stats coef
#' @importFrom stats lm
#' @export
#'
#' @examples
#' \dontrun{
#' # Load and prepare data
#' data <- read.csv("dose_levels_synthetic_data.csv")
#' data <- calculate_volume(data)
#' data <- calculate_dates(data, start_date = "24-Mar", 
#'                        date_format = "%d-%b", year = 2023)
#'                        
#' # Test for dose-response relationship
#' results <- dose_response_statistics(data, dose_column = "Dose")
#' }
dose_response_statistics <- function(df, 
                                   dose_column = "Dose", 
                                   treatment_column = "Treatment",
                                   volume_column = "Volume", 
                                   day_column = "Day", 
                                   id_column = "ID",
                                   cage_column = NULL,
                                   time_point = NULL,
                                   control_group_name = NULL,
                                   treatments = NULL,
                                   verbose = TRUE) {
  
  # Validate input
  required_columns <- c(dose_column, treatment_column, volume_column, day_column, id_column)
  missing_cols <- required_columns[!required_columns %in% colnames(df)]
  if (length(missing_cols) > 0) {
    stop("Missing required columns in data frame: ", paste(missing_cols, collapse = ", "))
  }

  # CODE_REVIEW.md T1 -- the per-animal reductions grouped on (treatment, dose,
  # ID) and, in the growth-rate step, on (dose, ID) alone, merging animals that
  # share an ear tag across cages or arms. Group on the animal key instead; the
  # user's columns keep their names, because the result's analysis_data is
  # indexed by them.
  has_cage <- !is.null(cage_column) && cage_column %in% colnames(df)
  df <- as.data.frame(df)

  # CODE_REVIEW.md R20.16 -- one agent's dose series and its control. Every
  # arm used to go onto one dose axis, so on the Master demo Drug_B (dose 10)
  # and the Drug_A Mid + Drug_B combination (dose 15) were fitted as doses of
  # Drug_A, and `control_group_name` changed nothing.
  series <- dr_select_series(df, dose_column, treatment_column,
                             control_group_name, treatments)
  df <- series$df
  # R20.14: the key includes the dose, so an ID reused across the dose groups
  # of a single-agent layout (one treatment name at every dose) is not one
  # animal.
  df$.mouse_key <- make_mouse_key(
    as.character(df[[treatment_column]]), format(df[[dose_column]]),
    as.character(df[[id_column]]),
    if (has_cage) as.character(df[[cage_column]]) else "")
  
  # CODE_REVIEW.md R20.17 / R20-K -- the default endpoint was each animal's own
  # last observation, so controls removed at the volume limit contributed
  # capped, earlier volumes (EC50 7.5 against a true 2.2), and a requested day
  # silently dropped any arm with no animals left. The analysis now uses one
  # day on which every dose group is evaluable (>= 50 % and >= 3 of its animals
  # on study): by default the last such day. A requested day that is not
  # evaluable is an error, which includes a control that has thinned out.
  dose_num <- suppressWarnings(as.numeric(df[[dose_column]]))
  vol_num  <- suppressWarnings(as.numeric(df[[volume_column]]))
  day_num  <- suppressWarnings(as.numeric(df[[day_column]]))
  keep <- is.finite(dose_num) & is.finite(vol_num) & is.finite(day_num)
  d_std <- data.frame(MouseKey  = df$.mouse_key[keep],
                      Treatment = format(dose_num[keep], trim = TRUE,
                                         drop0trailing = TRUE),
                      Day       = day_num[keep],
                      Volume    = vol_num[keep],
                      stringsAsFactors = FALSE)
  evaluability <- me_evaluability(d_std)
  eval_day <- me_resolve_eval_day(evaluability, time_point)

  # Prepare data for analysis: the animals measured on the evaluation day.
  analysis_data <- prepare_dose_data(df, dose_column = dose_column, treatment_column = treatment_column,
                                    volume_column = volume_column, day_column = day_column,
                                    id_column = id_column, time_point = eval_day)

  # TGI per dose group at the same day, from the endpoint model shared with
  # synergy and the therapeutic window (R20-K). The control is dose 0, by
  # construction of the series (R20.16).
  ctrl_dose <- "0"
  endpoint_model <- me_endpoint_model(d_std)
  tgi_table <- if (!is.na(ctrl_dose) && !is.null(endpoint_model)) {
    dr_tgi_table(endpoint_model, evaluability, eval_day, ctrl_dose)
  } else NULL
  
  # Generate summary statistics
  summary_stats <- generate_summary_statistics(analysis_data, dose_column = dose_column,
                                               volume_column = volume_column, verbose = verbose)
  
  # Create visualizations
  plots <- create_dose_plots(analysis_data, summary_stats, dose_column = dose_column, volume_column = volume_column)
  
  # Perform statistical analyses
  stats_results <- perform_statistical_analyses(analysis_data, dose_column = dose_column, volume_column = volume_column, 
                                              day_column = day_column, id_column = id_column, original_df = df,
                                              verbose = verbose)
  
  # Generate user-friendly report
  if (isTRUE(verbose)) generate_user_report(stats_results, plots)
  
  # Return all results
  return(list(
    dose_effect_test = list(
      linear_model = stats_results$linear_model,
      linear_pvalue = stats_results$statistics$linear_p_value,
      slope = stats_results$statistics$linear_slope
    ),
    trend_test = list(
      jonckheere_test = stats_results$jt_result,
      linear_trend_pvalue = stats_results$statistics$linear_trend_pvalue,
      quadratic_trend_pvalue = stats_results$statistics$quadratic_trend_pvalue
    ),
    linear_model = stats_results$linear_model,
    anova_model = stats_results$anova_model,
    plots = plots,
    summary_table = summary_stats,
    statistics = stats_results$statistics,
    # The evaluation day, the evaluable-day record (R20-K), and TGI per dose
    # group from the endpoint model with intervals from its draws.
    endpoint_day = eval_day,
    evaluability = evaluability,
    tgi_table    = tgi_table,
    control_dose = ctrl_dose,
    series       = list(control = series$control, treatments = series$treatments),
    endpoint_model = me_endpoint_model_info(endpoint_model),
    # R17.3: the per-observation frame the analysis actually ran on, at the
    # chosen time point. The dashboard needs it to rebuild its scatter, box and
    # trend plots, and had been calling `mouseExperiment::prepare_dose_data()` to
    # re-derive it -- a function that is internal and not exported, so the call
    # threw on every run, the tryCatch swallowed it, and all three plots rendered
    # blank. Returning what was analysed also removes the chance of a
    # re-derivation drifting from it (the R14.4 principle).
    analysis_data = analysis_data
  ))
}

# Removed validate_input function as it's now inlined in the main function

#' One agent's dose series and its control (CODE_REVIEW.md R20.16)
#'
#' @return A list: `df` (the rows of the series and the control, with the dose
#'   column numeric and a missing control dose set to 0), `control` (the
#'   control arm name(s)) and `treatments`.
#' @noRd
#' @keywords internal
dr_select_series <- function(df, dose_column, treatment_column,
                             control_group_name = NULL, treatments = NULL) {
  tx   <- as.character(df[[treatment_column]])
  dose <- suppressWarnings(as.numeric(df[[dose_column]]))
  if (!is.null(control_group_name)) {
    if (!control_group_name %in% tx) {
      stop("Control group '", control_group_name, "' is not in column '",
           treatment_column, "'.", call. = FALSE)
    }
    is_ctrl <- !is.na(tx) & tx == control_group_name
    cd <- unique(dose[is_ctrl])
    if (any(is.finite(cd) & cd != 0)) {
      stop("Control group '", control_group_name, "' has dose ",
           paste(cd[is.finite(cd) & cd != 0], collapse = ", "),
           "; the control must have dose 0.", call. = FALSE)
    }
    # A vehicle arm coded with a missing dose used to be dropped silently.
    dose[is_ctrl] <- 0
  } else {
    is_ctrl <- is.finite(dose) & dose == 0
    if (!any(is_ctrl)) {
      stop("No rows with dose 0 to serve as the control. Name the control arm ",
           "with control_group_name.", call. = FALSE)
    }
  }
  control <- sort(unique(tx[is_ctrl]))
  if (length(control) > 1L) {
    stop("Dose 0 holds more than one arm (", paste(control, collapse = ", "),
         "). Name the control with control_group_name.", call. = FALSE)
  }

  others <- sort(unique(tx[!is_ctrl & !is.na(tx)]))
  if (is.null(treatments)) {
    if (length(others) > 1L) {
      stop("The data hold ", length(others), " arms besides the control (",
           paste(others, collapse = ", "), "). A dose-response curve needs one ",
           "agent's dose series; choose its arms with `treatments =`, so that ",
           "different agents are not put on one dose axis.", call. = FALSE)
    }
    treatments <- others
  } else {
    treatments <- unique(as.character(treatments))
    unknown <- setdiff(treatments, tx)
    if (length(unknown)) {
      stop("Treatment(s) not in column '", treatment_column, "': ",
           paste(unknown, collapse = ", "), ".", call. = FALSE)
    }
  }
  # The control is a set of rows, not a name: in a single-agent layout the
  # dose-0 rows can carry the agent's own name.
  in_series <- !is_ctrl & !is.na(tx) & tx %in% treatments
  if (!any(in_series)) {
    stop("No treated arms besides the control.", call. = FALSE)
  }
  no_dose <- in_series & !is.finite(dose)
  if (any(no_dose)) {
    warning(sum(no_dose), " row(s) of the dose series have no numeric dose and ",
            "are left out.", call. = FALSE)
  }
  in_series <- in_series & is.finite(dose)
  arms_at <- tapply(tx[in_series], dose[in_series], function(x) sort(unique(x)),
                    simplify = FALSE)
  shared <- arms_at[vapply(arms_at, length, integer(1L)) > 1L]
  if (length(shared)) {
    stop(paste(sprintf("Dose %s holds more than one arm (%s)", names(shared),
                       vapply(shared, paste, character(1L), collapse = ", ")),
               collapse = "; "),
         ". A dose-response curve needs one agent's dose series; choose its ",
         "arms with `treatments =`.", call. = FALSE)
  }
  if (any(dose[in_series] == 0)) {
    stop("An arm of the dose series has dose 0, which is the control's dose.",
         call. = FALSE)
  }

  out <- df[is_ctrl | in_series, , drop = FALSE]
  out[[dose_column]] <- dose[is_ctrl | in_series]
  list(df = out, control = control,
       treatments = sort(unique(tx[in_series])))
}

#' TGI per dose group at the evaluation day, from the endpoint model
#' @noRd
#' @keywords internal
dr_tgi_table <- function(em, ev, t, ctrl, n_draws = 2000L, seed = 20260930L) {
  arms <- ev$arms
  arms <- arms[order(suppressWarnings(as.numeric(arms)))]
  lm_ <- me_endpoint_logmeans(em, arms, t)
  dr_ <- me_draw_logmeans(em, me_beta_draws(em, n_draws, seed), arms, t)
  on_study <- ev$table[ev$table$Day == t, ]
  do.call(rbind, lapply(seq_along(arms), function(i) {
    a <- arms[i]
    tgi_draws <- 100 * (1 - exp(dr_[, a] - dr_[, ctrl]))
    q <- if (identical(a, ctrl)) c(0, 0) else
      stats::quantile(tgi_draws, c(0.025, 0.975), names = FALSE)
    data.frame(
      Dose        = as.numeric(a),
      Mean_Volume = exp(lm_$log_mean[i]),
      TGI         = if (identical(a, ctrl)) 0 else
        100 * (1 - exp(lm_$log_mean[i] - lm_$log_mean[lm_$Treatment == ctrl])),
      TGI_Lower   = q[1],
      TGI_Upper   = q[2],
      N_Enrolled  = on_study$N_Enrolled[on_study$Treatment == a],
      N_On_Study  = on_study$N_On_Study[on_study$Treatment == a],
      Control     = identical(a, ctrl),
      stringsAsFactors = FALSE)
  }))
}

#' Prepare data for dose-response analysis
#' 
#' @param df Data frame
#' @param dose_column Dose column name
#' @param treatment_column Treatment column name
#' @param volume_column Volume column name
#' @param day_column Day column name
#' @param id_column ID column name
#' @param time_point Specific time point to analyze
#' 
#' @return Prepared data frame
#' @noRd
#' @keywords internal
prepare_dose_data <- function(df, dose_column = "Dose", treatment_column = "Treatment", 
                             volume_column = "Volume", day_column = "Day", 
                             id_column = "ID", time_point = NULL) {
  # Create working copy
  analysis_data <- df
  # dose_response_statistics() adds the animal key; build it here when a caller
  # reaches this helper directly.
  if (!".mouse_key" %in% names(analysis_data)) {
    analysis_data$.mouse_key <- make_mouse_key(
      as.character(analysis_data[[treatment_column]]),
      as.character(analysis_data[[id_column]]))
  }
  
  # Ensure dose is numeric
  analysis_data[[dose_column]] <- as.numeric(analysis_data[[dose_column]])
  # CODE_REVIEW.md R20.14 -- only measured volumes: a row with a missing
  # volume on an animal's last day used to drop the animal, and one on the
  # evaluation day entered the group counts.
  analysis_data <- analysis_data[
    is.finite(suppressWarnings(as.numeric(analysis_data[[volume_column]]))) &
      is.finite(analysis_data[[dose_column]]), , drop = FALSE]

  # Filter to specific time point or use each animal's last measured volume
  if (!is.null(time_point)) {
    analysis_data <- analysis_data[analysis_data[[day_column]] == time_point, ]
    if (nrow(analysis_data) == 0) {
      stop(paste("No data found for time point", time_point))
    }
  } else {
    # Get last measurement for each mouse
    analysis_data <- analysis_data %>%
      dplyr::group_by(.data[[".mouse_key"]]) %>%
      dplyr::filter(.data[[day_column]] == max(.data[[day_column]])) %>%
      dplyr::ungroup()
  }
  
  return(analysis_data)
}

#' Generate summary statistics for each dose level
#'
#' @param analysis_data Prepared data frame
#' @param dose_column Dose column name
#' @param volume_column Volume column name
#' @param verbose Logical; whether to print summary statistics
#'
#' @return Summary statistics table
#' @noRd
#' @keywords internal
generate_summary_statistics <- function(analysis_data, dose_column = "Dose", volume_column = "Volume",
                                        verbose = TRUE) {
  summary_stats <- analysis_data %>%
    dplyr::group_by(.data[[dose_column]]) %>%
    dplyr::summarize(
      mean_volume = mean(.data[[volume_column]], na.rm = TRUE),
      median_volume = median(.data[[volume_column]], na.rm = TRUE),
      sd_volume = sd(.data[[volume_column]], na.rm = TRUE),
      n = dplyr::n(),
      sem_volume = sd_volume / sqrt(n),
      ci95_lower = mean_volume - qt(0.975, n-1) * sem_volume,
      ci95_upper = mean_volume + qt(0.975, n-1) * sem_volume,
      .groups = "drop"
    )

  if (isTRUE(verbose)) {
    message("Summary statistics by dose level:")
    message(paste(utils::capture.output(print(summary_stats)), collapse = "\n"))
  }

  return(summary_stats)
}

#' Create plots for dose-response visualization
#' 
#' @param analysis_data Prepared data frame
#' @param summary_stats Summary statistics table
#' @param dose_column Dose column name
#' @param volume_column Volume column name
#' 
#' @return List of plot objects
#' @keywords internal
create_dose_plots <- function(analysis_data, summary_stats, dose_column = "Dose", volume_column = "Volume") {
  plots <- list()
  
  # Scatter plot with regression line
  plots$scatter <- ggplot2::ggplot(analysis_data, 
                                 ggplot2::aes(x = .data[[dose_column]], y = .data[[volume_column]])) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, color = "blue") +
    ggplot2::labs(title = "Linear Dose-Response Relationship",
                x = "Dose", y = "Tumor Volume") +
    ggplot2::theme_minimal()
  
  # Box plot by dose level
  plots$boxplot <- ggplot2::ggplot(analysis_data, 
                                 ggplot2::aes(x = factor(.data[[dose_column]]),
                                              y = .data[[volume_column]])) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.2, alpha = 0.6) +
    ggplot2::labs(title = "Tumor Volume by Dose Level",
                x = "Dose", y = "Tumor Volume") +
    ggplot2::theme_minimal()
  
  # Bar plot with error bars
  plots$barplot <- ggplot2::ggplot(summary_stats, 
                                 ggplot2::aes(x = .data[[dose_column]], y = .data[["mean_volume"]])) +
    ggplot2::geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = ci95_lower, ymax = ci95_upper), width = 0.2) +
    ggplot2::labs(title = "Mean Tumor Volume by Dose Level (with 95% CI)",
                x = "Dose", y = "Mean Tumor Volume") +
    ggplot2::theme_minimal()
  
  return(plots)
}

#' Perform statistical analyses for dose-response relationship
#' 
#' @param analysis_data Prepared data frame
#' @param dose_column Dose column name
#' @param volume_column Volume column name
#' @param day_column Day column name
#' @param id_column ID column name
#' @param original_df Original data frame
#' 
#' @return List of statistical analysis results
#' @keywords internal
perform_statistical_analyses <- function(analysis_data, dose_column = "Dose", volume_column = "Volume", 
                                        day_column = "Day", id_column = "ID", original_df = NULL,
                                        verbose = TRUE) {
  statistics <- list()
  
  # 1. Linear regression model
  linear_model <- stats::lm(paste(me_bt(volume_column), "~", me_bt(dose_column)), data = analysis_data)
  linear_summary <- summary(linear_model)
  
  if (isTRUE(verbose)) {
    message("Linear regression model:")
    message(paste(utils::capture.output(print(linear_summary)), collapse = "\n"))
  }
  
  # Store key statistics
  statistics$linear_p_value <- linear_summary$coefficients[2, 4]
  statistics$linear_r_squared <- linear_summary$r.squared
  statistics$linear_slope <- linear_summary$coefficients[2, 1]
  
  # 2. ANOVA test
  anova_model <- stats::aov(as.formula(paste(me_bt(volume_column), "~", paste0("factor(", me_bt(dose_column), ")"))), 
                           data = analysis_data)
  anova_summary <- summary(anova_model)
  
  if (isTRUE(verbose)) {
    message("ANOVA model:")
    message(paste(utils::capture.output(print(anova_summary)), collapse = "\n"))
  }
  
  # Store ANOVA p-value
  if (length(anova_summary) > 0 && nrow(anova_summary[[1]]) > 0) {
    statistics$anova_p_value <- anova_summary[[1]][1, "Pr(>F)"]
  } else {
    statistics$anova_p_value <- NA
  }
  
  # 3. Post-hoc Tukey test
  if (!is.na(statistics$anova_p_value) && statistics$anova_p_value < 0.05) {
    tukey_results <- stats::TukeyHSD(anova_model)
    if (isTRUE(verbose)) {
      message("Tukey HSD test:")
      message(paste(utils::capture.output(print(tukey_results)), collapse = "\n"))
    }
    statistics$tukey_results <- tukey_results
  }
  
  # 4. Try to fit non-linear dose-response models
  statistics <- try_nonlinear_models(
    analysis_data, dose_column, volume_column, statistics, linear_model,
    verbose = verbose
  )

  # 5. Growth rate analysis. R20.14: a failure here no longer ends the whole
  # dose-response analysis.
  statistics <- tryCatch(
    analyze_growth_rate(original_df, analysis_data, dose_column,
                        volume_column, day_column, id_column,
                        statistics, verbose = verbose),
    error = function(e) {
      warning("Growth-rate analysis failed: ", conditionMessage(e),
              call. = FALSE)
      statistics
    })

  # 6. Polynomial trend analysis
  statistics <- analyze_polynomial_trends(
    analysis_data, dose_column, volume_column, statistics,
    verbose = verbose
  )
  
  # 7. Jonckheere-Terpstra trend test
  #
  # CODE_REVIEW.md R3.31 — this was hardwired to NULL under a comment
  # referencing "mentioned issues" that appear nowhere in the source, leaving
  # `clinfun` as a declared dependency the package never loaded, the returned
  # `trend_test$jonckheere_test` field permanently NULL, and two reporting
  # branches unreachable — while the test suite verified the test by calling
  # clinfun directly, so the suite was green on disabled functionality.
  #
  # JT is the canonical permutation test for an ordered alternative across dose
  # groups. It depends only on the *ordering* of the dose levels, never their
  # spacing, so it is immune to the unequal-dose-spacing problem that affects
  # the orthogonal-polynomial decomposition in analyze_polynomial_trends()
  # (Round 2 G.7). It also makes no normality assumption.
  jt_result <- run_jonckheere_test(analysis_data, dose_column, volume_column,
                                   verbose = verbose)


  return(list(
    linear_model = linear_model,
    anova_model = anova_model,
    jt_result = jt_result,
    statistics = statistics
  ))
}

#' Try to fit non-linear dose-response models
#' 
#' @param analysis_data Prepared data frame
#' @param dose_column Dose column name
#' @param volume_column Volume column name
#' @param statistics Statistics list
#' @param linear_model Linear model
#' 
#' @return Updated statistics list
#' @keywords internal
try_nonlinear_models <- function(analysis_data, dose_column = "Dose",
                              volume_column = "Volume",
                              statistics = list(), linear_model = NULL,
                              verbose = TRUE) {
  if (requireNamespace("drc", quietly = TRUE)) {
    tryCatch({
      # Prepare data for drc
      drc_data <- analysis_data[!is.na(analysis_data[[dose_column]]) & 
                              !is.na(analysis_data[[volume_column]]), ]
      
      # Handle zero doses
      if (any(drc_data[[dose_column]] == 0)) {
        min_non_zero_dose <- min(drc_data[[dose_column]][drc_data[[dose_column]] > 0])
        epsilon <- min(min_non_zero_dose/10, 0.01)
        drc_data[[dose_column]] <- ifelse(drc_data[[dose_column]] == 0, 
                                       epsilon, 
                                       drc_data[[dose_column]])
      }
      
      # Fit two log-logistic shapes: 4-parameter (symmetric) and 5-parameter
      # (asymmetric — adds an asymmetry parameter). NEITHER constrains the
      # direction of effect; both fit inhibition or stimulation via the sign
      # of the Slope parameter. The previous variable names "decr"/"incr"
      # and labels "inhibition"/"stimulation" were therefore misleading
      # (CODE_REVIEW.md G.2).
      #
      # CODE_REVIEW.md R20.18 -- a tumour volume cannot be negative, so the
      # lower asymptote is constrained to be >= 0. The 5-parameter curve is
      # considered only with at least 6 dose levels: with fewer, AIC chose it
      # in 15 % of LL.4 truths, and its EC50 interval then went below 0 in
      # 82 % of them. L-BFGS-B does not converge reliably for LL.5, so an LL.5
      # fit whose lower asymptote is negative is rejected instead.
      n_dose_levels <- length(unique(analysis_data[[dose_column]][
        is.finite(analysis_data[[dose_column]])]))
      fml <- as.formula(paste(me_bt(volume_column), "~", me_bt(dose_column)))
      dr_model_4p <- drc::drm(
        fml, data = drc_data,
        fct  = drc::LL.4(names = c("Slope", "Lower Limit", "Upper Limit",
                                   "EC50")),
        lowerl = c(-Inf, 0, -Inf, -Inf)
      )
      dr_model_5p <- if (n_dose_levels >= 6L) tryCatch({
        m5 <- drc::drm(
          fml, data = drc_data,
          fct  = drc::LL.5(names = c("Slope", "Lower Limit", "Upper Limit",
                                     "EC50", "Asymmetry")))
        if (stats::coef(m5)[2] >= 0) m5 else NULL
      }, error = function(e) NULL)

      # Compare models on shape (symmetric vs asymmetric); select lower AIC.
      if (is.null(dr_model_5p) || AIC(dr_model_4p) <= AIC(dr_model_5p)) {
        dr_model   <- dr_model_4p
        model_type <- "symmetric"      # LL.4 = symmetric 4-parameter
      } else {
        dr_model   <- dr_model_5p
        model_type <- "asymmetric"     # LL.5 = asymmetric 5-parameter
      }

      # Derive direction independently from the sign of the Slope parameter
      # in the selected model. Slope > 0 (LL.4/LL.5 convention) implies
      # decreasing response with dose → "inhibition"; Slope < 0 → "stimulation".
      slope_est <- tryCatch(
        unname(stats::coef(dr_model)["Slope:(Intercept)"]),
        error = function(e) NA_real_
      )
      model_direction <- if (is.finite(slope_est)) {
        if (slope_est > 0) "inhibition" else "stimulation"
      } else NA_character_

      # Get model summary
      dr_summary <- summary(dr_model)

      if (isTRUE(verbose)) {
        message("Selected dose-response shape: ", model_type)
        message("Inferred direction (from Slope sign): ", model_direction)
        message("Non-linear dose-response model:")
        message(paste(utils::capture.output(print(dr_summary)), collapse = "\n"))
      }
      
      # Extract parameters.
      #
      # R14.1: these four lookups used the drc default names b/c/d/e, but the
      # LL.4()/LL.5() calls above pass `names = c("Slope", "Lower Limit",
      # "Upper Limit", "EC50")`, so the fitted object carries those instead.
      # Zero of the four matched, and every dose-response run reported EC50,
      # Hill slope and both asymptotes as NA — which the dashboard rendered
      # verbatim. The direction check at the Slope line above was updated when
      # the names were introduced (§G.2); this block was not.
      #
      # The `exp()` was wrong too. In LL.4 the `e` parameter is the EC50 on the
      # natural dose scale; exponentiating it would have returned ~1e13 once the
      # name lookup started working. (LL2.4 is the log-parameterised variant
      # where exp() would be correct.)
      #
      # Read by position rather than by name so a future relabelling cannot
      # silently reintroduce this: LL.4 is always ordered b, c, d, e.
      params  <- stats::coef(dr_model)
      pick <- function(pattern, index) {
        hit <- params[grep(pattern, names(params))]
        if (length(hit) == 1L) unname(hit[1]) else unname(params[index])
      }
      statistics$hill_slope  <- pick("^Slope",       1L)
      statistics$lower_limit <- pick("^Lower Limit", 2L)
      statistics$upper_limit <- pick("^Upper Limit", 3L)

      # CODE_REVIEW.md R20.18 -- the EC50 is the dose giving half the fitted
      # response range, ED(model, 50). It was the `e` parameter, which equals
      # it under LL.4 but not under LL.5 (5.28 against an ED50 of 52.28 in one
      # fit), while the interval described ED50. The interval is now computed
      # on log dose, where the delta method is reasonable and the bounds stay
      # positive; on the dose scale it reached [-394.7, 592.5] on the
      # dashboard demo.
      ed <- tryCatch(drc::ED(dr_model, 50, interval = "delta", display = FALSE),
                     error = function(e) NULL)
      statistics$ec50 <- if (!is.null(ed)) unname(ed[1, "Estimate"]) else
        pick("^EC50", 4L)
      statistics$ec50_ci <- if (!is.null(ed) && is.finite(ed[1, "Std. Error"]) &&
                                statistics$ec50 > 0) {
        se_log <- unname(ed[1, "Std. Error"]) / statistics$ec50
        q <- stats::qt(0.975, stats::df.residual(dr_model))
        c(lower = statistics$ec50 * exp(-q * se_log),
          upper = statistics$ec50 * exp(q * se_log),
          se = unname(ed[1, "Std. Error"]))
      }

      # Read the EC50 as descriptive when it lies outside the tested doses or
      # the curve has a parameter for every dose level (nothing then checks
      # its shape).
      tested <- range(analysis_data[[dose_column]][analysis_data[[dose_column]] > 0])
      n_par  <- length(stats::coef(dr_model))
      statistics$n_dose_levels <- n_dose_levels
      statistics$n_curve_parameters <- n_par
      statistics$ec50_in_range <- is.finite(statistics$ec50) &&
        statistics$ec50 >= tested[1] && statistics$ec50 <= tested[2]
      notes <- c(
        if (!statistics$ec50_in_range)
          sprintf("The EC50 lies outside the tested doses (%s to %s).",
                  format(tested[1]), format(tested[2])),
        if (n_par >= n_dose_levels)
          sprintf(paste("The curve has %d parameters for %d dose levels, so",
                        "nothing checks its shape."), n_par, n_dose_levels))
      statistics$ec50_note <- if (length(notes)) paste(notes, collapse = " ") else ""
      
      # Store model. dr_model_type now reports the *shape* (symmetric or
      # asymmetric); dr_model_direction reports the inferred direction
      # (inhibition / stimulation) derived from the Slope sign.
      statistics$dr_model <- dr_model
      statistics$dr_model_type      <- model_type        # shape
      statistics$dr_model_direction <- model_direction   # direction
      
      # Store per-model information criteria.
      # WARNING: linear_aic and nonlinear_aic are NOT directly comparable.
      # stats::lm and drc::drm use different likelihood parameterisations
      # (drc omits the log(2π) constant and estimates a separate variance
      # parameter), so their AIC values are on different scales. Use them
      # within each model family only; do not use delta-AIC for model selection.
      statistics$linear_aic <- AIC(linear_model)
      statistics$nonlinear_aic <- AIC(dr_model)
      statistics$linear_bic <- BIC(linear_model)
      statistics$nonlinear_bic <- BIC(dr_model)
      statistics$aic_comparison_note <-
        "AIC/BIC values from lm and drc are not on the same scale and cannot be directly compared for model selection."

      # CODE_REVIEW.md DIAGNOSTICS gap (5) — frequentist dose-response had no
      # goodness-of-fit diagnostics. Add lack-of-fit test, residuals plot, and
      # a residual-summary table so the dashboard's DR Diagnostics tab has
      # something to render.
      # With a parameter per dose level the lack-of-fit test has no degrees
      # of freedom (it returned NaN with a warning); it is then not run.
      statistics$dr_lack_of_fit <- if (n_par < n_dose_levels) tryCatch({
        # drc::modelFit returns a data frame with one row per nested model:
        # ANOVA F vs smoother (default), one DF per dose level. A small
        # p-value indicates the parametric model fits worse than a one-mean-
        # per-dose smoother — i.e., lack of fit.
        as.data.frame(drc::modelFit(dr_model))
      }, error = function(e) NULL)

      statistics$dr_residuals_df <- tryCatch({
        data.frame(
          Dose      = drc_data[[dose_column]],
          Observed  = drc_data[[volume_column]],
          Fitted    = stats::fitted(dr_model),
          Residual  = stats::residuals(dr_model),
          stringsAsFactors = FALSE
        )
      }, error = function(e) NULL)

      statistics$dr_residuals_plot <- tryCatch({
        if (!requireNamespace("ggplot2", quietly = TRUE)) return(NULL)
        rd <- statistics$dr_residuals_df
        if (is.null(rd)) return(NULL)
        ggplot2::ggplot(rd, ggplot2::aes(x = .data[["Dose"]],
                                         y = .data[["Residual"]])) +
          ggplot2::geom_point(alpha = 0.6) +
          ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                              colour = "red") +
          ggplot2::geom_smooth(method = "loess", se = FALSE,
                               colour = "steelblue", formula = y ~ x) +
          ggplot2::theme_classic() +
          ggplot2::labs(title = "Dose-response residuals",
                        subtitle = "Systematic curvature → wrong functional form",
                        x = "Dose", y = "Residual")
      }, error = function(e) NULL)

    }, error = function(e) {
      message("Non-linear regression failed: ", e$message)
      message("Continuing with linear analysis only.")
    })
  } else {
    message("Package 'drc' not available. Skipping non-linear regression analysis.")
  }
  
  return(statistics)
}

#' Analyze growth rate vs dose relationship
#' 
#' @param df Original data frame
#' @param analysis_data Prepared data frame
#' @param dose_column Dose column name
#' @param volume_column Volume column name
#' @param day_column Day column name
#' @param id_column ID column name
#' @param statistics Statistics list
#' 
#' @return Updated statistics list
#' @keywords internal
analyze_growth_rate <- function(df, analysis_data, dose_column = "Dose", volume_column = "Volume", 
                               day_column = "Day", id_column = "ID", statistics = list(),
                               verbose = FALSE) {
  # Called directly (not through dose_response_statistics()), the data carry no
  # animal key; fall back to the ID alone.
  if (!".mouse_key" %in% names(df)) {
    df$.mouse_key <- make_mouse_key(as.character(df[[id_column]]))
  }
  # CODE_REVIEW.md R20.14 -- growth rates from measured volumes only. Missing
  # and non-positive volumes were set to half the animal's smallest volume, so
  # the NA rows after an animal's death in a full-grid export became
  # "shrinkage" (true rates +0.083 and +0.038 came out -0.048 and -0.046), and
  # an animal with no positive volume at all (a non-take, common) made lm()
  # fail and took the whole analysis with it. Each animal now needs three
  # measured days with a positive volume, or it is left out.
  v    <- suppressWarnings(as.numeric(df[[volume_column]]))
  day  <- suppressWarnings(as.numeric(df[[day_column]]))
  dose <- suppressWarnings(as.numeric(df[[dose_column]]))
  ok   <- is.finite(v) & v > 0 & is.finite(day) & is.finite(dose)
  if (length(unique(day[ok])) < 2L) return(statistics)
  d <- data.frame(key = df$.mouse_key[ok], dose = dose[ok], day = day[ok],
                  lv = log(v[ok]), stringsAsFactors = FALSE)
  rates <- lapply(split(d, paste(d$dose, d$key, sep = "\r")), function(s) {
    if (length(unique(s$day)) < 3L) return(NULL)
    data.frame(dose = s$dose[1], growth_rate = unname(stats::coef(
      stats::lm(lv ~ day, data = s))[2]))
  })
  growth_rates <- do.call(rbind, rates)
  n_left_out <- length(rates) - if (is.null(growth_rates)) 0L else nrow(growth_rates)
  n_no_positive <- length(setdiff(unique(df$.mouse_key), unique(d$key)))
  statistics$growth_rate_animals_left_out <- n_left_out + n_no_positive

  if (!is.null(growth_rates) && length(unique(growth_rates$dose)) >= 2L) {
    names(growth_rates)[1] <- dose_column
    # Test relationship between dose and growth rate
    growth_model <- stats::lm(paste("growth_rate ~", me_bt(dose_column)), data = growth_rates)
    growth_summary <- summary(growth_model)

    if (isTRUE(verbose)) {
      message("Growth rate vs dose model:")
      message(paste(utils::capture.output(print(growth_summary)), collapse = "\n"))
    }

    statistics$growth_dose_p_value <- growth_summary$coefficients[2, 4]
    statistics$growth_dose_r_squared <- growth_summary$r.squared
    statistics$growth_model <- growth_model
  }

  return(statistics)
}

#' Jonckheere-Terpstra trend test across ordered dose groups
#'
#' Permutation test for the ordered alternative "response changes monotonically
#' with dose". Uses only the *ordering* of the dose levels, so unlike the
#' orthogonal-polynomial decomposition in \code{analyze_polynomial_trends()} it
#' is unaffected by unequal dose spacing (e.g. 0 / 10 / 30 / 100 mg/kg), and it
#' assumes no particular response distribution.
#'
#' The test is two-sided (CODE_REVIEW.md R20.20). It chose its alternative
#' from the direction of the dose-group means, which is a look at the data
#' before the test: on null data P(p < 0.05) was 0.096. The direction is still
#' reported, in \code{direction}, as a description.
#'
#' @param analysis_data Prepared data frame.
#' @param dose_column,volume_column Column names.
#' @param verbose Print progress messages.
#' @return An \code{htest}-like list from \code{clinfun::jonckheere.test()} with
#'   added \code{alternative_used} ("two.sided") and \code{direction} fields,
#'   or \code{NULL} when the test cannot be run (clinfun absent, fewer than
#'   three dose levels, or an error).
#' @noRd
#' @keywords internal
run_jonckheere_test <- function(analysis_data, dose_column = "Dose",
                                volume_column = "Volume", verbose = TRUE) {

  doses   <- as.numeric(analysis_data[[dose_column]])
  volumes <- as.numeric(analysis_data[[volume_column]])
  ok      <- is.finite(doses) & is.finite(volumes)
  doses   <- doses[ok]
  volumes <- volumes[ok]

  # JT tests a trend across ordered groups; it needs at least three of them to
  # say anything a two-group test would not.
  if (length(unique(doses)) < 3L) {
    if (isTRUE(verbose)) {
      message("Jonckheere-Terpstra trend test skipped: needs >= 3 distinct ",
              "dose levels (found ", length(unique(doses)), ").")
    }
    return(NULL)
  }

  # The direction of the dose-group means, reported as a description only:
  # choosing the alternative from it doubled the false-positive rate (R20.20).
  grp_means <- tapply(volumes, doses, mean, na.rm = TRUE)
  grp_means <- grp_means[order(as.numeric(names(grp_means)))]
  slope     <- stats::coef(stats::lm(
    grp_means ~ as.numeric(names(grp_means))
  ))[2]
  direction <- if (is.finite(slope) && slope > 0) "increasing" else "decreasing"
  alternative <- "two.sided"

  jt <- tryCatch(
    clinfun::jonckheere.test(volumes, doses, alternative = alternative),
    error = function(e) {
      warning("Jonckheere-Terpstra test failed: ", conditionMessage(e),
              call. = FALSE)
      NULL
    }
  )
  if (is.null(jt)) return(NULL)

  jt$alternative_used <- alternative
  jt$direction <- direction
  if (isTRUE(verbose)) {
    message("Jonckheere-Terpstra trend test (two-sided; means ", direction,
            " with dose): p = ", format.pval(jt$p.value, digits = 3))
  }
  jt
}

#' Analyze polynomial trends in dose-response
#'
#' @param analysis_data Prepared data frame
#' @param dose_column Dose column name
#' @param volume_column Volume column name
#' @param statistics Statistics list
#'
#' @return Updated statistics list
#' @keywords internal
analyze_polynomial_trends <- function(analysis_data, dose_column = "Dose",
                                 volume_column = "Volume",
                                 statistics = list(), verbose = TRUE) {
  tryCatch({
    # Check if we have enough dose levels
    if (length(unique(analysis_data[[dose_column]])) >= 3) {
      # Create categorical factor for dose. CODE_REVIEW.md G.7: when doses
      # are unequally spaced (the common preclinical pattern 0, 10, 30, 100
      # mg/kg etc.), the default contr.poly() generates orthogonal contrasts
      # for evenly-spaced indices 1, 2, 3, … — NOT for the actual dose
      # scores. Reported linear / quadratic / cubic p-values would then be
      # the orthogonal decomposition on the *index axis*, not on the true
      # dose scale.
      #
      # Pass the numeric dose levels as `scores` so contr.poly() builds the
      # contrast matrix on the actual dose values.
      sorted_doses <- sort(unique(as.numeric(analysis_data[[dose_column]])))
      analysis_data$dose_factor <- factor(
        analysis_data[[dose_column]],
        levels = sorted_doses
      )

      stats::contrasts(analysis_data$dose_factor) <- stats::contr.poly(
        n      = length(sorted_doses),
        scores = sorted_doses
      )
      
      # Fit model with polynomial contrasts
      poly_model <- stats::lm(as.formula(paste(me_bt(volume_column), "~ dose_factor")), data = analysis_data)
      poly_summary <- summary(poly_model)
      poly_anova <- stats::anova(poly_model)
      
      if (isTRUE(verbose)) {
        message("Polynomial contrasts for dose-response trends:")
        message(paste(utils::capture.output(print(summary(poly_model))), collapse = "\n"))
        message(paste(utils::capture.output(print(poly_anova)), collapse = "\n"))
      }
      
      # Extract p-values
      coef_table <- coef(summary(poly_model))
      
      if (nrow(coef_table) >= 2) { # Linear term
        statistics$linear_trend_pvalue <- coef_table[2, 4]
      }
      if (nrow(coef_table) >= 3) { # Quadratic term
        statistics$quadratic_trend_pvalue <- coef_table[3, 4]
      }
      if (nrow(coef_table) >= 4) { # Cubic term
        statistics$cubic_trend_pvalue <- coef_table[4, 4]
      }
      
      # Store overall trend significance
      statistics$overall_trend_pvalue <- poly_anova[1, "Pr(>F)"]
      
      # Store model
      statistics$poly_model <- poly_model
    } else {
      message("Not enough dose levels for polynomial trend analysis (need at least 3).")
    }
  }, error = function(e) {
    message("Polynomial trend analysis failed: ", e$message)
  })
  
  return(statistics)
}

#' Generate user-friendly report of dose-response analysis
#' 
#' @param stats_results Statistical analysis results
#' @param plots List of plots
#' 
#' @return Invisible NULL
#' @keywords internal
generate_user_report <- function(stats_results, plots) {
  statistics <- stats_results$statistics
  jt_result <- stats_results$jt_result
  
  if (length(statistics) > 0) {
    message("\n========== DOSE-RESPONSE RELATIONSHIP ANALYSIS ==========\n")
    
    message("QUESTION: Is there a dose-response relationship?\n")
    
    # Summarize key findings
    dose_effect_detected <- FALSE
    evidence_strength <- "No evidence"
    
    # Check linear regression
    if (!is.null(statistics$linear_p_value) && statistics$linear_p_value < 0.05) {
      dose_effect_detected <- TRUE
      if (statistics$linear_p_value < 0.001) {
        evidence_strength <- "Strong evidence"
      } else if (statistics$linear_p_value < 0.01) {
        evidence_strength <- "Good evidence"
      } else {
        evidence_strength <- "Some evidence"
      }
    }
    
    # Check Jonckheere-Terpstra test
    if (!is.null(jt_result) && !is.null(jt_result$p.value) && jt_result$p.value < 0.05) {
      dose_effect_detected <- TRUE
      if (jt_result$p.value < 0.001) {
        evidence_strength <- "Strong evidence"
      } else if (jt_result$p.value < 0.01) {
        evidence_strength <- "Good evidence"
      } else if (evidence_strength == "No evidence") {
        evidence_strength <- "Some evidence"
      }
    }
    
    # Check linear trend
    if (!is.null(statistics$linear_trend_pvalue) && statistics$linear_trend_pvalue < 0.05) {
      dose_effect_detected <- TRUE
      if (statistics$linear_trend_pvalue < 0.001 && evidence_strength != "Strong evidence") {
        evidence_strength <- "Strong evidence"
      } else if (statistics$linear_trend_pvalue < 0.01 && 
                evidence_strength != "Strong evidence" && 
                evidence_strength != "Good evidence") {
        evidence_strength <- "Good evidence"
      } else if (evidence_strength == "No evidence") {
        evidence_strength <- "Some evidence"
      }
    }
    
    # Report conclusion
    message("CONCLUSION: ", evidence_strength, " of a dose-response relationship.\n")
    
    message("TEST RESULTS:\n")
    
    # Linear dose effect
    message("1. Linear Dose Effect Test:")
    message("   Slope = ", round(statistics$linear_slope, 4))
    message("   R^2 = ", round(statistics$linear_r_squared, 4))
    message("   p-value = ", format.pval(statistics$linear_p_value, digits = 3))
    if (statistics$linear_p_value < 0.05) {
      message("   INTERPRETATION: Significant linear relationship between dose and tumor volume.")
      direction <- ifelse(statistics$linear_slope < 0, "decreasing", "increasing")
      message("   As dose increases, tumor volume tends to be ", direction)
    } else {
      message("   INTERPRETATION: No significant linear relationship detected.")
    }
    
    # Trend tests
    message("\n2. Trend Tests (for ordered dose-response):")
    if (!is.null(jt_result) && !is.null(jt_result$p.value)) {
      message("   Jonckheere-Terpstra test p-value = ", format.pval(jt_result$p.value, digits = 3))
      if (jt_result$p.value < 0.05) {
        message("   INTERPRETATION: Significant monotonic trend across dose levels detected.")
      } else {
        message("   INTERPRETATION: No significant monotonic trend across dose levels.")
      }
    }
    
    if (!is.null(statistics$linear_trend_pvalue)) {
      message("\n   Polynomial trend test results:")
      message("   - Linear trend p-value = ", format.pval(statistics$linear_trend_pvalue, digits = 3))
      
      if (!is.null(statistics$quadratic_trend_pvalue)) {
        message("   - Quadratic trend p-value = ", format.pval(statistics$quadratic_trend_pvalue, digits = 3))
      }
      
      if (!is.null(statistics$cubic_trend_pvalue)) {
        message("   - Cubic trend p-value = ", format.pval(statistics$cubic_trend_pvalue, digits = 3))
      }
      
      if (statistics$linear_trend_pvalue < 0.05) {
        message("   INTERPRETATION: Significant linear trend component detected.")
      }
      
      if (!is.null(statistics$quadratic_trend_pvalue) && statistics$quadratic_trend_pvalue < 0.05) {
        message("   INTERPRETATION: Significant quadratic (curved) component in the dose-response relationship.")
      }
    }
    
    # ANOVA results
    message("\n3. Group Differences (ANOVA):")
    message("   p-value = ", format.pval(statistics$anova_p_value, digits = 3))
    if (!is.na(statistics$anova_p_value) && statistics$anova_p_value < 0.05) {
      message("   INTERPRETATION: Significant differences detected between dose groups.")
    } else {
      message("   INTERPRETATION: No significant differences detected between dose groups.")
    }
    
    # Growth rate analysis
    if (!is.null(statistics$growth_dose_p_value)) {
      message("\n4. Growth Rate Analysis:")
      message("   p-value = ", format.pval(statistics$growth_dose_p_value, digits = 3))
      if (statistics$growth_dose_p_value < 0.05) {
        message("   INTERPRETATION: Dose significantly affects tumor growth rate.")
        direction <- ifelse(coef(statistics$growth_model)[2] < 0, "decreases", "increases")
        message("   Higher doses ", direction, " tumor growth rate.")
      } else {
        message("   INTERPRETATION: No significant effect of dose on tumor growth rate detected.")
      }
    }
    
    message("\n=====================================================")
    
    # Plot guide
    message("\nPlease check the returned plots to visualize the dose-response relationship.")
  }
  
  invisible(NULL)
}