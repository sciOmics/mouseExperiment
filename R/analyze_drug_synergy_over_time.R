#' Analyze Drug Combination Synergy Over Time
#'
#' This function tests for synergistic effects of drug combinations in tumor growth data
#' over multiple time points. It extends the analyze_drug_synergy function by evaluating
#' how synergy metrics change over the course of the experiment, providing insights
#' into when synergy may be strongest.
#'
#' @param df A data frame containing tumor growth data.
#' @param treatment_column A character string specifying the column name for treatment groups. Default is "Treatment".
#' @param volume_column A character string specifying the column name for tumor volume measurements. Default is "Volume".
#' @param time_column A character string specifying the column name for time points. Default is "Day".
#' @param drug_a_name A character string specifying the name of the first single agent treatment group.
#' @param drug_b_name A character string specifying the name of the second single agent treatment group.
#' @param combo_name A character string specifying the name of the combination treatment group.
#' @param control_name A character string specifying the name of the control/vehicle group. Default is "Control".
#' @param min_time_point Optional. A numeric value specifying the minimum time point to include in analysis. Default is NULL.
#' @param max_time_point Optional. A numeric value specifying the maximum time point to include in analysis. Default is NULL.
#' @param id_column Column identifying individual animals (required). Before
#'   v0.25.0 this function had no such argument: each day's analysis looked for a
#'   column literally named "ID", and without one every animal in an arm shared a
#'   key, so every day failed (CODE_REVIEW.md R20.15).
#' @param cage_column Optional cage column, part of the animal key.
#' @param endpoint_method,n_boot,boot_seed As in
#'   \code{\link{analyze_drug_synergy}}. Under the model estimand one endpoint
#'   model is fitted and evaluated on every day.
#' @param verbose Logical. If TRUE, prints detailed results to the console.
#'        Default is TRUE for interactive use; set to FALSE for programmatic/dashboard use.
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{timepoint_results}{A list of results at each time point, each containing the outputs from analyze_drug_synergy}
#'   \item{synergy_summary}{One row per evaluable day: TGIs, the Bliss
#'     expectation and difference, their 95 % intervals (\code{*_Lower},
#'     \code{*_Upper}), p-values of the combination against each agent, and
#'     \code{Bliss_Applies} (FALSE when a single agent did not inhibit growth that
#'     day, so Bliss does not apply; named \code{Evaluable} before v0.26.0).}
#'   \item{evaluability}{The evaluable-day record: only days on which all four
#'     arms have at least 50 % of their enrolled animals, and at least 3, on
#'     study are analysed; the others are listed in \code{excluded} with the
#'     arms that failed (CODE_REVIEW.md R20.1, R20-K).}
#'   \item{peak_bliss_synergy}{The row of \code{synergy_summary} with the largest
#'     Bliss difference among evaluable days; a 0-row data frame when no day is
#'     evaluable (a single agent never inhibited growth).}
#'   \item{drug_a_name, drug_b_name, combo_name}{Names of treatment groups for plotting}
#' }
#'
#' @details
#' This function applies the analyze_drug_synergy approach to each time point in the data,
#' allowing for the assessment of how synergy develops over time. It calculates key metrics
#' including:
#'
#' 1. Bliss Independence effect differences at each time point
#' 3. Statistical significance of combination advantage over monotherapies
#'
#' To visualize the results, use the plot_synergy_trend
#' function with the output from this function.
#'
#' @section Assumptions and Limitations:
#' \strong{Bliss Independence applied to TGI:} Bliss Independence was formulated for the
#' probability of cell death, not for proportional growth inhibition. Applying it to TGI is a
#' common pragmatic choice but carries a ceiling effect: when individual drug TGIs are large
#' (each > 50%), the Bliss expected combined TGI approaches 100%, making it nearly impossible
#' to demonstrate synergy by this criterion regardless of the true biological interaction.
#' Interpret Bliss results cautiously when individual-agent TGIs exceed 50% at a given time point.
#'
#' \strong{No Combination Index.} The Loewe / CI path was removed in v0.21.0:
#' what it computed was response additivity rather than Loewe additivity, and it
#' labelled Bliss-additive combinations antagonistic across most of the range
#' where active single agents sit. See \code{\link{analyze_drug_synergy}}.
#'
#' @examples
#' # Analyze synergy over all available time points
#' data(combo_treatment_synthetic_data)
#' data_processed <- calculate_volume(combo_treatment_synthetic_data)
#' data_processed <- calculate_dates(data_processed, start_date = "03/24/2025")
#' 
#' synergy_results <- analyze_drug_synergy_over_time(
#'   df = data_processed,
#'   drug_a_name = "Drug A",
#'   drug_b_name = "Drug B", 
#'   combo_name = "Combo",
#'   control_name = "Control"
#' )
#' 
#' # Access synergy metrics at each time point
#' print(head(synergy_results$synergy_summary))
#' trend_plot <- plot_synergy_trend(synergy_results)
#' print(trend_plot)
#' 
#' # Check when synergy was strongest
#' print(synergy_results$peak_bliss_synergy)
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom ggpubr ggarrange
#' @export
analyze_drug_synergy_over_time <- function(df, 
                                      treatment_column = "Treatment",
                                      volume_column = "Volume",
                                      time_column = "Day",
                                      drug_a_name,
                                      drug_b_name,
                                      combo_name,
                                      control_name = "Control",
                                      min_time_point = NULL,
                                      max_time_point = NULL,
                                      id_column = "ID",
                                      cage_column = NULL,
                                      endpoint_method = c("model", "last_obs", "survivors"),
                                      n_boot = 2000L,
                                      boot_seed = NULL,
                                      verbose = TRUE) {

  endpoint_method <- match.arg(endpoint_method)

  # Input validation
  required_columns <- c(treatment_column, volume_column, time_column, id_column)
  missing_cols <- required_columns[!required_columns %in% colnames(df)]
  
  if (length(missing_cols) > 0) {
    stop("Missing required columns in the data frame: ", paste(missing_cols, collapse = ", "))
  }
  
  # Check that the specified groups exist in the data
  all_groups <- c(drug_a_name, drug_b_name, combo_name, control_name)
  missing_groups <- all_groups[!all_groups %in% unique(df[[treatment_column]])]
  
  if (length(missing_groups) > 0) {
    stop("The following specified groups do not exist in the treatment column: ", 
         paste(missing_groups, collapse = ", "))
  }
  
  if (!is.null(cage_column) && !cage_column %in% colnames(df)) cage_column <- NULL
  arms <- c(control_name, drug_a_name, drug_b_name, combo_name)

  # CODE_REVIEW.md R20.1 / R20-K -- only evaluable days: every one of the four
  # arms has >= 50 % and >= 3 of its animals on study. Days after the control
  # arm has thinned out used to be evaluated anyway, by extrapolation.
  d_std <- data.frame(
    MouseKey  = if (!is.null(cage_column)) {
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
    stringsAsFactors = FALSE)
  d_std <- d_std[is.finite(d_std$Day) & is.finite(d_std$Volume), , drop = FALSE]
  evaluability <- me_evaluability(d_std, arms)
  all_timepoints <- evaluability$days

  # Apply time point filters if specified
  if (!is.null(min_time_point)) {
    all_timepoints <- all_timepoints[all_timepoints >= min_time_point]
  }
  if (!is.null(max_time_point)) {
    all_timepoints <- all_timepoints[all_timepoints <= max_time_point]
  }
  if (length(all_timepoints) == 0) {
    stop(if (!length(evaluability$days)) me_no_evaluable_day_msg(evaluability)
         else "No evaluable day lies within min_time_point and max_time_point.",
         call. = FALSE)
  }

  # One endpoint model for every day, so the days are estimates from the same
  # fit rather than from a refit per day.
  endpoint_model <- if (endpoint_method == "model") me_endpoint_model(d_std) else NULL

  if (isTRUE(verbose)) {
    message("Analyzing drug synergy across ", length(all_timepoints),
            " evaluable days...")
  }

  timepoint_results <- list()
  synergy_rows <- list()

  for (tp in all_timepoints) {
    tryCatch({
      ep <- endpoint_volumes(
        df, id_column = id_column, treatment_column = treatment_column,
        time_column = time_column, volume_column = volume_column,
        cage_column = cage_column, endpoint_day = tp,
        endpoint_method = endpoint_method, arms = arms,
        model = endpoint_model)
      synergy_results <- me_synergy_from_endpoint(
        ep, control_name = control_name, drug_a_name = drug_a_name,
        drug_b_name = drug_b_name, combo_name = combo_name,
        n_boot = n_boot, boot_seed = boot_seed, verbose = FALSE,
        warn_not_evaluable = FALSE)
      timepoint_results[[as.character(tp)]] <- synergy_results

      bliss_result <- synergy_results$bliss_independence
      stat_tests   <- synergy_results$statistical_tests
      summary_df   <- synergy_results$summary
      ci           <- synergy_results$synergy_ci
      ci_of <- function(metric, col) {
        if (is.null(ci)) return(NA_real_)
        v <- ci[[col]][ci$Metric == metric]
        if (length(v)) v else NA_real_
      }

      synergy_rows[[length(synergy_rows) + 1]] <- data.frame(
        Time_Point = tp,
        TGI_Drug_A = summary_df$TGI_Percent[summary_df$Treatment == drug_a_name],
        TGI_Drug_B = summary_df$TGI_Percent[summary_df$Treatment == drug_b_name],
        TGI_Combo  = summary_df$TGI_Percent[summary_df$Treatment == combo_name],
        Bliss_Expected_TGI = summary_df$TGI_Percent[summary_df$Treatment == "Bliss Expected"],
        Bliss_Difference = bliss_result$difference * 100,
        # Intervals for the reported estimates (R20.2), in the same units.
        TGI_Drug_A_Lower = ci_of("TGI_A_pct", "CI_Lower"),
        TGI_Drug_A_Upper = ci_of("TGI_A_pct", "CI_Upper"),
        TGI_Drug_B_Lower = ci_of("TGI_B_pct", "CI_Lower"),
        TGI_Drug_B_Upper = ci_of("TGI_B_pct", "CI_Upper"),
        TGI_Combo_Lower  = ci_of("TGI_Combo_pct", "CI_Lower"),
        TGI_Combo_Upper  = ci_of("TGI_Combo_pct", "CI_Upper"),
        Bliss_Difference_Lower = 100 * ci_of("Bliss_Excess_FE", "CI_Lower"),
        Bliss_Difference_Upper = 100 * ci_of("Bliss_Excess_FE", "CI_Upper"),
        P_Value_vs_Drug_A = stat_tests$P_Value[1],
        P_Value_vs_Drug_B = stat_tests$P_Value[2],
        Synergy_Assessment = synergy_results$overall_assessment,
        # FALSE when a single agent did not inhibit growth that day, so the
        # Bliss columns are NA (R20.8).
        Bliss_Applies = isTRUE(synergy_results$bliss_applies),
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      warning(paste("Error analyzing time point", tp, ":", e$message), call. = FALSE)
    })
  }

  # Combine accumulated rows into summary data frame
  synergy_summary <- if (length(synergy_rows) > 0) do.call(rbind, synergy_rows) else data.frame()

  # Check if we have any successful results
  if (nrow(synergy_summary) == 0) {
    stop("Could not calculate synergy for any time points.")
  }
  
  # Order the summary by time point
  synergy_summary <- synergy_summary[order(synergy_summary$Time_Point), ]

  # One warning for all the days on which Bliss does not apply (R14.2, R20.8).
  not_bliss <- synergy_summary$Time_Point[!synergy_summary$Bliss_Applies]
  if (length(not_bliss)) {
    warning("Bliss synergy is not evaluated on ", length(not_bliss), " day(s) (",
            paste(not_bliss, collapse = ", "), "): a single agent did not inhibit ",
            "growth relative to control on those days.", call. = FALSE)
  }
  
  # Find when Bliss synergy was strongest -- over evaluable days only (R20.8).
  # A day on which an agent accelerated growth has no Bliss difference, so it
  # cannot be the peak. With no evaluable day there is no peak: a 0-row frame.
  evaluable_days <- synergy_summary[is.finite(synergy_summary$Bliss_Difference), ,
                                    drop = FALSE]
  peak_bliss_synergy <- evaluable_days[which.max(evaluable_days$Bliss_Difference), ,
                                       drop = FALSE]
  
  # Print summary of findings (only when verbose)
  if (isTRUE(verbose)) {
    message("\n=== Drug Combination Synergy Analysis Over Time ===")
    message("Analysis performed across ", nrow(synergy_summary), " time points from ", 
        min(synergy_summary$Time_Point), " to ", max(synergy_summary$Time_Point), "\n")
    
    message("Peak Synergy Findings:")
    if (nrow(peak_bliss_synergy) > 0L) {
      message("Strongest Bliss Synergy at Day ", peak_bliss_synergy$Time_Point,
              " (Difference = ", round(peak_bliss_synergy$Bliss_Difference, 1), "%)\n")
    } else {
      message("No evaluable day: a single agent did not inhibit growth on any day.\n")
    }
    
    message("Synergy Summary by Time Point:")
    message(paste(utils::capture.output(
      print(synergy_summary[, c("Time_Point", "TGI_Combo", "Bliss_Expected_TGI",
                             "Bliss_Difference", "Synergy_Assessment")])
    ), collapse = "\n"))
  }
  
  # Return comprehensive results
  return(list(
    timepoint_results = timepoint_results,
    synergy_summary = synergy_summary,
    peak_bliss_synergy = peak_bliss_synergy,
    # The rule, the evaluable days and the days it excluded, with the arms
    # that failed it (R20-K).
    evaluability = evaluability,
    endpoint_method = endpoint_method,
    endpoint_model = me_endpoint_model_info(endpoint_model),
    # Add these for plotting functions
    drug_a_name = drug_a_name,
    drug_b_name = drug_b_name,
    combo_name = combo_name
  ))
}

#' Plot Drug Synergy Trend Over Time
#'
#' Creates a line plot visualizing tumor growth inhibition (TGI) trends over time
#' for different treatment groups, highlighting synergy and antagonism regions.
#'
#' @param synergy_results Results object from analyze_drug_synergy_over_time function
#' @param custom_title Optional custom title for the plot
#' @param custom_colors Optional named vector of custom colors for treatment groups
#'
#' @return A ggplot2 object visualizing TGI trends over time
#' @export
#'
#' @examples
#' \dontrun{
#' # First run the analysis
#' results <- analyze_drug_synergy_over_time(
#'   df = tumor_data,
#'   drug_a_name = "Drug A",
#'   drug_b_name = "Drug B", 
#'   combo_name = "Drug A + Drug B"
#' )
#' 
#' # Then create the trend plot
#' plot_synergy_trend(results)
#' }
plot_synergy_trend <- function(synergy_results, custom_title = NULL, custom_colors = NULL) {
  # Validate input
  if (!is.list(synergy_results) || is.null(synergy_results$synergy_summary)) {
    stop("Input must be a valid result object from analyze_drug_synergy_over_time()")
  }
  
  # Extract necessary data
  synergy_summary <- synergy_results$synergy_summary
  drug_a_name <- synergy_results$drug_a_name
  drug_b_name <- synergy_results$drug_b_name
  combo_name <- synergy_results$combo_name
  
  # Set title
  if (is.null(custom_title)) {
    title <- "Tumor Growth Inhibition and Synergy Over Time"
  } else {
    title <- custom_title
  }
  
  # Create a named color vector to ensure correct assignment
  if (is.null(custom_colors)) {
    color_values <- c("red", "blue", "purple", "gray50")
    names(color_values) <- c(drug_a_name, drug_b_name, combo_name, "Bliss Expected")
  } else {
    color_values <- custom_colors
  }
  
  # Annotation coordinates: additive offsets so they work when Time_Point starts at 0
  .tp_lo   <- min(synergy_summary$Time_Point)
  .tp_span <- max(max(synergy_summary$Time_Point) - .tp_lo, 1)
  .tp_x1   <- .tp_lo + .tp_span * 0.05
  .tp_x2   <- .tp_lo + .tp_span * 0.15
  .tp_x3   <- .tp_lo + .tp_span * 0.20
  .tgi_max <- max(synergy_summary$TGI_Combo, na.rm = TRUE)

  # Create the plot
  trend_plot <- ggplot2::ggplot(synergy_summary, ggplot2::aes(x = Time_Point)) +
    # Treatment TGI lines - use named mapping for consistency
    ggplot2::geom_line(ggplot2::aes(y = TGI_Drug_A, color = drug_a_name)) +
    ggplot2::geom_line(ggplot2::aes(y = TGI_Drug_B, color = drug_b_name)) +
    ggplot2::geom_line(ggplot2::aes(y = TGI_Combo, color = combo_name), linewidth = 1.2) +
    ggplot2::geom_line(ggplot2::aes(y = Bliss_Expected_TGI, color = "Bliss Expected"), linetype = "dashed") +
    # Synergy area (when combo effect > bliss expected)
    # Ribbons only where Bliss applies: non-evaluable days have NA expectations
    # and are dropped explicitly (R20.8).
    ggplot2::geom_ribbon(data = subset(synergy_summary, is.finite(Bliss_Expected_TGI) &
                                         TGI_Combo > Bliss_Expected_TGI),
                      ggplot2::aes(ymin = Bliss_Expected_TGI, ymax = TGI_Combo),
                      fill = "lightgreen", alpha = 0.4) +
    # Antagonism area (when combo effect < bliss expected)
    ggplot2::geom_ribbon(data = subset(synergy_summary, is.finite(Bliss_Expected_TGI) &
                                         TGI_Combo < Bliss_Expected_TGI),
                      ggplot2::aes(ymin = TGI_Combo, ymax = Bliss_Expected_TGI),
                      fill = "pink", alpha = 0.4) +
    # Formatting
    ggplot2::scale_color_manual(
      name = "Treatment",
      values = color_values,
      labels = c(drug_a_name, drug_b_name, combo_name, "Bliss Expected")
    ) +
    # Manual legend for ribbon areas
    ggplot2::annotate("rect",
              xmin = .tp_x1, xmax = .tp_x2,
              ymin = .tgi_max * 0.85, ymax = .tgi_max * 0.90,
              fill = "lightgreen", alpha = 0.4) +
    ggplot2::annotate("text",
              x = .tp_x3, y = .tgi_max * 0.875,
              label = "Synergy", hjust = 0) +
    ggplot2::annotate("rect",
              xmin = .tp_x1, xmax = .tp_x2,
              ymin = .tgi_max * 0.75, ymax = .tgi_max * 0.80,
              fill = "pink", alpha = 0.4) +
    ggplot2::annotate("text",
              x = .tp_x3, y = .tgi_max * 0.775,
              label = "Antagonism", hjust = 0) +
    ggplot2::labs(
      title = title,
      subtitle = paste0("Green area = synergy (", combo_name, " > Bliss Expected)\n",
                       "Pink area = antagonism (", combo_name, " < Bliss Expected)"),
      x = "Time Point",
      y = "Tumor Growth Inhibition (%)"
    ) +
    ggplot2::theme_minimal()
  
  return(trend_plot)
}