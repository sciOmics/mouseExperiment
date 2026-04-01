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
#' @param time_point Optional specific time point (day) to analyze. Default: NULL (uses last time point).
#' @param control_group_name Name of the control group. Default: "Control".
#'
#' @return A list containing:
#'   \item{dose_effect_test}{Statistical test results for dose-dependency}
#'   \item{trend_test}{Results of trend tests}
#'   \item{linear_model}{Linear regression model}
#'   \item{anova_model}{ANOVA model comparing dose groups}
#'   \item{plots}{List of data visualizations}
#'   \item{summary_table}{Data frame summarizing results for each dose level}
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
                                   time_point = NULL,
                                   control_group_name = "Control",
                                   verbose = TRUE) {
  
  # Validate input
  required_columns <- c(dose_column, treatment_column, volume_column, day_column, id_column)
  missing_cols <- required_columns[!required_columns %in% colnames(df)]
  if (length(missing_cols) > 0) {
    stop("Missing required columns in data frame: ", paste(missing_cols, collapse = ", "))
  }
  
  # Prepare data for analysis
  analysis_data <- prepare_dose_data(df, dose_column = dose_column, treatment_column = treatment_column, 
                                    volume_column = volume_column, day_column = day_column, 
                                    id_column = id_column, time_point = time_point)
  
  # Generate summary statistics
  summary_stats <- generate_summary_statistics(analysis_data, dose_column = dose_column, volume_column = volume_column)
  
  # Create visualizations
  plots <- create_dose_plots(analysis_data, summary_stats, dose_column = dose_column, volume_column = volume_column)
  
  # Perform statistical analyses
  stats_results <- perform_statistical_analyses(analysis_data, dose_column = dose_column, volume_column = volume_column, 
                                              day_column = day_column, id_column = id_column, original_df = df,
                                              verbose = verbose)
  
  # Generate user-friendly report
  if (verbose) generate_user_report(stats_results, plots)
  
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
    statistics = stats_results$statistics
  ))
}

# Removed validate_input function as it's now inlined in the main function

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
#' @export
prepare_dose_data <- function(df, dose_column = "Dose", treatment_column = "Treatment", 
                             volume_column = "Volume", day_column = "Day", 
                             id_column = "ID", time_point = NULL) {
  # Create working copy
  analysis_data <- df
  
  # Ensure dose is numeric
  analysis_data[[dose_column]] <- as.numeric(analysis_data[[dose_column]])
  
  # Filter to specific time point or use last time point
  if (!is.null(time_point)) {
    analysis_data <- analysis_data[analysis_data[[day_column]] == time_point, ]
    if (nrow(analysis_data) == 0) {
      stop(paste("No data found for time point", time_point))
    }
  } else {
    # Get last measurement for each mouse
    analysis_data <- analysis_data %>%
      dplyr::group_by(.data[[treatment_column]], .data[[dose_column]], .data[[id_column]]) %>%
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
#' 
#' @return Summary statistics table
#' @export
generate_summary_statistics <- function(analysis_data, dose_column = "Dose", volume_column = "Volume") {
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
  
  message("Summary statistics by dose level:")
  message(paste(utils::capture.output(print(summary_stats)), collapse = "\n"))
  
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
  linear_model <- stats::lm(paste(volume_column, "~", dose_column), data = analysis_data)
  linear_summary <- summary(linear_model)
  
  if (verbose) {
    message("Linear regression model:")
    message(paste(utils::capture.output(print(linear_summary)), collapse = "\n"))
  }
  
  # Store key statistics
  statistics$linear_p_value <- linear_summary$coefficients[2, 4]
  statistics$linear_r_squared <- linear_summary$r.squared
  statistics$linear_slope <- linear_summary$coefficients[2, 1]
  
  # 2. ANOVA test
  anova_model <- stats::aov(as.formula(paste(volume_column, "~", paste0("factor(", dose_column, ")"))), 
                           data = analysis_data)
  anova_summary <- summary(anova_model)
  
  if (verbose) {
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
    if (verbose) {
      message("Tukey HSD test:")
      message(paste(utils::capture.output(print(tukey_results)), collapse = "\n"))
    }
    statistics$tukey_results <- tukey_results
  }
  
  # 4. Try to fit non-linear dose-response models
  statistics <- try_nonlinear_models(analysis_data, dose_column, volume_column, statistics, linear_model)
  
  # 5. Growth rate analysis
  statistics <- analyze_growth_rate(original_df, analysis_data, dose_column, volume_column, 
                                  day_column, id_column, statistics)
  
  # 6. Polynomial trend analysis
  statistics <- analyze_polynomial_trends(analysis_data, dose_column, volume_column, statistics)
  
  # Jonckheere-Terpstra test not run due to mentioned issues
  jt_result <- NULL
  
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
try_nonlinear_models <- function(analysis_data, dose_column = "Dose", volume_column = "Volume", 
                              statistics = list(), linear_model = NULL) {
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
      
      # Fit models — decreasing (LL.4) vs increasing (LL.5, 5-param log-logistic)
      dr_model_decr <- drc::drm(as.formula(paste(volume_column, "~", dose_column)), 
                              data = drc_data, 
                              fct = drc::LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "EC50")))
      
      dr_model_incr <- drc::drm(as.formula(paste(volume_column, "~", dose_column)), 
                              data = drc_data, 
                              fct = drc::LL.5(names = c("Slope", "Lower Limit", "Upper Limit", "EC50", "Asymmetry")))
      
      # Compare models
      model_aic_decr <- AIC(dr_model_decr)
      model_aic_incr <- AIC(dr_model_incr)
      
      # Select better model
      if (model_aic_decr < model_aic_incr) {
        dr_model <- dr_model_decr
        model_type <- "inhibition"
      } else {
        dr_model <- dr_model_incr
        model_type <- "stimulation"
      }
      
      # Get model summary
      dr_summary <- summary(dr_model)
      
      if (verbose) {
        message("Selected dose-response model type: ", model_type)
        message("Non-linear dose-response model:")
        message(paste(utils::capture.output(print(dr_summary)), collapse = "\n"))
      }
      
      # Extract parameters
      params <- dr_model$coefficients
      statistics$ec50 <- exp(params["e:(Intercept)"])
      statistics$hill_slope <- params["b:(Intercept)"]
      statistics$lower_limit <- params["c:(Intercept)"]
      statistics$upper_limit <- params["d:(Intercept)"]
      
      # Store model
      statistics$dr_model <- dr_model
      statistics$dr_model_type <- model_type
      
      # Compare models
      statistics$linear_aic <- AIC(linear_model)
      statistics$nonlinear_aic <- AIC(dr_model)
      statistics$linear_bic <- BIC(linear_model)
      statistics$nonlinear_bic <- BIC(dr_model)
      
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
                               day_column = "Day", id_column = "ID", statistics = list()) {
  # Only run if we have multiple time points
  if (length(unique(df[[day_column]])) > 1) {
    # Calculate growth rate for each mouse
    growth_rates <- df %>%
      dplyr::group_by(.data[[dose_column]], .data[[id_column]]) %>%
      dplyr::mutate(log_volume = log1p(.data[[volume_column]])) %>%
      dplyr::arrange(.data[[day_column]]) %>%
      dplyr::summarize(
        growth_rate = if(dplyr::n() >= 3) {
          # Linear regression to estimate growth rate
          tmp_df <- data.frame(lv = log_volume, dv = .data[[day_column]])
          model <- stats::lm(lv ~ dv, data = tmp_df)
          coef(model)[2] # Slope coefficient = growth rate
        } else {
          NA
        },
        .groups = "drop"
      ) %>%
      dplyr::filter(!is.na(growth_rate))
    
    if (nrow(growth_rates) > 0) {
      # Test relationship between dose and growth rate
      growth_model <- stats::lm(paste("growth_rate ~", dose_column), data = growth_rates)
      growth_summary <- summary(growth_model)
      
      message("Growth rate vs dose model:")
      message(paste(utils::capture.output(print(growth_summary)), collapse = "\n"))
      
      statistics$growth_dose_p_value <- growth_summary$coefficients[2, 4]
      statistics$growth_dose_r_squared <- growth_summary$r.squared
      statistics$growth_model <- growth_model
    }
  }
  
  return(statistics)
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
analyze_polynomial_trends <- function(analysis_data, dose_column = "Dose", volume_column = "Volume", 
                                 statistics = list()) {
  tryCatch({
    # Check if we have enough dose levels
    if (length(unique(analysis_data[[dose_column]])) >= 3) {
      # Create categorical factor for dose
      analysis_data$dose_factor <- factor(analysis_data[[dose_column]], 
                                       levels = sort(unique(analysis_data[[dose_column]])))
      
      # Set up polynomial contrasts
      stats::contrasts(analysis_data$dose_factor) <- stats::contr.poly(levels(analysis_data$dose_factor))
      
      # Fit model with polynomial contrasts
      poly_model <- stats::lm(as.formula(paste(volume_column, "~ dose_factor")), data = analysis_data)
      poly_summary <- summary(poly_model)
      poly_anova <- stats::anova(poly_model)
      
      if (verbose) {
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