#' Bayesian Tumor Growth Analysis with Multi-Core Support
#'
#' This function performs a Bayesian analysis of tumor growth using the `brms` package,
#' incorporating multi-core processing for faster model fitting. The function fits a
#' Bayesian model to the tumor growth data, allowing for random and fixed effects,
#' and computes the Bayes factor for model comparison.
#'
#' @param df A data frame containing the tumor growth data.
#' @param time_column A string specifying the column name for the time variable (default: "Day").
#' @param volume_column A string specifying the column name for the tumor volume (default: "Volume").
#' @param group_column A string specifying the column name for the treatment group (default: "Group").
#' @param id_column A string specifying the column name for the unique mouse ID (default: "ID").
#' @param cores An integer specifying the number of CPU cores to use for parallel processing (default: 4).
#'
#' @return A list containing:
#' \item{model}{The fitted Bayesian model object.}
#' \item{posterior_samples}{A data frame containing posterior samples of the model parameters.}
#' \item{pp_check}{Posterior predictive check object (if available).}
#' \item{bayes_factor}{The computed Bayes factor comparing the full model to a time-only null model.}
#' \item{log_bayes_factor}{The natural logarithm of the Bayes factor (useful for very large values).}
#' \item{bayes_factor_comparison}{Detailed results from bridge sampling and model comparison.}
#' \item{null_model}{The fitted null model (time-only) used for Bayes factor calculation.}
#' \item{priors}{A list containing information about the priors used in the full model.}
#' \item{null_model_priors}{A list containing information about the priors used in the null model.}
#' \item{bayes_factor_interpretation}{Text interpretation of the Bayes factor strength of evidence.}
#' \item{plot_pp_check}{A function that safely plots the posterior predictive check when called.}
#'
#' @details 
#' This function uses weakly informative priors that allow the data to drive parameter estimation:
#' \itemize{
#'   \item Fixed effects: Normal(0, 5) prior for all fixed effects (intercept, time, group, interactions)
#'   \item Random effects: Half-Cauchy(0, 2) prior for the standard deviation of random effects
#' }
#' 
#' The Bayes factor calculation compares the full model (with treatment group effects) to a 
#' simpler time-only model to directly test whether group differences have explanatory power
#' beyond the general effect of time on tumor growth.
#'
#' @import brms
#' @import future
#' @import future.apply
#' @import bridgesampling
#'
#' @examples
#' \dontrun{
#' # Load and prepare data
#' data(synthetic_data)
#' df <- calculate_volume(synthetic_data)
#' df <- calculate_dates(df, start_date = "2022-02-24")
#' 
#' # Run Bayesian analysis
#' results <- tumor_growth_statistics_bayes(df)
#' 
#' # Access key components
#' summary(results$model)
#' 
#' # View posterior samples
#' head(results$posterior_samples)
#' 
#' # Check Bayes factor if available
#' if (!is.null(results$bayes_factor)) {
#'   # For very large Bayes factors, use log scale
#'   if (results$bayes_factor > 1e10) {
#'     cat("Log Bayes Factor:", results$log_bayes_factor, "\n")
#'     cat("Bayes Factor: approximately 10^", 
#'         round(results$log_bayes_factor/log(10), 1), "\n", sep="")
#'   } else {
#'     cat("Bayes Factor:", formatC(results$bayes_factor, format="f", digits=2), "\n")
#'   }
#'   
#'   cat("Interpretation:", results$bayes_factor_interpretation, "\n")
#'   
#'   # Compare the null and full models
#'   cat("\nFull model formula:", as.character(results$model$formula)[1], "\n")
#'   cat("Null model formula:", as.character(results$null_model$formula)[1], "\n")
#'   
#'   # View prior information
#'   cat("\nPriors used in full model:\n")
#'   cat("- Fixed effects:", results$priors$fixed_effects, "\n")
#'   cat("- Random effects:", results$priors$random_effects, "\n")
#' }
#' 
#' # Plot posterior predictive check
#' pp_plot <- results$plot_pp_check()
#' if (!is.null(pp_plot)) print(pp_plot)
#'
#' # To see all prior specifications directly from the model object
#' print(results$model$prior)
#' 
#' # Create caterpillar plot of parameter estimates
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   if (exists("plot_caterpillar_bayes")) {
#'     cat_plot <- plot_caterpillar_bayes(results$model)
#'     print(cat_plot)
#'   }
#' }
#' }
#' 
#' @export
tumor_growth_statistics_bayes <- function(df, time_column = "Day", volume_column = "Volume", group_column = "Group", id_column = "ID", cores = 4) {

  # Ensure required columns exist
  if (!all(c(time_column, volume_column, group_column, id_column) %in% colnames(df))) {
    stop("One or more specified columns are not found in the dataframe.")
  }

  # Create a unique identifier for each mouse
  df$Mouse_ID <- factor(paste(df[[group_column]], df[[id_column]], sep = "_"))

  # Specify the Bayesian model formula
  formula <- as.formula(paste(volume_column, "~", time_column, "*", group_column, "+ (1 | Mouse_ID)"))

  # Set up multi-core processing in a platform-independent way
  if (.Platform$OS.type == "windows") {
    future::plan(future::multisession, workers = cores)
  } else {
    future::plan(future::multicore, workers = cores)
  }

  # Define prior distributions
  fixed_effects_prior <- brms::prior(normal(0, 5), class = "b")            # Prior for fixed effects (normal with SD=5)
  random_effects_prior <- brms::prior(cauchy(0, 2), class = "sd", group = "Mouse_ID")  # Half-Cauchy prior for random effects SD
  
  # Store prior information for reference
  prior_info <- list(
    fixed_effects = "Normal(0, 5) for all fixed effects (intercept, time, group, interactions)",
    random_effects = "Half-Cauchy(0, 2) for the standard deviation of the random mouse effects",
    description = paste(
      "The Normal(0, 5) prior for fixed effects is weakly informative, allowing for reasonable effect sizes without",
      "being too constraining. The Half-Cauchy(0, 2) prior for random effects standard deviation is recommended by",
      "Gelman (2006) as a weakly informative prior that allows for both small and large variance components."
    )
  )
  
  # Combine priors
  priors <- c(fixed_effects_prior, random_effects_prior)
  
  # Print prior information
  message("Using the following priors:")
  message("- Fixed effects: ", prior_info$fixed_effects)
  message("- Random effects: ", prior_info$random_effects)
  
  # Fit the Bayesian model using the brms package with multi-core support
  bayesian_model <- brms::brm(formula, data = df, family = gaussian(),
                        prior = priors,
                        chains = cores, iter = 2000, warmup = 1000, thin = 1,
                        control = list(adapt_delta = 0.95), cores = cores,
                        save_pars = brms::save_pars(all = TRUE),   # Save all parameters for bridgesampling
                        silent = 2)

  # Print the model summary
  print(summary(bayesian_model))

  # Posterior predictive checks (PPC) with error handling and safer plotting
  ppcheck <- NULL
  tryCatch({
    message("Generating posterior predictive checks...")
    
    # Create PPC object without trying to plot it immediately
    ppcheck <- brms::pp_check(bayesian_model, type = "dens_overlay", ndraws = 50, plot = FALSE)
    
    # Store the PPC object but don't try to print/plot it here - it can be plotted later
    # This avoids device/file handling errors in notebooks/markdown documents
    
    message("Posterior predictive checks object created successfully.")
  }, error = function(e) {
    warning(paste("Posterior predictive checks failed:", e$message))
  })

  # Compute the Bayes factor for model comparison using the bridgesampling package
  bayes_factor <- NULL
  bayes_factor_result <- NULL
  
  # Try to calculate bridge sampling with error handling and better default model
  tryCatch({
    message("Calculating Bayes factor using bridge sampling...")
    
    # 1. Create a more reasonable null model that includes time but not group or interaction
    # This is a better comparison model than intercept-only for tumor growth data
    null_formula <- as.formula(paste(volume_column, "~", time_column, "+ (1 | Mouse_ID)"))
    
    # Define priors for the null model (consistent with the full model)
    null_intercept_prior <- brms::prior(normal(0, 5), class = "Intercept")  # Prior for intercept
    null_time_prior <- brms::prior(normal(0, 2), class = "b")               # Prior for time effect
    null_random_prior <- brms::prior(cauchy(0, 2), class = "sd", group = "Mouse_ID")  # Prior for random effects
    
    # Combine priors for null model
    null_priors <- c(null_intercept_prior, null_time_prior, null_random_prior)
    
    # Print null model prior information
    message("Null model priors:")
    message("- Intercept: Normal(0, 5)")
    message("- Time effect: Normal(0, 2)")
    message("- Random effects: Half-Cauchy(0, 2)")
    
    # 2. Fit the null model with the same settings as the full model
    null_model <- brms::brm(null_formula, data = df, family = gaussian(),
                          prior = null_priors,
                          chains = 2, iter = 1000, warmup = 500, 
                          control = list(adapt_delta = 0.9),
                          cores = min(2, cores),  # Use fewer cores for null model
                          save_pars = brms::save_pars(all = TRUE),
                          silent = 2)
    
    # 3. Calculate log marginal likelihood for both models
    full_bridge <- bridgesampling::bridge_sampler(bayesian_model, silent = TRUE)
    null_bridge <- bridgesampling::bridge_sampler(null_model, silent = TRUE)
    
    # 4. Compare the models to get Bayes factor
    bf_result <- bridgesampling::bayes_factor(full_bridge, null_bridge)
    
    # 5. Store results
    bayes_factor <- bf_result$bf
    bayes_factor_result <- list(
      full_model_bridge = full_bridge,
      null_model_bridge = null_bridge,
      comparison = bf_result
    )
    
    message("Bridge sampling calculation completed successfully.")
    
    # Calculate the log Bayes factor - better for very large values
    log_bf <- tryCatch({
      bridgesampling::bf(full_bridge, null_bridge, log = TRUE)
    }, error = function(e) {
      warning("Error calculating log Bayes factor: ", e$message)
      return(NA_real_)
    })
    
    # Only proceed with numeric calculations if log_bf is a number
    if (!is.na(log_bf) && is.numeric(log_bf)) {
      # Calculate the regular Bayes factor
      bayes_factor <- exp(log_bf)
      
      # Check for potentially inaccurate extremely high values
      if (bayes_factor > 1e10) {
        message(paste("Log Bayes Factor (full vs. time-only model):", round(log_bf, 2)))
        message(paste("Bayes Factor is approximately 10^", round(log_bf/log(10), 1), sep=""))
        message("NOTE: Very high Bayes factors (>10^10) may be numerically unstable.")
        message("      Consider the log Bayes factor or relative model performance instead.")
      } else {
        message(paste("Bayes Factor (full model vs. time-only model):", formatC(bayes_factor, format = "f", digits = 2)))
      }
    } else {
      warning("Could not calculate log Bayes factor, using direct Bayes factor comparison instead.")
      # If log version failed, try the direct version which is sometimes more stable
      bayes_factor <- tryCatch({
        bf_result$bf
      }, error = function(e) {
        return(NA_real_)
      })
      
      if (!is.na(bayes_factor) && is.numeric(bayes_factor)) {
        message(paste("Bayes Factor (full model vs. time-only model):", formatC(bayes_factor, format = "f", digits = 2)))
      } else {
        warning("Bayes factor calculation failed entirely. Model comparison results may be unreliable.")
        bayes_factor <- NA_real_
      }
    }
    
    # Add clearer interpretation
    message("COMPARISON: Full model with group effect and interaction vs. null model with only time effect")
    
    # Interpret either using log Bayes factor or regular Bayes factor, depending on what's available
    if (!is.na(log_bf) && is.numeric(log_bf)) {
      # Interpret using log Bayes factor
      if (log_bf > log(100)) {
        message("Interpretation: Decisive evidence that treatment groups affect tumor growth rates")
      } else if (log_bf > log(30)) {
        message("Interpretation: Very strong evidence that treatment groups affect tumor growth rates")
      } else if (log_bf > log(10)) {
        message("Interpretation: Strong evidence that treatment groups affect tumor growth rates")
      } else if (log_bf > log(3)) {
        message("Interpretation: Substantial evidence that treatment groups affect tumor growth rates")
      } else if (log_bf > 0) {
        message("Interpretation: Weak evidence that treatment groups affect tumor growth rates")
      } else {
        message("Interpretation: No evidence that treatment groups affect tumor growth differently")
      }
    } else if (!is.na(bayes_factor) && is.numeric(bayes_factor)) {
      # Interpret using direct Bayes factor
      if (bayes_factor > 100) {
        message("Interpretation: Decisive evidence that treatment groups affect tumor growth rates")
      } else if (bayes_factor > 30) {
        message("Interpretation: Very strong evidence that treatment groups affect tumor growth rates")
      } else if (bayes_factor > 10) {
        message("Interpretation: Strong evidence that treatment groups affect tumor growth rates")
      } else if (bayes_factor > 3) {
        message("Interpretation: Substantial evidence that treatment groups affect tumor growth rates")
      } else if (bayes_factor > 1) {
        message("Interpretation: Weak evidence that treatment groups affect tumor growth rates")
      } else if (!is.na(bayes_factor)) {
        message("Interpretation: No evidence that treatment groups affect tumor growth differently")
      }
    } else {
      message("Interpretation: Unable to interpret Bayes factor (calculation failed)")
    }
    
  }, error = function(e) {
    warning(paste("Bridge sampling failed:", e$message, 
                  "\nReturning model without Bayes factor."))
  })

  # Extract posterior samples
  posterior_samples <- as.data.frame(bayesian_model)

  # Prepare results list with all components, using NULL for those that weren't calculated
  results <- list(
    model = bayesian_model,
    posterior_samples = posterior_samples,
    pp_check = ppcheck,
    bayes_factor = bayes_factor,
    bayes_factor_comparison = bayes_factor_result,
    
    # Add an explicit NULL model component for easier access
    null_model = if(exists("null_model")) null_model else NULL,
    
    # Add prior information for reference
    priors = if(exists("prior_info")) prior_info else list(
      fixed_effects = "Normal(0, 5) for all fixed effects",
      random_effects = "Half-Cauchy(0, 2) for random effects SD",
      description = "Weakly informative priors that allow the data to drive parameter estimation"
    ),
    
    # Add null model prior information if available
    null_model_priors = if(exists("null_priors")) list(
      intercept = "Normal(0, 5)",
      time_effect = "Normal(0, 2)",
      random_effects = "Half-Cauchy(0, 2)"
    ) else NULL
  )
  
  # Store and interpret Bayes factor information
  # First add the log_bayes_factor to results if available
  if (exists("log_bf") && !is.na(log_bf) && is.numeric(log_bf)) {
    results$log_bayes_factor <- log_bf
  } else if (!is.null(bayes_factor) && !is.na(bayes_factor) && is.numeric(bayes_factor)) {
    # If we have a valid Bayes factor but not the log, calculate it
    results$log_bayes_factor <- log(bayes_factor)
  } else {
    # If neither is valid, set to NA
    results$log_bayes_factor <- NA_real_
  }
  
  # Add interpretation based on what's available
  if (!is.null(bayes_factor) && !is.na(bayes_factor) && is.numeric(bayes_factor)) {
    # Create interpretation based on the Bayes factor value
    if (bayes_factor > 1e10 && !is.na(results$log_bayes_factor)) {
      # For extremely large values, use scientific notation if log is available
      results$bayes_factor_interpretation <- paste0(
        "Very high Bayes factor (approximately 10^", 
        round(results$log_bayes_factor/log(10), 1), 
        "), suggesting decisive evidence that treatment groups affect tumor growth. ",
        "Note: Very high values may be numerically unstable."
      )
    } else if (bayes_factor > 100) {
      results$bayes_factor_interpretation <- "Decisive evidence that treatment groups affect tumor growth rates differently than the time-only model"
    } else if (bayes_factor > 30) {
      results$bayes_factor_interpretation <- "Very strong evidence that treatment groups affect tumor growth rates"
    } else if (bayes_factor > 10) {
      results$bayes_factor_interpretation <- "Strong evidence that treatment groups affect tumor growth rates"
    } else if (bayes_factor > 3) {
      results$bayes_factor_interpretation <- "Substantial evidence that treatment groups affect tumor growth rates"
    } else if (bayes_factor > 1) {
      results$bayes_factor_interpretation <- "Weak evidence that treatment groups affect tumor growth rates"
    } else {
      results$bayes_factor_interpretation <- "No evidence that treatment groups affect tumor growth differently from time alone"
    }
  } else {
    # If no valid Bayes factor is available
    results$bayes_factor_interpretation <- "Bayes factor could not be calculated or is invalid"
  }
  
  # Add a helper function to safely display the posterior predictive check
  results$plot_pp_check <- function() {
    if (!is.null(results$pp_check)) {
      tryCatch({
        # Create a fresh plotting environment
        plot_obj <- brms::pp_check(results$model, type = "dens_overlay", ndraws = 50)
        return(plot_obj)
      }, error = function(e) {
        message("Error plotting posterior predictive check: ", e$message)
        return(NULL)
      })
    } else {
      message("No posterior predictive check object available")
      return(NULL)
    }
  }
  
  # Return the results as a list
  return(results)
}
