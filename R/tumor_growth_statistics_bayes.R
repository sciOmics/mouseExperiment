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
#' \item{bayes_factor}{The computed Bayes factor for model comparison.}
#'
#' @import brms
#' @import future
#' @import future.apply
#' @import bridgesampling
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

  # Set up multi-core processing (future package)
  plan(multicore, workers = cores)  # Specify number of cores to use

  # Fit the Bayesian model using the brms package with multi-core support
  bayesian_model <- brm(formula, data = df, family = gaussian(),
                        prior = c(
                          prior(normal(0, 5), class = "b"),           # Priors for fixed effects
                          prior(cauchy(0, 2), class = "sd")           # Priors for random effects
                        ),
                        chains = cores, iter = 2000, warmup = 1000, thin = 1,
                        control = list(adapt_delta = 0.95), cores = cores)

  # Print the model summary
  print(summary(bayesian_model))

  # Posterior predictive checks (PPC)
  pp_check(bayesian_model)

  # Compute the Bayes factor for model comparison using the bridgesampling package
  # Explicitly call `bf` from bridgesampling
  bayes_factor_result <- bridgesampling::bridge_sampler(bayesian_model)
  bayes_factor <- bayes_factor_result$bf

  # Extract posterior samples
  posterior_samples <- as.data.frame(bayesian_model)

  # Return the results as a list
  return(list(model = bayesian_model, posterior_samples = posterior_samples, bayes_factor = bayes_factor))
}
