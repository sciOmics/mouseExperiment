#' Synthetic data for tumor growth analysis
#'
#' A dataset containing simulated measurements from a mouse tumor growth experiment
#'
#' @format A data frame with 120 observations of 6 variables:
#'   \describe{
#'     \item{Date}{Date of measurement}
#'     \item{Group}{Treatment group identifier}
#'     \item{ID}{Mouse ID number}
#'     \item{Length}{Tumor length measurement in mm}
#'     \item{Width}{Tumor width measurement in mm}
#'     \item{Survival_Censor}{Binary indicator for survival analysis (1=event occurred, 0=censored)}
#'   }
#'
#' @source {Simulated data for illustration purposes}
#'
"synthetic_data"