# Copyright (c) 2025 Insight BioAnalytics. All rights reserved.
# Proprietary and confidential.

#' Test dataset
#'
#' A simple dataset for testing the package
#'
#' @format A data frame with 10 rows and 2 columns:
#' \describe{
#'   \item{ID}{Row identifier}
#'   \item{Value}{Test value}
#' }
#' @source Simulated data
#' @export
"my_data"

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
#' @export
"synthetic_data"

#' Treatment schedule for combination treatment experiments
#'
#' A dataset containing the dosing schedule for the combination treatment experiment.
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{Treatment}{Treatment group name}
#'   \item{Day}{Study day when treatment was administered}
#'   \item{Dose}{Dose level administered}
#'   \item{Route}{Route of administration}
#' }
#' @source Simulated data
#' @export
"combo_treatment_schedule"

#' Treatment schedule for dose-response experiments
#'
#' A dataset containing the dosing schedule for the dose-response experiment.
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{Treatment}{Treatment group name}
#'   \item{Day}{Study day when treatment was administered}
#'   \item{Dose}{Dose level administered in mg/kg}
#'   \item{Route}{Route of administration}
#' }
#' @source Simulated data
#' @export
"dose_levels_treatment_schedule"

#' Synthetic data for tumor growth experiments with combination treatment
#'
#' A dataset containing synthetic measurements of tumor growth in mice treated
#' with different combinations of treatments (Control, Drug A, Drug B, and combination).
#'
#' @format A data frame with rows and variables:
#' \describe{
#'   \item{Date}{Date of measurement in MM/DD/YYYY format}
#'   \item{Cage}{Cage number (1-8)}
#'   \item{ID}{Mouse identifier within cage}
#'   \item{Treatment}{Treatment group (Control, Drug A, Drug B, Combo)}
#'   \item{Length}{Tumor length measurement in mm}
#'   \item{Width}{Tumor width measurement in mm}
#'   \item{Survival_Censor}{Censoring indicator (1 = event/death, 0 = censored/alive)}
#' }
#' @source Synthetic data created for demonstration purposes
#' @export
"combo_treatment_synthetic_data"

#' Synthetic data for dose-response experiments
#'
#' A dataset containing synthetic measurements of tumor growth in mice treated
#' with different dose levels of a hypothetical compound.
#'
#' @format A data frame with rows and variables:
#' \describe{
#'   \item{Date}{Date of measurement in DD-MMM format}
#'   \item{Cage}{Cage number}
#'   \item{ID}{Mouse identifier within cage}
#'   \item{Treatment}{Treatment group}
#'   \item{Dose}{Dose level in arbitrary units (0, 10, 50, 100)}
#'   \item{Length}{Tumor length measurement in mm}
#'   \item{Width}{Tumor width measurement in mm}
#'   \item{Survival_Censor}{Censoring indicator (1 = event/death, 0 = censored/alive)}
#' }
#' @source Synthetic data created for demonstration purposes
#' @export
"dose_levels_synthetic_data"