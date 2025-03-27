
#' Synthetic tumor growth data
#' 
#' A dataset containing synthetic tumor volume measurements over time
#' for different treatment groups.
#' 
#' @format A data frame with 120 rows and 7 variables:
#' \describe{
#'   \item{Mouse_ID}{Mouse identifier, with format M followed by a number}
#'   \item{Day}{Day of measurement, starting from day 0}
#'   \item{Treatment}{Treatment group (Control, Treatment A, Treatment B)}
#'   \item{Volume}{Tumor volume measurement}
#'   \item{Group}{Same as Treatment, an alternative name for the treatment group}
#'   \item{ID}{Numeric identifier for the mouse, without the "M" prefix}
#'   \item{Cage}{Cage number where the mouse is housed}
#' }
#' @source Synthetic data generated using random number generation
#' @examples
#' data(synthetic_data)
#' head(synthetic_data)
"synthetic_data"

#' Example data
#' 
#' A simple example dataset with x and y coordinates.
#' 
#' @format A data frame with 10 rows and 2 variables:
#' \describe{
#'   \item{x}{x-coordinate, sequence from 1 to 10}
#'   \item{y}{y-coordinate, x plus random noise}
#' }
#' @source Synthetic data generated using random number generation
#' @examples
#' data(my_data)
#' plot(my_data$x, my_data$y)
"my_data"

#' Combination treatment synthetic data
#' 
#' A dataset containing synthetic tumor volume measurements over time
#' for different combination treatment groups.
#' 
#' @format A data frame with rows for each mouse at each timepoint and 6 variables:
#' \describe{
#'   \item{Mouse_ID}{Mouse identifier, with format M followed by a number}
#'   \item{Day}{Day of measurement, starting from day 0}
#'   \item{Treatment}{Treatment group (Control, aPD1, HDACi, HDACi + PD1)}
#'   \item{Volume}{Tumor volume measurement}
#'   \item{ID}{Numeric identifier for the mouse, without the "M" prefix}
#'   \item{Cage}{Cage number where the mouse is housed}
#' }
#' @source Synthetic data generated using random number generation with treatment-specific effects
#' @examples
#' data(combo_treatment_synthetic_data)
#' head(combo_treatment_synthetic_data)
"combo_treatment_synthetic_data"

#' Combination treatment schedule
#' 
#' A dataset specifying the dosing schedule for combination treatments.
#' 
#' @format A data frame with 20 rows and 3 variables:
#' \describe{
#'   \item{Treatment}{Treatment group (Control, aPD1, HDACi, HDACi + PD1)}
#'   \item{Day}{Day of dose administration}
#'   \item{Dose}{Dose amount}
#' }
#' @source Synthetic treatment schedule
#' @examples
#' data(combo_treatment_schedule)
#' head(combo_treatment_schedule)
"combo_treatment_schedule"

#' Dose levels synthetic data
#' 
#' A dataset containing synthetic tumor volume measurements over time
#' for different dose levels of a single drug.
#' 
#' @format A data frame with rows for each mouse at each timepoint and 7 variables:
#' \describe{
#'   \item{Mouse_ID}{Mouse identifier, with format M followed by a number}
#'   \item{Day}{Day of measurement, starting from day 0}
#'   \item{Treatment}{Treatment name (always "Drug X")}
#'   \item{Dose}{Dose level (0, 10, 25, 50, 100)}
#'   \item{Volume}{Tumor volume measurement}
#'   \item{ID}{Numeric identifier for the mouse, without the "M" prefix}
#'   \item{Cage}{Cage number where the mouse is housed}
#' }
#' @source Synthetic data generated using random number generation with dose-dependent effects
#' @examples
#' data(dose_levels_synthetic_data)
#' head(dose_levels_synthetic_data)
"dose_levels_synthetic_data"

#' Dose levels treatment schedule
#' 
#' A dataset specifying the dosing schedule for different dose levels.
#' 
#' @format A data frame with 20 rows and 4 variables:
#' \describe{
#'   \item{Treatment}{Treatment name (always "Drug X")}
#'   \item{Dose}{Dose group level (0, 10, 25, 50, 100)}
#'   \item{Day}{Day of dose administration}
#'   \item{Administered_Dose}{Actual dose administered}
#' }
#' @source Synthetic treatment schedule
#' @examples
#' data(dose_levels_treatment_schedule)
#' head(dose_levels_treatment_schedule)
"dose_levels_treatment_schedule"
