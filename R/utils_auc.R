#' Calculate Area Under the Curve (Trapezoidal Rule)
#'
#' Computes the area under the curve for paired time-value vectors using the
#' trapezoidal rule. Data are automatically sorted by time before calculation.
#'
#' @param time_values Numeric vector of time points.
#' @param volume_values Numeric vector of corresponding measurements (e.g. tumor volume).
#'
#' @return A single numeric AUC value, or \code{NA} if fewer than 2 valid points.
#'
#' @details
#' The trapezoidal rule approximates the integral by summing trapezoids formed
#' between consecutive time points:
#' \deqn{AUC = \sum_{i=2}^{n}{\frac{(v_{i-1} + v_i)(t_i - t_{i-1})}{2}}}
#'
#' @examples
#' calculate_auc(c(0, 3, 7, 14), c(100, 150, 300, 800))
#'
#' @export
calculate_auc <- function(time_values, volume_values) {
  # Input validation

  if (length(time_values) != length(volume_values)) {
    stop("time_values and volume_values must have the same length")
  }

  # Remove NA pairs
  valid <- !is.na(time_values) & !is.na(volume_values)
  time_values <- time_values[valid]
  volume_values <- volume_values[valid]

  # Need at least 2 points to calculate AUC
  if (length(time_values) < 2) {
    return(NA_real_)
  }

  # Sort by time
  ord <- order(time_values)
  time_values <- time_values[ord]
  volume_values <- volume_values[ord]

  # Vectorised trapezoidal rule
  dt <- diff(time_values)
  avg_height <- (volume_values[-length(volume_values)] + volume_values[-1]) / 2
  sum(dt * avg_height)
}
