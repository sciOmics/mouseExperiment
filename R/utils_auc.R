#' Build a composite mouse key from two or three fields
#'
#' Uses \code{"|||"} as the separator, which is safe for IDs, treatment names,
#' and cage labels containing letters, digits, spaces, hyphens, underscores, or
#' dots. Pass arguments in a consistent order (e.g. Treatment, ID, Cage) and
#' use \code{split_mouse_key()} to reverse the operation.
#'
#' @param ... Character vectors to paste together (recycled).
#' @return Character vector of composite keys.
#' @noRd
#' @keywords internal
make_mouse_key <- function(...) paste(..., sep = "|||")

#' Split a composite mouse key back into its components
#'
#' @param key Single character string produced by \code{make_mouse_key()}.
#' @return Character vector of components.
#' @noRd
#' @keywords internal
split_mouse_key <- function(key) strsplit(key, "|||", fixed = TRUE)[[1]]

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
