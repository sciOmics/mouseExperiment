#' Load namespace on package attach
#'
#' Ensures required packages are attached when mouseExperimentDashboard is loaded
#'
#' @keywords internal
.onLoad <- function(libname, pkgname) {
  # Load shiny and other required packages silently
  requireNamespace("shiny", quietly = TRUE)
  requireNamespace("mouseExperiment", quietly = TRUE)
}

#' Null-coalescing operator
#'
#' Returns \code{lhs} unless it is \code{NULL}, in which case returns \code{rhs}.
#' @keywords internal
`%||%` <- function(lhs, rhs) if (!is.null(lhs)) lhs else rhs
