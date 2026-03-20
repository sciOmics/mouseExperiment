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
