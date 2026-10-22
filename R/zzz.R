
#' Null-coalescing operator
#'
#' Returns \code{lhs} unless it is \code{NULL}, in which case returns \code{rhs}.
#' @noRd
`%||%` <- function(lhs, rhs) if (!is.null(lhs)) lhs else rhs
