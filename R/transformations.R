#' Geometric mean of a vector
#'
#' Calculate the geometric mean of a vector
#' Stolen from \url{https://stackoverflow.com/a/25555105}.
#'
#' @param x numeric vector of non-negative values
#' @param na.rm logical. Remove NAs from calculation
#' @param zero_propagate logical. Should zeros be included in the calculation. Default FALSE.
#' @return geometric mean of the vector
#' @export
#' @examples
#' x <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 1500)
#'
#' gm <- geometric_mean(x)
geometric_mean <- function(x, zero_propagate = FALSE, na.rm = TRUE) {
  if (any(x < 0, na.rm = TRUE)) {
    return(NaN)
  }
  if (zero_propagate) {
    if (any(x == 0, na.rm = TRUE)) {
      return(0)
    }
    exp(mean(log(x), na.rm = na.rm))
  } else {
    exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
  }
}


#' Centered Log-ratio transformation
#'
#' Calculate the centered log-ratio transformation of a vector
#'
#' @param x numeric vector
#' @param base integer base of the log function. Default 2.
#' @export
#' @examples
#' x <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 1500)
#'
#' clr_transformed <- clr(x)
clr <- function(x, base = 2) {
  x <- log((x / geometric_mean(x)), base)
  x[!is.finite(x) | is.na(x)] <- 0.0

  return(x)
}
