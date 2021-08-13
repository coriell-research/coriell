#' Calculate threshold value on ranked input
#' 
#' Calculate a threshold value of a vector based on the unimodal threshold method
#' described by \href{https://users.cs.cf.ac.uk/Paul.Rosin/resources/papers/unimodal2.pdf"}{Rosin 2001} 
#' for image thresholding. Conceptually, the method attempts to draw a line from 
#' the peak of the distribution to the tail and calculates the maximum distance
#' from a point to that line. The point where the distance is maximized is the
#' threshold.
#' 
#' The function takes a numeric vector and ranks the values using \code{data.table::frank}.
#' with \code{ties.method = "first"}. The returned values for the threshold "x" and "y", 
#' indicate the rank of the data point and the value at that rank, respectively.
#' 
#' @param x numeric vector of values.
#' @param show_plot logical. default (FALSE). Show the scatter plot for the given distribution. 
#' @return named vector containing "x" and "y" position of calculated threshold
#' @import data.table
#' @export
#' @examples
#' # simulate exponential counts
#' x <- (1:100)^exp(3)
#' rank_threshold(x, show = TRUE)
rank_threshold <- function(x, show = FALSE) {
  stopifnot("NA values present in x" = all(!is.na(x)))
  stopifnot("x must be numeric" = is.numeric(x))
  rnk <- data.table::frank(x, ties.method = "first")
  
  # find the min and max points
  x1 <- min(rnk)
  y1 <- x[which.min(rnk)]
  x2 <- rnk[which.max(x)]
  y2 <- max(x)
  
  # generate slope and intercept values from min and max points
  X <- c(x1, x2)
  Y <- c(y1, y2)
  
  # use lm to get coefs for line
  line <- lm(Y ~ X)
  A <- coef(line)[[2]] # slope
  C <- coef(line)[[1]] # intercept
  B <- -1
  
  # distance from a point(m, n) to the line given by Y ~ X
  d <- function(m, n) abs((A * m) + (B * n) + C) / sqrt(A^2 + B^2)
  distances <- d(rnk, x)
  x_thresh <- rnk[which.max(distances)]
  y_thresh <- x[which.max(distances)]
  
  if (show) {
    plot(rnk, x, xlab = "rank", ylab = "value", main = "Estimated Thresholds (red lines)")
    abline(reg = line, col = "black", lty = 2)
    abline(v = x_thresh, col = "red", lty = 3)
    abline(h = y_thresh, col = "red", lty = 3)
  }
  
  return(c("x" = x_thresh, "y" = y_thresh))

}
