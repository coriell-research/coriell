#' Calculate threshold value for unimodal distribution
#' 
#' Calculate a threshold value of a unimodal distribution based on the method
#' described by \href{https://users.cs.cf.ac.uk/Paul.Rosin/resources/papers/unimodal2.pdf"}{Rosin 2001} 
#' for image thresholding. Conceptually, the method attempts to draw a line from 
#' the peak of the distribution to the tail and calculates the maximum distance
#' from a point to that line. The point where the distance is maximixed is the
#' threshold. 
#' @param x numeric vector of x-coordinates
#' @param y numeric vector of y-coordinates
#' @return named vector containing "x" and "y" position of calculated threshold
#' @export
#' @examples 
#' # increasing exponential
#' x <- 1:100
#' y <- x^exp(3)
#' plot(x, y)
#' abline(v = unimodal_threshold(x, y)[1])
#' 
#' # decreasing exponential
#' x2 <- 1:100
#' y2 <- rev(x2)^exp(3)
#' plot(x2, y2)
#' abline(v = unimodal_threshold(x2, y2)[1])
#'
#' # with bump at right side
#' x3 <- 1:10
#' y3 <- c(1, 1, 1, 1, 2, 3, 5, 6, 4, 1)
#' plot(x3, y3)
#' abline(v = unimodal_threshold(x3, y3)[1])
#'
#' # with bump on left
#' x4 <- 1:10
#' y4 <- rev(y3)
#' plot(x4, y4)
#' abline(v = unimodal_threshold(x4, y4)[1])
#'
unimodal_threshold <- function(x, y) {
  stopifnot("x and y must be same length" = length(x) == length(y))
  stopifnot("Duplicated values in x" = !anyDuplicated(x))
  stopifnot("NA values in x" = all(!is.na(x)))
  stopifnot("NA values in y" = all(!is.na(y)))
  
  # check in which direction the distribution lies
  test_m <- coef(lm(y ~ x))[[2]]
  
  # if slope is positive then lower x is min(x) and upper x is x[max(y)]
  # if slope is negative then lower x is x[max(y)] and upper x is max(x)
  if (test_m > 0) {
    x1 <- min(x)
    y1 <- y[which.min(x)]
    x2 <- x[which.max(y)]
    y2 <- max(y)
  } else if (test_m < 0) {
    x1 <- x[which.max(y)]
    y1 <- max(y)
    x2 <- max(x)
    y2 <- y[which.max(x)]
  } else {
    stop("Slope of line is zero")
  }
  
  # points used to generate slope and intercept values
  X <- c(x1, x2)
  Y <- c(y1, y2)
  
  # use lm to get coefs for line
  line <- lm(Y ~ X)
  A <- coef(line)[[2]] # slope
  C <- coef(line)[[1]] # intercept
  B <- -1
  
  # distance from a point(m, n) to the line given by Y ~ X
  d <- function(m, n) abs((A * m) + (B * n) + C) / sqrt(A^2 + B^2)
  
  # limit to points in range
  xr <- x[x >= x1 & x <= x2]
  yr <- y[x >= x1 & x <= x2]
  
  distances <- d(xr, yr)

  return(c("x" = xr[which.max(distances)],
           "y" = yr[which.max(distances)]))
}