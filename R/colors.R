#' Generate random RGB color palette
#'
#' This function will randomly generate N colors from the rgb colorspace. It
#' makes no attempt at ensuring the colors are distinct or well separated.
#'
#' @param n numeric. Number of random colors to generate
#' @param alpha numeric. Transparency level of color palette (0-1). Default 1.0
#' @export
#' @return vector of hex values for RGB colors of length n
random_rgb_palette <- function(n, alpha = 1.0) {
  stopifnot("n must be numeric" = is.numeric(n))
  stopifnot("n must be greater than 0" = n > 0)
  stopifnot("length(n) must be equal to 1" = length(n) == 1)
  stopifnot("alpha must be between 0 and 1" = alpha >= 0 & alpha <= 1)

  pal <- grDevices::rgb(
    red = stats::runif(n),
    green = stats::runif(n),
    blue = stats::runif(n),
    alpha = alpha
    )
  
  return(pal)
}

#' Generate a distinct RGB color palette
#'
#' This function uses \code{kmeans()} over the RGB colorspace to generate N distinct
#' RGB colors.
#'
#' This function uses \emph{very} inaccurate defaults for the \code{kmeans()} function in
#' the interest of speed. It's usually not a problem if \code{kmeans()} does not
#' converge (colors are distinct enough for most purposes). If you get warnings,
#' or you don't like the colors produced, you can modify the default arguments to the
#' \code{kmeans()} function by passing in additional arguments to \code{...}. For example,
#' increasing iterations can be done by passing \code{iter.max = 100}. Changing the
#' kmeans algorithm can be done by specifying \code{algorithm = "MacQueen"}. Changing
#' these arguments may eliminate warnings or produce more distinct colors.
#'
#' @param n numeric. Number of colors to generate
#' @param alpha numeric. Transparency level of the color palette (0-1). Default 1.0
#' @param ... arguments passed to \code{kmeans()}
#' @export
#' @return vector of distinct RGB colors
distinct_rgb_palette <- function(n, alpha = 1.0, ...) {
  stopifnot("n must be numeric" = is.numeric(n))
  stopifnot("n must be greater than 0" = n > 0)
  stopifnot("length(n) must be equal to 1" = length(n) == 1)
  stopifnot("alpha must be between 0 and 1" = alpha >= 0 & alpha <= 1)

  default_args <- list(algorithm = "Lloyd", iter.max = 1)
  user_args <- list(...)
  default_args[names(user_args)] <- user_args

  space <- expand.grid(rep(list(0:255), 3))
  res <- do.call(stats::kmeans, c(list(x = space, centers = n), default_args))
  centers <- floor(res$centers)
  pal <- apply(centers, 1, function(x) rgb(x[1], x[2], x[3], maxColorValue = 255, alpha = (255 * alpha)))

  return(pal)
}
