#' Generate Random RGB Color Palette
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
  stopifnot("length(n) must be equal to 1" = length(n) == 1)
  stopifnot("alpha must be between 0 and 1" = alpha >= 0 & alpha <= 1)
  
  grDevices::rgb(red = stats::runif(n), 
                 green = stats::runif(n), 
                 blue = stats::runif(n), 
                 alpha = alpha)
}
