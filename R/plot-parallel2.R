#' Parallel coordinates plot
#'
#' Create a parallel coordinates (scaled expression data on y-axis, samples on x-axis)
#' for a matrix of data. This function calls \code{MASS::parcoord()} after optionally
#' removing the low variance features and scaling the data.
#'
#' @param x feature x sample matrix
#' @param remove_var numeric. What proportion of low variance features to remove from the matrix
#' before plotting. Default 0.9
#' @param scale bool. Should the data be standardized prior to plotting. Default TRUE
#' @param color character. Color of the lines. Default "black"
#' @param alpha numeric. Alpha value used to add transparency to lines. Default 0.1
#' @param ... Additional variables passed to \code{MASS:parcoord()}
#' @returns Parallel coordinate plot of scaled matrix data
#' @export
#' @examples
#'
#' plot_parallel2(GSE161650_lc, main = "10% Most Variable Genes")
#'
plot_parallel2 <- function(
  x,
  remove_var = 0.9,
  scale = TRUE,
  color = "black",
  alpha = 0.1,
  ...
) {
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("MASS package is required.")
  }

  m <- coriell::remove_var(x, p = remove_var)

  if (isTRUE(scale)) {
    m <- t(scale(t(m), center = TRUE, scale = TRUE))
  }

  MASS::parcoord(
    as.matrix(m),
    col = grDevices::adjustcolor(color, alpha.f = alpha),
    ...
  )
}
