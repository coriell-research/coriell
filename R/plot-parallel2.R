#' Parallel coordinates plot
#'
#' Create a parallel coordinates (scaled expression data on y-axis, samples on x-axis)
#' for a matrix of data.
#'
#' @param x feature x sample matrix
#' @param remove_var numeric. What proportion of low variance features to remove from the matrix
#' before plotting. Default 0.9
#' @param scale bool. Should the data be standardized prior to plotting. Default TRUE
#' @param color character. Color of the lines. Default "black"
#' @param alpha numeric. Alpha value used to add transparency to lines. Default 0.1
#' @param plot_title character. title of the plot. Default NULL
#' @param x_label character. x-axis label
#' @param y_label character. y-axis label
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
  plot_title = NULL,
  x_label = NULL,
  y_label = NULL,
  ...
) {
  m <- coriell::remove_var(x, p = remove_var)

  if (isTRUE(scale)) {
    m <- t(scale(t(m), center = TRUE, scale = TRUE))
  }

  plot.new()
  plot.window(xlim = c(1, ncol(m)), ylim = range(m, na.rm = TRUE))
  for (i in 1:nrow(m)) {
    lines(
      x = 1:ncol(m),
      y = as.numeric(m[i, ]),
      col = grDevices::adjustcolor(color, alpha.f = alpha)
    )
  }
  axis(side = 1, at = 1:ncol(m), labels = colnames(m), lty = 3, las = 2)
  axis(side = 2)
  title(
    main = plot_title,
    ylab = y_label,
    xlab = x_label,
  )
}
