#' Get the axis limits of a ggplot2 object
#'
#' Internal helper function used to get the axis limits of a plot
#' @param p \code{ggplot} plot object
#' @keywords internal
#' @export
#' @return list with elements for x and y limits
get_axis_limits <- function(p) {
  
  b <- ggplot2::ggplot_build(p)
  
  x_range <- b$layout$panel_scales_x[[1]]$range$range
  y_range <- b$layout$panel_scales_y[[1]]$range$range
  xmin <- x_range[1]
  xmax <- x_range[2]
  ymin <- y_range[1]
  ymax <- y_range[2]
  
  # Get limits if the ranges have been reset by coord_
  x_limits <- b$layout$coord$limits$x
  y_limits <- b$layout$coord$limits$y
  
  if (!is.null(x_limits)) {
    xmin <- x_limits[1]
    xmax <- x_limits[2]
  }

  if (!is.null(y_limits)) {
    ymin <- y_limits[1]
    ymax <- y_limits[2]
  }

  list(
    "x_min" = xmin,
    "x_max" = xmax,
    "y_min" = ymin,
    "y_max" = ymax
  )
}
