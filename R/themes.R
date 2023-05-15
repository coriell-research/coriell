#' ggplot2 theme for coriell package
#'
#' Based on \code{ggplot2::theme_light()} this theme will adjust the plot title,
#' legend position (bottom), axis titles, axis ticks, and strip text.
#' @export
#' @examples
#' library(ggplot2)
#'
#' df <- data.frame(
#'   X = 1:100,
#'   Y = rnorm(100),
#'   C = sample(LETTERS[1:4], 100, replace = TRUE)
#' )
#'
#' # Default ggplot2
#' ggplot(df, aes(X, Y, color = C)) +
#'   geom_point() +
#'   labs(
#'     title = "ggplot2 default plot",
#'     subtitle = "subtitle for default plot"
#'   )
#'
#' # With theme_coriell()
#' ggplot(df, aes(X, Y, color = C)) +
#'   geom_point() +
#'   labs(
#'     title = "Plot with theme_coriell()",
#'     subtitle = "Subtitle for plot with theme_coriell()"
#'   ) +
#'   theme_coriell()
#'
theme_coriell <- function() {
  ggplot2::theme_light() +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title = ggplot2::element_text(size = 22, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 16),
      axis.title = ggplot2::element_text(size = 14),
      strip.text = ggplot2::element_text(size = 12, face = "bold"),
      axis.text.x = ggplot2::element_text(size = 14, angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 14)
    )
}
