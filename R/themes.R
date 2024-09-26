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
#'   C = sample(c("Group 1", "Group 2", "Group 3"), 100, replace = TRUE)
#' )
#'
#' # Default ggplot2
#' ggplot(df, aes(X, Y, color = C)) +
#'   geom_point() +
#'   labs(
#'     title = "ggplot2 default plot",
#'     subtitle = "subtitle for default plot",
#'     x = "X Label",
#'     y = "Y Label"
#'   ) + 
#'   facet_wrap(~C)
#'
#' # With theme_coriell()
#' ggplot(df, aes(X, Y, color = C)) +
#'   geom_point() +
#'   labs(
#'     title = "Plot using theme_coriell()",
#'     subtitle = "Show facets and plot legend",     
#'     x = "X Label",
#'     y = "Y Label"
#'   ) +
#'   facet_wrap(~C) +
#'   theme_coriell()
#'
theme_coriell <- function() {
  ggplot2::theme_light() +
    ggplot2::theme(

      text = ggplot2::element_text(family = "Helvetica"),
      
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = 16, face = "bold"),
      legend.text = ggplot2::element_text(size = 14),
      
      plot.title = ggplot2::element_text(size = 22, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 18),
      
      axis.title = ggplot2::element_text(size = 16),
      axis.text.x = ggplot2::element_text(size = 16),
      axis.text.y = ggplot2::element_text(size = 16),
      
      strip.text.x = ggplot2::element_text(size = 16, color = "black", face = "bold"),
      strip.text.y = ggplot2::element_text(size = 16, color = "black", face = "bold"),
      strip.background = ggplot2::element_rect(fill = "white", linewidth = 2),
      
      panel.border = ggplot2::element_blank()
      
    )
}
