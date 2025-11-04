#' Create a density plot from the columns of a matrix
#'
#' This is an alternate to the \code{plot_density()} function using base graphics.
#' The advantage of this function is that it should work with most matrix-like
#' objects. This is especially helpful for \code{DelayedArrays}, which
#' \code{plot_density()} does not handle.
#'
#' @param x feature x sample matrix or matrix-like object
#' @param metadata data.frame containing metadata per sample. rownames of metadata
#' must match the colnames of the input matrix. Default NULL, each sample in the
#' matrix will be plotted.
#' @param col_by metadata column used to color density lines. Default NULL,
#' each sample in the matrix will be plotted.
#' @param hcl_palette color palette applied to 'col_by' variable. One of the
#' \code{hcl.pals()}. Default "Dark 3"
#' @param plot_title title of the plot. Default NULL
#' @param x_label x-axis title. Default NULL
#' @param y_label y-axis title. Default "Density"
#' @param line_weight line weight of the density lines. Default 1
#' @param line_type line type of the density lines. Default 1
#' @param show_legend should the legend be displayed on the plot if fill_by is set? Default TRUE
#' @param legend_position location keyword for the legend. One of "bottomright",
#'  "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and
#'  "center". Default "topright"
#' @param legend_horiz logical; if TRUE, set the legend horizontally rather than vertically.
#' Default FALSE
#' @param legend_cex character expansion factor for legend text relative to current par. Default 0.8
#' @param legend_ncol number of columns to set the legend items. Default 1. This is not used if
#' legend_horiz=TRUE
#' @param ... Additional arguments not currently used.
#'
#' @returns density plot of matrix columns
#' @export
#'
#' @examples
#'
#' # Create metadata for plotting
#' metadata <- data.frame(row.names = colnames(GSE161650_lc))
#' metadata$Group <- rep(c("DMSO", "THZ1"), each = 3)
#'
#' # Plot the density by sample
#' plot_density2(GSE161650_lc)
#'
#' # Color each sample by their Group in metadata
#' plot_density2(
#'   GSE161650_lc,
#'   metadata,
#'   col_by = "Group",
#'   x_label = "log2 CPMs",
#'   plot_title = "Distribution of logCPM values"
#'   )
#'
plot_density2 <- function(
  x,
  metadata = NULL,
  col_by = NULL,
  hcl_palette = "Dark 3",
  plot_title = NULL,
  x_label = NULL,
  y_label = "Density",
  line_weight = 1,
  line_type = 1,
  show_legend = TRUE,
  legend_position = "topright",
  legend_horiz = FALSE,
  legend_cex = 0.8,
  legend_ncol = 1,
  ...
) {
  stopifnot(
    "colnames(x) do not match rownames(metadata)" = all(
      colnames(x) == rownames(metadata)
    )
  )
  stopifnot(
    "col_by must be a column in metadata" = col_by %in% colnames(metadata)
  )

  rng <- range(x[is.finite(x)], na.rm = TRUE)

  # Compute densities for each column
  d <- apply(
    x,
    2,
    function(j) {
      stats::density(as.numeric(j), from = rng[1], to = rng[2], na.rm = TRUE)$y
    },
    simplify = TRUE
  )
  x_values <- seq(from = rng[1], to = rng[2], length.out = nrow(d))

  # Determine colors
  if (!is.null(col_by)) {
    g <- as.factor(metadata[[col_by]])
    cols <- grDevices::hcl.colors(n = nlevels(g), palette = hcl_palette)
    names(cols) <- levels(g)
  }

  plot.new()
  plot.window(xlim = rng, ylim = range(d))
  grid(col = "lightgrey", lty = 3)
  for (i in 1:ncol(d)) {
    line_color <- if (is.null(col_by)) "black" else cols[g[i]]
    lines(
      x_values,
      d[, i],
      col = line_color,
      lwd = line_weight,
      lty = line_type
    )
  }
  axis(1, lwd = 0)
  axis(2, lwd = 0)
  title(main = plot_title, xlab = x_label, ylab = y_label)

  if (!is.null(col_by)) {
    legend(
      legend_position,
      legend = levels(g),
      col = cols,
      lty = line_type,
      lwd = line_weight,
      bty = "n",
      horiz = legend_horiz,
      ncol = legend_ncol,
      cex = legend_cex
    )
  }
}
