#' Show boxplots for columns of data in a matrix
#'
#' This is an alternative to \code{plot_boxplot()} using base R. This function
#' should work with most matrix or matrix-like objects (which is especially
#' useful for \code{DelayedArray} objects). The function also includes an
#' option for computing the relative-log-expression of the input data, given the
#' input data is on the log-scale.
#'
#' @param x feature x sample matrix or matrix-like object.
#' @param metadata data.frame containing metadata per sample. rownames of metadata
#' must match the colnames of the input matrix. Default NULL, each sample in the
#' matrix will be plotted.
#' @param fill_by metadata column used to color boxplots by the grouping variable.
#' Default NULL, each sample in the matrix will be plotted
#' @param rle should the relative-log-expression value be plotted. Requires
#' input matrix to be on the log-scale. Default = FALSE
#' @param hcl_palette color palette applied to 'fill_by' variable. One of the
#' \code{hcl.pals()}. Default "Dark 3"
#' @param plot_title title of the plot. Default NULL
#' @param x_label title of th x-axis. Default NULL
#' @param y_label title of the y-axis. Default NULL
#' @param show_outliers should boxplot outliers be shown. Default TRUE
#' @param show_legend should the legend be drawn on the plot when fill_by is set? Default TRUE
#' @param legend_position location keyword for the legend. One of "bottomright",
#'  "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and
#'  "center". Default "topright"
#' @param legend_horiz logical; if TRUE, set the legend horizontally rather than vertically.
#' Default TRUE
#' @param legend_cex character expansion factor for legend text relative to current par. Default 0.8
#' @param legend_ncol number of columns to set the legend items. Default 1. This is not used if
#' legend_horiz=TRUE
#' @param ... additional arguments passed to \code{bxp()}
#'
#' @returns boxplots of the column data
#' @export
#'
#' @examples
#'
#' # Create metadata for plotting
#' metadata <- data.frame(row.names = colnames(GSE161650_lc))
#' metadata$Group <- rep(c("DMSO", "THZ1"), each = 3)
#'
#' # Plot the boxplot by sample
#' plot_boxplot2(GSE161650_lc)
#'
#' # Plot the boxplot by coloring each Group and transforming values
#' #  to relative-log-expression
#' plot_boxplot2(
#'   GSE161650_lc,
#'   metadata,
#'   fill_by = "Group",
#'   rle = TRUE,
#'   plot_title = "Relative log-expression",
#'   y_label = "RLE",
#'   show_outliers = FALSE
#'   )
#'
plot_boxplot2 <- function(
  x,
  metadata = NULL,
  fill_by = NULL,
  rle = FALSE,
  hcl_palette = "Dark 3",
  plot_title = NULL,
  x_label = NULL,
  y_label = NULL,
  show_outliers = TRUE,
  show_legend = TRUE,
  legend_position = "top",
  legend_horiz = TRUE,
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
    "fill_by must be a column in metadata" = fill_by %in% colnames(metadata)
  )

  if (isTRUE(rle)) {
    if (is(x, "DelayedMatrix")) {
      if (requireNamespace("DelayedMatrixStats", quietly = TRUE)) {
        message("Running DelayedMatrixStats...")
        m <- DelayedMatrixStats::rowMedians(x, na.rm = TRUE)
        x <- x - m
      }
    } else {
      if (requireNamespace("matrixStats", quietly = TRUE)) {
        message("running matrixStats...")
        m <- matrixStats::rowMedians(x, na.rm = TRUE)
        x <- x - m
      } else {
        m <- apply(x, 1, median, na.rm = TRUE)
        x <- x - m
      }
    }
  }

  # Manually compute boxplot stats (for DelayedArrays mostly)
  stats_list <- apply(x, 2, grDevices::boxplot.stats, simplify = FALSE)
  stats_matrix <- sapply(stats_list, '[[', 'stats')
  n_vec <- sapply(stats_list, '[[', 'n')
  out_vec <- unlist(sapply(stats_list, '[[', 'out'))
  out_group_vec <- rep(
    1:ncol(x),
    times = sapply(stats_list, function(x) length(x$out))
  )
  col_names <- colnames(x)

  bxp_data <- list(
    stats = stats_matrix,
    n = n_vec,
    out = out_vec,
    group = out_group_vec,
    names = col_names
  )

  if (!is.null(fill_by)) {
    g <- as.factor(metadata[[fill_by]])
    cols <- grDevices::hcl.colors(n = nlevels(g), palette = hcl_palette)
    names(cols) <- levels(g)
  }

  bxp(
    bxp_data,
    outline = show_outliers,
    main = plot_title,
    ylab = y_label,
    xlab = x_label,
    frame.plot = FALSE,
    boxfill = if (is.null(fill_by)) "lightgrey" else cols[g],
    ...
  )

  if (isTRUE(rle)) {
    abline(h = 0, lty = 2, col = "black")
  }

  if (!is.null(fill_by)) {
    if (isTRUE(show_legend)) {
      legend(
        legend_position,
        legend = levels(g),
        fill = cols,
        bty = "n",
        horiz = legend_horiz,
        ncol = legend_ncol,
        cex = legend_cex
      )
    }
  }
}
