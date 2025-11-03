#' Plot results of UMAP
#'
#' This function creates a simple plot of the data.frame returned by the
#' \code{umap()} function.
#' @param df data.frame of UMAP embeddings and metadata
#' @param x First UMAP component. Default "UMAP1"
#' @param y Second UMAP component. Default "UMAP2"
#' @param colBy Column name of data.frame to color points by. Default NULL
#' @param shapeBy Column name of data.frame to shape the points by. Default NULL
#' @param pointSize Size of the points. Default 3
#' @param pointAlpha Alpha level of the points. Default 1
#' @param hline y-position of horizontal line. Default 0
#' @param vline x-position of vertical line. Default 0
#' @param hlineType linetype of the horizontal line. Default 2
#' @param vlineType linetype of the vertical line. Default 2
#' @export
#' @return ggplot object
#' @examples
#' # Create metadata for plotting
#' metadata <- data.frame(row.names = colnames(GSE161650_lc))
#' metadata$Group <- rep(c("DMSO", "THZ1"), each = 3)
#'
#' # PCA with PCAtools
#' p <- PCAtools::pca(GSE161650_lc, metadata, center = TRUE, scale = TRUE)
#' udata <- coriell::UMAP(p, n_neighbors = 2)
#' plot_umap(udata, colBy = "Group")
plot_umap <- function(
  df,
  x = "UMAP1",
  y = "UMAP2",
  colBy = NULL,
  shapeBy = NULL,
  pointSize = 3,
  pointAlpha = 1,
  hline = 0,
  vline = 0,
  hlineType = 2,
  vlineType = 2
) {
  if (!is.null(colBy)) {
    if (!(colBy %in% colnames(df))) {
      stop(colBy, " is not a column in df")
    }
  }
  if (!is.null(shapeBy)) {
    if (!(shapeBy %in% colnames(df))) {
      stop(shapeBy, " is not a column in df")
    }
  }

  ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data[[x]],
      y = .data[[y]],
      color = if (is.null(colBy)) NULL else .data[[colBy]],
      shape = if (is.null(shapeBy)) NULL else .data[[shapeBy]]
    )
  ) +
    ggplot2::geom_point(size = pointSize, alpha = pointAlpha) +
    ggplot2::geom_hline(yintercept = hline, linetype = hlineType) +
    ggplot2::geom_vline(xintercept = vline, linetype = vlineType) +
    ggplot2::labs(
      color = if (is.null(colBy)) NULL else colBy,
      shape = if (is.null(shapeBy)) NULL else shapeBy
    ) +
    theme_coriell() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank()
    )
}
