#' Parallel coordinates plot of expression data
#' 
#' The parallel coordinates plot will display a line plot showing the expression
#' value for gene on the y-axis by each sample on the x-axis.
#' 
#' @param x gene by sample matrix or \code{SummarizedExperiment} object
#' @param metadata data.frame containing metadata per sample. rownames of metadata
#' must match the colnames of the input matrix. Default NULL, each sample in the 
#' matrix will be plotted.
#' @param assay assay of \code{SummarizedExperiment} object to be plotted. Default "counts".
#' @param colBy metadata column used to color lines. Default NULL, every sample 
#' will get its own color.
#' @param removeVar Remove this proportion of features based on the variance 
#' across rows. Default NULL, all features are plotted. 
#' @param ... Additional parameters passed to \code{ggplot2::geom_line()}
#' @import data.table
#' @return ggplot object
#' @export
#' @examples
#' # Create metadata for plotting
#' metadata <- data.frame(row.names = colnames(GSE161650_lc))
#' metadata$Group <- rep(c("DMSO", "THZ1"), each = 3)
#' 
#' # Plot the PCP for each sample -- passing alpha value to geom_line()
#' plot_parallel(GSE161650_lc, alpha = 0.01) + 
#'   theme_coriell()
#'
#' # Plot the PCP by coloring each sample by Group from metadata
#' plot_parallel(GSE161650_lc, metadata, colBy = "Group", alpha = 0.01) + 
#'   theme_coriell()
plot_parallel <- function(x, ...) UseMethod("plot_parallel")


#' @rdname plot_parallel
#' @export
plot_parallel.matrix <- function(x, metadata = NULL, colBy = NULL, 
                                 removeVar = NULL, ...) {
  stopifnot("colnames(x) do not match rownames(metadata)" = all(colnames(x) == rownames(metadata)))
  stopifnot("colBy must be a column in metadata" = colBy %in% colnames(metadata))
  stopifnot("non-numeric columns in x" = all(apply(x, 2, is.numeric)))
  
  mat <- x
  if (!is.null(removeVar)) {
    stopifnot("RemoveVar must be between 0 and 1" = removeVar > 0 & removeVar < 1)
    v <- apply(mat, 1, var)
    o <- order(v, decreasing = TRUE)
    mat <- head(mat[o, ], n = nrow(mat) * (1 - removeVar))
  }
  
  # Coerce data into plot-friendly shape
  dt <- data.table::as.data.table(mat, keep.rownames = ".feature")
  dt.m <- data.table::melt(
    dt, 
    id.vars = ".feature", 
    variable.name = ".sample",
    value.name = ".value",
    variable.factor = FALSE
  )
  
  if (!is.null(metadata)) {
    meta <- data.table::as.data.table(metadata, keep.rownames = ".sample")
    dt.m <- dt.m[meta, on = ".sample", nomatch = NULL]
  }
  
  if (is.null(colBy))
    colBy <- ".sample"
  
  ggplot2::ggplot(dt.m, ggplot2::aes_string(x = ".sample", y = ".value", group = ".feature")) +
    ggplot2::geom_line(ggplot2::aes_string(color = colBy), ...) +
    ggplot2::labs(
      x = NULL, 
      y = NULL,
      color = if (colBy == ".sample") "Sample" else colBy
      ) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(alpha = 1)))
}


#' @rdname plot_parallel
#' @export
plot_parallel.SummarizedExperiment <- function(x, assay = "counts", colBy = NULL, 
                                               removeVar = NULL, ...) {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE))
    stop("SummarizedExperiment package is not installed.")
  
  M <- SummarizedExperiment::assay(x, assay)
  mat <- M
  
  if (!is.null(removeVar)) {
    stopifnot("RemoveVar must be between 0 and 1" = removeVar > 0 & removeVar < 1)
    v <- apply(mat, 1, var)
    o <- order(v, decreasing = TRUE)
    mat <- head(mat[o, ], n = nrow(mat) * (1 - removeVar))
  }
  
  # Coerce data into plot-friendly shape
  dt <- data.table::as.data.table(mat, keep.rownames = ".feature")
  dt.m <- data.table::melt(
    dt, 
    id.vars = ".feature", 
    variable.name = ".sample",
    value.name = ".value",
    variable.factor = FALSE
  )
  
  meta <- data.table::setDT(
    data.frame(SummarizedExperiment::colData(x)), 
    keep.rownames = ".sample"
  )
  dt.m <- dt.m[meta, on = ".sample", nomatch = NULL]
  
  if (is.null(colBy))
    colBy <- ".sample"
  
  ggplot2::ggplot(dt.m, ggplot2::aes_string(x = ".sample", y = ".value", group = ".feature")) +
    ggplot2::geom_line(ggplot2::aes_string(color = colBy), ...) +
    ggplot2::labs(
      x = NULL, 
      y = NULL,
      color = if (colBy == ".sample") "Sample" else colBy
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(alpha = 1)))
}
