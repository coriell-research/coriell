#' Parallel coordinates plot of expression data
#'
#' The parallel coordinates plot will display a line plot showing the expression
#' value for gene on the y-axis by each sample on the x-axis.
#'
#' @param x gene by sample matrix or \code{SummarizedExperiment} object
#' @param metadata data.frame containing metadata per sample. rownames of metadata
#' @param assay If \code{SummarizedExperiment} is supplied what assay to plot. Default "counts"
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
plot_parallel.default <- function(x) {
  stop("Object of type ", class(x), " is not supported by this function")
}

#' @rdname plot_parallel
#' @export
plot_parallel.matrix <- function(x, metadata = NULL, colBy = NULL,
                                 removeVar = NULL, ...) {
  stopifnot("colnames(x) do not match rownames(metadata)" = all(colnames(x) == rownames(metadata)))
  stopifnot("colBy must be a column in metadata" = colBy %in% colnames(metadata))
  stopifnot("non-numeric columns in x" = all(apply(x, 2, is.numeric)))
  
  if (!is.null(removeVar))
    x <- remove_var(x, p = removeVar)

  # Coerce data into plot-friendly shape
  dt <- data.table::as.data.table(x, keep.rownames = ".feature")
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

  if (is.null(colBy)) {
    colBy <- ".sample"
  }

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
plot_parallel.data.frame <- function(x, metadata = NULL, colBy = NULL,
                                     removeVar = NULL, ...) {
  if (is(x, "tbl_df") || is(x, "data.table")) {
    stop("You supplied a tibble or a data.table. Please use base::data.frame objects with rownames(x) == colnames(metadata)")
  }
  m <- data.matrix(x)
  plot_parallel.matrix(m, metadata, colBy, removeVar, ...)
}


#' @rdname plot_parallel
#' @export
plot_parallel.SummarizedExperiment <- function(x, assay = "counts", colBy = NULL,
                                               removeVar = NULL, ...) {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("SummarizedExperiment package is not installed.")
  }
  m <- SummarizedExperiment::assay(x, assay)
  plot_parallel.matrix(m, metadata, colBy, removeVar, ...)
}
