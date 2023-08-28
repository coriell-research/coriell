#' Create boxplot from expression data
#'
#' Create a boxplot (or violin plot) of expression distributions for the given
#' expression matrix. Optionally plot the relative log expression of the matrix.
#'
#' @param x matrix of expression values or \code{SummarizedExperiment} object
#' @param assay If a SummarizedExperiment is supplied what assay is used. Default "counts"
#' @param metadata data.frame containing metadata per sample. rownames of metadata
#' must match the colnames of the input matrix. Default NULL, each sample in x
#' will be plotted individually.
#' @param fillBy metadata column used to fill boxplots. Default NULL,
#' each sample will be a distinct color.
#' @param rle Plot the relative log expression value. This option assumes
#' that the input matrix has already been logged. Be sure to take the log of the
#' input matrix prior to setting this option.
#' @param violin Plot the data as a violin plot instead of a boxplot. Default FALSE.
#' @param ... Additional parameters passed to \code{ggplot2::geom_boxplot()} or
#' \code{ggplot2::geom_violin()}
#' @return ggplot object
#' @export
#' @examples
#' # Create metadata for plotting
#' metadata <- data.frame(row.names = colnames(GSE161650_lc))
#' metadata$Group <- rep(c("DMSO", "THZ1"), each = 3)
#'
#' # Plot the boxplot by sample
#' plot_boxplot(GSE161650_lc) +
#'   theme_coriell()
#'
#' # Plot the boxplot by coloring each Group
#' plot_boxplot(GSE161650_lc, metadata, fillBy = "Group") +
#'   theme_coriell()
#'
#' # Create a violin plot after RLE transformation
#' plot_boxplot(GSE161650_lc, metadata, fillBy = "Group", rle = TRUE, violin = TRUE) +
#'   theme_coriell()
plot_boxplot <- function(x, ...) UseMethod("plot_boxplot")

#' @export
#' @rdname plot_boxplot
plot_boxplot.default <- function(x) {
  stop("Object of type ", class(x), " is not supported by this function")
}

#' @export
#' @rdname plot_boxplot
plot_boxplot.matrix <- function(x, metadata = NULL, fillBy = NULL, rle = FALSE,
                                violin = FALSE, ...) {
  stopifnot("colnames(x) do not match rownames(metadata)" = all(colnames(x) == rownames(metadata)))
  stopifnot("fillBy must be a column in metadata" = fillBy %in% colnames(metadata))
  stopifnot("non-numeric columns in x" = all(apply(x, 2, is.numeric)))

  if (rle) {
    m <- apply(x, 1, median, na.rm = TRUE)
    x <- x - m
  }

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

  if (is.null(fillBy)) {
    fillBy <- ".sample"
  }

  # Force x axis to be in Group order
  fct_levels <- dt.m[order(get(fillBy))][, unique(.sample)]
  dt.m[, .sample := factor(.sample, levels = fct_levels)]

  p <- ggplot2::ggplot(dt.m, ggplot2::aes(x = .data[[".sample"]], y = .data[[".value"]]))
  if (violin) {
    p <- p + ggplot2::geom_violin(ggplot2::aes(fill = .data[[fillBy]]), ...)
  } else {
    p <- p + ggplot2::geom_boxplot(ggplot2::aes(fill = .data[[fillBy]]), ...)
  }

  p + ggplot2::geom_hline(
    yintercept = median(x, na.rm = TRUE),
    color = "red", linetype = 2
  ) +
    ggplot2::labs(
      x = NULL,
      y = NULL,
      fill = if (fillBy == ".sample") "Sample" else fillBy
    )
}

#' @export
#' @rdname plot_boxplot
plot_boxplot.data.frame <- function(x, metadata = NULL, fillBy = NULL,
                                    rle = FALSE, violin = FALSE, ...) {
  if (is(x, "tbl_df") || is(x, "data.table")) {
    stop("You supplied a tibble or a data.table. Please use base::data.frame objects with rownames(x) == colnames(metadata)")
  }
  m <- data.matrix(x)
  plot_boxplot.matrix(m, metadata, fillBy, rle, violin, ...)
}

#' @export
#' @rdname plot_boxplot
plot_boxplot.SummarizedExperiment <- function(x, assay = "counts", fillBy = NULL,
                                              rle = FALSE, violin = FALSE, ...) {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("SummarizedExperiment package is not installed.")
  }
  
  m <- SummarizedExperiment::assay(x, assay)
  metadata <- data.frame(SummarizedExperiment::colData(x))
  plot_boxplot.matrix(m, metadata, fillBy, rle, violin, ...)
}
