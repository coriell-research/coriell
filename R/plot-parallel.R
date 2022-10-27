#' Parallel coordinates plot of expression data
#' 
#' The parallel coordinates plot will display a line plot showing the expression
#' value for gene by each sample or summarized expression value of a gene for the 
#' group given by the "colBy" argument on the y-axis and the sample or group 
#' label on the x-axis.
#' 
#' @param x gene by sample matrix or data.frame of numeric values to be plotted.
#' @param metadata data.frame containing metadata per sample. rownames of metadata
#' must match the colnames of the input matrix. Default NULL, each sample in the 
#' matrix will be plotted.
#' @param stat Summary statistic to plot for each gene. One of c("mean", "median", "min", "max").
#' Default "median"
#' @param colBy metadata column used to color lines. Default NULL, every sample 
#' will get its own color.
#' @param removeVar If not NULL remove this proportion of features based on the 
#' variance across rows. Default NULL. 
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
#' # Plot the PCP by Group in metadata
#' plot_parallel(GSE161650_lc, metadata, colBy = "Group", alpha = 0.01) + 
#'   theme_coriell()
plot_parallel <- function(x, metadata = NULL, stat = "median",
                          colBy = NULL, removeVar = NULL, ...) {
  stopifnot("colnames(x) do not match rownames(metadata)" = all(colnames(x) == rownames(metadata)))
  stopifnot("colBy must be a column in metadata" = colBy %in% colnames(metadata))
  stopifnot("non-numeric columns in x" = all(apply(x, 2, is.numeric)))
  
  # Remove low variance features from the input matrix
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
  
  # Summarize the samples/features
  by_group <- dt.m[, .(.mean_val = mean(.value, na.rm = TRUE),
                       .median_val = median(.value, na.rm = TRUE),
                       .min_val = min(.value, na.rm = TRUE),
                       .max_val = max(.value, na.rm = TRUE)), 
                   by = c(".feature", colBy)]
  
  # Determine the summary stat to plot
  stat <- match.arg(stat)
  stat_val <- switch (stat,
    mean = ".mean_val",
    median = ".median_val",
    min = ".min_val",
    max = ".max_val"
  )
  
  # Create a y-axis label based on the summary stat used
  ylabel <- switch (stat,
    mean = "Mean Expr.",
    median = "Median Expr.",
    min = "Min Expr.",
    max = "Max Expr."
  )
  
  ggplot2::ggplot(by_group, ggplot2::aes_string(x = colBy, y = stat_val, group = ".feature")) +
    ggplot2::geom_line(ggplot2::aes_string(color = colBy), ...) +
    ggplot2::labs(x = NULL, y = ylabel) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(alpha = 1)))
}
