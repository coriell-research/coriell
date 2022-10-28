#' Create boxplot from an expression matrix
#' 
#' Create a boxplot of expression values for the given expression matrix. 
#' 
#' @param x gene by sample matrix or data.frame of numeric values to be plotted.
#' @param metadata data.frame containing metadata per sample. rownames of metadata
#' must match the colnames of the input matrix. Default NULL, each sample in the 
#' matrix will be plotted.
#' @param fillBy metadata column used to fill boxplots. Default NULL, 
#' each sample will be a distinct color.
#' @param ... Additional parameters passed to \code{ggplot2::geom_boxplot()}
#' @import data.table
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
plot_boxplot <- function(x, metadata = NULL, fillBy = NULL, ...) {
  stopifnot("colnames(x) do not match rownames(metadata)" = all(colnames(x) == rownames(metadata)))
  stopifnot("fillBy must be a column in metadata" = fillBy %in% colnames(metadata))
  stopifnot("non-numeric columns in x" = all(apply(x, 2, is.numeric)))

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
  
  if (is.null(fillBy)) 
    fillBy <- ".sample"
  
  ggplot2::ggplot(dt.m, ggplot2::aes_string(x = ".sample", y = ".value")) +
    ggplot2::geom_boxplot(ggplot2::aes_string(fill = fillBy), ...) +
    ggplot2::labs(
      x = NULL, 
      y = NULL, 
      fill = if (fillBy == ".sample") "Sample" else fillBy
      )
}
