#' Density plot of expression values
#' 
#' Create a density of expression values for the given expression matrix.
#' 
#' @param x gene by sample matrix or data.frame of numeric values to be plotted.
#' @param metadata data.frame containing metadata per sample. rownames of metadata
#' must match the colnames of the input matrix. Default NULL, each sample in the 
#' matrix will be plotted.
#' @param colBy metadata column used to color density lines. Default NULL, 
#' each sample in the matrix will be plotted.
#' @param ... Additional parameters passed to \code{ggplot2::geom_density()}
#' @import data.table
#' @return ggplot object
#' @export
#' @examples
#' # Create metadata for plotting
#' metadata <- data.frame(row.names = colnames(GSE161650_lc))
#' metadata$Group <- rep(c("DMSO", "THZ1"), each = 3)
#' 
#' # Plot the density by sample
#' plot_density(GSE161650_lc) + 
#'   theme_coriell()
#'   
#' # Color each sample by their Group in metadata
#' plot_density(GSE161650_lc, metadata, colBy = "Group") + 
#'   theme_coriell()
plot_density <- function(x, metadata = NULL, colBy = NULL, ...) {
  stopifnot("colnames(x) do not match rownames(metadata)" = all(colnames(x) == rownames(metadata)))
  stopifnot("colBy must be a column in metadata" = colBy %in% colnames(metadata))
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
  
  if (is.null(colBy)) 
    colBy <- ".sample"
  
  ggplot2::ggplot(dt.m, ggplot2::aes_string(x = ".value", group = ".sample")) +
    ggplot2::geom_density(ggplot2::aes_string(color = colBy), ...) +
    ggplot2::labs(
      x = NULL, 
      y = "Density", 
      color = if (colBy == ".sample") "Sample" else colBy
      )
}