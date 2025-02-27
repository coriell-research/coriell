#' Plot the distance between all columns of a matrix
#'
#' Create a column-vs-column heatmap of a matrix
#' @param x numeric matrix or data.frame that can be converted to one. Samples
#' in columns, features in rows. 
#' @param metadata data.frame containing metadata to be used as row labels.
#' Default NULL.
#' @param method distance method to be used. Must be one of the distance
#' measure to be used. This must be one of "euclidean", "maximum", "manhattan",
#' "canberra", "binary" or "minkowski"
#' @param ... additional values passed to \code{pheatmap::pheatmap()}
#' @export
#' @return heatmap with distance values displayed in cells.
#' @examples
#' # Generate a matrix of simulated counts
#' counts <- simulate_counts()$table
#' 
#' # Create annotation metadata
#' df <- data.frame(
#'   row.names = colnames(counts), 
#'   Group = rep(c("Control", "Treatment"), each = 3)
#'   )
#' 
#' # Show the sample-vs-sample distances
#' plot_dist(counts)
#' 
#' # Additional arguments can be passed to the function
#' plot_dist(
#'   x = counts, 
#'   metadata = df, 
#'   main = "Sample-vs-Sample Distance", 
#'   color = viridisLite::rocket(n = 50),
#'   annotation_colors = list(Group = c("Treatment" = "firebrick2", "Control" = "grey80"))
#'   )
plot_dist <- function(x, metadata = NULL, method = "euclidean", ...) {
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("pheatmap package is required.")
  }

  if (!requireNamespace("viridisLite", quietly = TRUE)) {
    stop("viridisLite package is required.")
  }

  default_args <- list(
    scale = "none",
    show_rownames = TRUE,
    show_colnames = TRUE,
    angle_col = 315,
    border_color = "grey20",
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    treeheight_row = 0,
    display_numbers = TRUE,
    number_format = "%.1f",
    number_color = "black",
    color = viridisLite::magma(n = 50)
  )
  user_args <- list(...)
  default_args[names(user_args)] <- user_args

  if (requireNamespace("rdist", quietly = TRUE)) {
    d <- rdist::rdist(t(x), metric = method)
  } else {
    d <- dist(t(x), method = method)
  }

  # Create annotations if metadata is present
  if (!is.null(metadata)) {
    stopifnot("rownames(metadata) != colnames(x)" = rownames(metadata) == colnames(x))
    default_args[["annotation_row"]] <- metadata
    default_args[["annotation_col"]] <- metadata
  }

  default_args[["mat"]] <- as.matrix(d)
  do.call(pheatmap::pheatmap, default_args)
}
