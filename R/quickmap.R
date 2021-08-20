#' Heatmap with sensible defaults for RNA-seq expression data
#'
#' Generate a heatmap using \code{pheatmap} with sensible defaults for RNA-seq.
#' 
#' @details The default arguments to \code{pheatmap::pheatmap()} are:
#' \itemize{
#'   * \code{scale = "row"}
#'   * \code{show_rownames = FALSE}
#'   * \code{border_col = NA}
#'   * \code{cluster_rows = TRUE}
#'   * \code{cluster_cols = TRUE}
#'   * \code{color = colorRampPalette(c("dodgerblue3", "grey99", "firebrick3"))(50)} for diverging_palette = TRUE
#'   * \code{color = rev(viridisLite::magma(n = 50))} for diverging_palette = FALSE
#'   * \code{treeheight_row = 0}
#'   * \code{clustering_distance_rows = "correlation"}
#'   * \code{clustering_distance_cols = "euclidean"}
#'   * \code{clustering_method = "ward.D2"}
#'   * \code{angle_col = 315}
#' }
#' You can pass in additional arguments or simply override the defaults as well.
#'
#' @md
#' @param mat numeric matrix to be passed onto pheatmap function
#' @param diverging_palette logical. Default(TRUE). Sets the color scale to a diverging palette (blue -> white -> red). If FALSE, set the
#' color to a continuous color palette \code{viridis::magma}, useful for un-scaled expression data.
#' @param ... args to be passed to \code{pheatmap} function
#' @return pheatmap object. See \code{?pheatmap::pheatmap} for details
#' @examples
#' # generate fake count data
#' X <- coriell::simulate_counts(n_genes = 100)$table
#' 
#' # display heatmap of scaled count data and add title to the plot
#' quickmap(X, main = "Control vs Treatment")
#' @export
quickmap <- function(mat, diverging_palette = TRUE, ...) {
  diverging_pal <- grDevices::colorRampPalette(c("dodgerblue3", "grey99", "firebrick3"))(50)
  continuous_pal <- rev(viridisLite::magma(n = 50))
  col_code <- if (diverging_palette) diverging_pal else continuous_pal

  default_args <- list(
    scale = "row", show_rownames = FALSE, border_color = NA, cluster_rows = TRUE,
    cluster_cols = TRUE, color = col_code, treeheight_row = 0, clustering_distance_rows = "correlation",
    clustering_distance_cols = "euclidean", clustering_method = "ward.D2", angle_col = 315
  )
  user_args <- list(...)
  default_args[names(user_args)] <- user_args
  do.call(pheatmap::pheatmap, c(list(mat = mat), default_args))
}
