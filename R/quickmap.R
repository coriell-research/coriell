#' Heatmap with sensible defaults for RNA-seq expression data
#'
#' Generate a heatmap using pheatmap with sensible defaults for RNA-seq.
#' Utilizes a default blue to red diverging color palette, removes rownames
#' and borders, clusters rows and columns, sets row-wise clustering is euclidean
#' column-wise clustering is correlation with clustering method set to complete.
#' Display font size is reduced and column labels are set to an angle. Defaults can
#' be overridden and/or any additional arguments can be passed to the pheatmap
#' function via ...
#'
#' @param mat numeric matrix to be passed onto pheatmap function
#' @param diverging_palette logical. Default(TRUE). Sets the color scale to a diverging palette (blue -> white -> red). If FALSE, set the
#' color to a continuous color palette (viridis::magma), useful for un-scaled expression data.
#' @param ... args to be passed to pheatmap function
#' @return pheatmap object. See ?pheatmap::pheatmap for details
#' @export
quickmap <- function(mat, diverging_palette = TRUE, ...) {
  diverging_pal <- grDevices::colorRampPalette(c("dodgerblue3", "grey99", "firebrick3"))(50)
  continuous_pal <- rev(viridisLite::magma(n = 50))
  col_code <- if (diverging_palette) diverging_pal else continuous_pal

  default_args <- list(
    scale = "row", show_rownames = FALSE, border_color = NA, cluster_rows = TRUE,
    cluster_cols = TRUE, color = col_code, cellwidth = 20, treeheight_row = 0, clustering_distance_rows = "euclidean",
    clustering_distance_cols = "correlation", clustering_method = "complete", fontsize_col = 8, angle_col = 45
  )
  user_args <- list(...)
  default_args[names(user_args)] <- user_args
  do.call(pheatmap::pheatmap, c(list(mat = mat), default_args))
}
