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
#'   * \code{color = colorRampPalette(c("dodgerblue3", "grey99", "firebrick3"))(n_breaks)} for diverging_palette = TRUE
#'   * \code{color = rev(viridisLite::magma(n = n_breaks))} for diverging_palette = FALSE
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
#' @param n_breaks numeric. The number of breaks to use in the color palette. Default 50.
#' @param fix_extreme logical. Should the extreme values at the ends of a diverging palette be fixed colors? Default FALSE.
#' @param thresh. If fix_extreme is TRUE then what proportion of the extreme values should be filled. For example if the range
#' of the final scaled values is from -20 to 20 and thresh is set to 0.5 then -20 to -10 will be colored blue and 10 to 20
#' will be colored red. Default 0.5. See examples below.
#' @param removeVar If not NULL remove this proportion of features based on the variance across rows. Default NULL. 
#' @param ... args to be passed to \code{pheatmap} function
#' @return pheatmap object. See \code{?pheatmap::pheatmap} for details
#' @examples
#' # generate fake count data
#' X <- coriell::simulate_counts(n_genes = 100)$table
#' 
#' # display heatmap of scaled count data and add title to the plot
#' quickmap(X, main = "Control vs Treatment")
#' 
#' # See effect of fix_extreme -- set thresh high to illustrate effect
#' quickmap(X, main = "Control vs Treatment", fix_extreme = TRUE, thresh = 0.9)
#' @export
quickmap <- function(mat, diverging_palette = TRUE, n_breaks = 50, 
                     fix_extreme = FALSE, thresh = 0.5, removeVar = NULL, ...) {
  stopifnot("thresh must be between 0 and 1" = thresh > 0 & thresh < 1)
  diverging_pal <- grDevices::colorRampPalette(c("dodgerblue3", "grey99", "firebrick3"))(n_breaks)
  continuous_pal <- rev(viridisLite::magma(n = n_breaks))
  col_code <- if (diverging_palette) diverging_pal else continuous_pal

  # Define function for fixing colors at the extreme ends of the color scale
  getBreaks <- function(m, margin, thresh) {
    rng <- range(apply(m, MARGIN = margin, scale))
    ends <- ceiling(max(abs(rng)))
    bk <- seq(-ends, ends, length.out = n_breaks)
    val <- ends - (ends * thresh)
    n_low <- length(bk[bk <= -val])
    n_high <- length(bk[bk >= val])
    n_mid <- length(bk[bk > -val & bk < val])
    pal <- c(
      rep("dodgerblue3", n_low),
      colorRampPalette(c("dodgerblue3", "grey99", "firebrick3"))(n_mid),
      rep("firebrick3", n_high)
    )
    list(breaks = bk, palette = pal)
  }
  
  default_args <- list(
    scale = "row", show_rownames = FALSE, border_color = NA, cluster_rows = TRUE,
    cluster_cols = TRUE, color = col_code, treeheight_row = 0, clustering_distance_rows = "correlation",
    clustering_distance_cols = "euclidean", clustering_method = "ward.D2", angle_col = 315
  )
  user_args <- list(...)
  default_args[names(user_args)] <- user_args
  
  if (!is.null(removeVar)) {
    stopifnot("RemoveVar must be between 0 and 1" = removeVar > 0 & removeVar < 1)
    v <- apply(mat, 1, var)
    o <- order(v, decreasing = TRUE)
    mat <- head(mat[o, ], n = nrow(mat) * (1 - removeVar))
  }
  
  if (fix_extreme) {
    if (default_args[["scale"]] == "none") {
      message("There are no scaled values to fix the extremes (i.e. scale = 'none'). Ignoring")
      return(do.call(pheatmap::pheatmap, c(list(mat = mat), default_args)))
    } else if (default_args[["scale"]] == "row") {
      b <- getBreaks(mat, 1, thresh)
    } else {
      b <- getBreaks(mat, 2, thresh)
    }
    default_args[["color"]] <- b$palette
    default_args[["breaks"]] <- b$breaks
    return(do.call(pheatmap::pheatmap, c(list(mat = mat), default_args)))
  }
  do.call(pheatmap::pheatmap, c(list(mat = mat), default_args))
}
