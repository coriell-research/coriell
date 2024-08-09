# Fix colors at the extreme ends of the color scale
#
.getBreaks <- function(m, thresh, n_breaks=30) {
  rng <- range(m[!is.na(m)])
  
  if (-thresh < rng[1] || thresh > rng[2]) {
    thresh <- floor(min(abs(rng)))
    msg <- paste0("thresh must be within the range of the matrix. low=", 
                  round(rng[1], 2), " high=", round(rng[2], 2), 
                  ". Setting threshold to +/- ", thresh)
    warning(msg)
  }
  
  bk <- seq(rng[1], rng[2], length.out = n_breaks)
  n_low <- sum(bk < -thresh)
  n_high <- sum(bk > thresh)
  n_mid <- n_breaks - (n_low + n_high)

  pal <- c(
    rep("dodgerblue3", n_low),
    colorRampPalette(c("dodgerblue3", "grey99", "firebrick3"))(n_mid),
    rep("firebrick3", n_high)
  )
  
  list(breaks = bk, palette = pal)
}

# Cluster a matrix using fastcluster, Rfast, and rdist
#
# Overide pheatmap's default function to use faster methods if available
.clusterMatrix <- function(x, distance, method) {
  if (distance == "correlation") {
    if (requireNamespace("Rfast", quietly = TRUE)) {
      d <- as.dist(1 - Rfast::cora(t(x), large = TRUE))
    } else {
      d <- as.dist(1 - cor(t(x)))
    }
  } else {
    if (requireNamespace("rdist", quietly = TRUE)) {
      d <- rdist::rdist(x, metric = distance)
    } else {
      d <- dist(x, method = distance)
    }
  }

  if (requireNamespace("fastcluster", quietly = TRUE)) {
    return(fastcluster::hclust(d, method = method))
  }
  hclust(d, method = method)
}

#' Heatmap with sensible defaults for RNA-seq expression data
#'
#' Generate a heatmap using \code{pheatmap} with sensible defaults for RNA-seq.
#' \code{quickmap()} will also attempt to perform vectorized scaling, clustering,
#' and distance calculations with functions from the \code{Rfast},
#' \code{fastcluster} and \code{rdist} packages in order to speed up
#' calculations for large gene expression matrices.
#'
#' @details The default arguments to \code{pheatmap::pheatmap()} are:
#' \itemize{
#'   * \code{scale = "row"}
#'   * \code{show_rownames = FALSE}
#'   * \code{border_col = NA}
#'   * \code{cluster_rows = TRUE}
#'   * \code{cluster_cols = TRUE}
#'   * \code{color = colorRampPalette(c("dodgerblue3", "grey99", "firebrick3"))(30)} for diverging_palette = TRUE
#'   * \code{color = rev(viridisLite::magma(n = 30))} for diverging_palette = FALSE
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
#' @param diverging_palette logical. Default(TRUE). Sets the color scale to a 
#' diverging palette (blue -> white -> red). If FALSE, set the color to a 
#' continuous color palette \code{viridis::magma()}, useful for un-scaled 
#' expression data.
#' @param fix_extreme logical. Should the extreme values at the ends of a 
#' diverging palette be fixed colors? Default FALSE.
#' @param thresh If fix_extreme=TRUE then at what value should the color scale 
#' be fixed? Since scaling is implied by the use of this argument, the default
#' value is set to 1.96 (positive and negative).
#' @param removeVar If not NULL, remove this proportion of features based on 
#' the variance across rows. Default NULL. Note, if removeVar is used the 
#' lowest variance features are removed prior to scaling the input data
#' @param ... args to be passed to \code{pheatmap()} function
#' @return pheatmap object. See \code{?pheatmap::pheatmap()} for details
#' @export
#' @examples
#' # display heatmap of scaled count data and add title to the plot
#' quickmap(GSE161650_lc, main = "THZ1 vs DMSO")
#'
#' # Remove 90% lowest variance features and fix color scale
#' quickmap(
#'   GSE161650_lc,
#'   removeVar = 0.9,
#'   main = "THZ1 vs DMSO",
#'   fix_extreme = TRUE,
#'   thresh = 1
#' )
quickmap <- function(mat, diverging_palette = TRUE, fix_extreme = FALSE, 
                     thresh = 1.96, removeVar = NULL, ...) {
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("pheatmap package is required.")
  }

  if (!requireNamespace("viridisLite", quietly = TRUE)) {
    stop("viridisLite package is required.")
  }

  diverging_pal <- grDevices::colorRampPalette(c("dodgerblue3", "grey99", "firebrick3"))(30)
  continuous_pal <- rev(viridisLite::magma(n = 30))
  col_code <- if (diverging_palette) diverging_pal else continuous_pal

  default_args <- list(
    scale = "row",
    show_rownames = FALSE,
    border_color = NA,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    color = col_code,
    treeheight_row = 0,
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "euclidean",
    clustering_method = "ward.D2",
    angle_col = 315
  )
  user_args <- list(...)
  default_args[names(user_args)] <- user_args

  if (!is.null(removeVar)) 
    mat <- remove_var(mat, p = removeVar)

  # Override pheatmap scaling, uses slower apply call internally
  if (default_args[["scale"]] == "row") {
    mat <- t(scale(t(mat), center = TRUE, scale = TRUE))
    default_args[["scale"]] <- "none"
  } else if (default_args[["scale"]] == "column") {
    mat <- scale(mat, center = TRUE, scale = TRUE)
    default_args[["scale"]] <- "none"
  }
  
  # Fix extreme values at ends of color scale
  if (isTRUE(fix_extreme)) {
    b <- .getBreaks(mat, thresh)
    default_args[["color"]] <- b$palette
    default_args[["breaks"]] <- b$breaks
  }
  
  # Calculate distances and clustering using faster functions if available
  drow <- default_args[["clustering_distance_rows"]]
  dcol <- default_args[["clustering_distance_cols"]]
  meth <- default_args[["clustering_method"]]

  if (isTRUE(default_args[["cluster_rows"]])) {
    hc_row <- .clusterMatrix(mat, distance = drow, method = meth)
    default_args[["cluster_rows"]] <- hc_row
  }

  if (isTRUE(default_args[["cluster_cols"]])) {
    hc_col <- .clusterMatrix(t(mat), distance = dcol, method = meth)
    default_args[["cluster_cols"]] <- hc_col
  }
  
  do.call(pheatmap::pheatmap, modifyList(default_args, list(mat = mat)))
}
