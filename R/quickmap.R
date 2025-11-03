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
#'   * \code{color = colorRampPalette(c("dodgerblue3", "grey99", "firebrick3"))(30)}
#'   * \code{treeheight_row = 0}
#'   * \code{clustering_distance_rows = "correlation"}
#'   * \code{clustering_distance_cols = "euclidean"}
#'   * \code{clustering_method = "ward.D2"}
#'   * \code{angle_col = 315}
#' }
#' @md
#'
#' @param x numeric matrix to be passed onto pheatmap function
#' @param metadata column metadata to add to heatmap annotation. Default NULL
#' @param remove_var Remove this proportion of features based on variance
#' across rows. Default NULL.
#' @param ... args to be passed to \code{pheatmap()} function
#' @return pheatmap object. See \code{?pheatmap::pheatmap()} for details
#' @export
#' @examples
#'
#' # Display heatmap of scaled count data and add title to the plot
#' quickmap(GSE161650_lc, main = "THZ1 vs DMSO")
#'
#' # Override quickmap defaults by supplying additional args
#' quickmap(GSE161650_lc, cluster_rows = FALSE, cluster_cols = FALSE)
#'
quickmap <- function(x, metadata = NULL, remove_var = NULL, ...) {
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("pheatmap package is required.")
  }

  default_args <- list(
    scale = "row",
    show_rownames = FALSE,
    border_color = NA,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    color = grDevices::colorRampPalette(c(
      "dodgerblue3",
      "grey99",
      "firebrick3"
    ))(30),
    treeheight_row = 0,
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "euclidean",
    clustering_method = "ward.D2",
    angle_col = 315
  )
  user_args <- list(...)
  default_args[names(user_args)] <- user_args

  if (!is.null(metadata)) {
    stopifnot(
      "rownames(metadata) != colnames(x)" = rownames(metadata) == colnames(x)
    )
    default_args[["annotation_col"]] <- metadata
  }

  if (!is.null(remove_var)) {
    x <- coriell::remove_var(x, p = remove_var)
  }

  do.call(pheatmap::pheatmap, modifyList(default_args, list(mat = x)))
}
