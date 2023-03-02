#' Compute outliers by IQR method
#'
#' Compute outliers for all columns of a numeric matrix using the IQR method. 
#' For each column of the input matrix, a value is called as an outlier if the 
#' value is less than the first quartile i.e. Q1 - IQR * scale factor or greater 
#' than Q3 + IQR * scale factor.
#'
#' @param X numeric matrix or data.frame that can be converted to a numeric 
#' matrix with variables in the columns and sample names in the rows.
#' @param scale_factor numeric. Factor to scale the outlier range. default 1.5.
#' @return numeric matrix where values of 1L indicates an outlier and 0L indicates 
#' non-outlier
#' @export
#' @examples
#' set.seed(12345)
#' M <- matrix(
#'   data = c(rnorm(10, 10, 1), rnorm(10, 100, 15)), 
#'   ncol = 2, 
#'   dimnames = list(paste0("sample", 1:10), c("var1", "var2"))
#'   )
#' 
#' # Create one outlier in first and last rows
#' M[1, 1] <- 100
#' M[1, 2] <- 1000
#' M[10, 1] <- -10
#' M[10, 2] <- 0
#' 
#' # Show outliers on boxplot
#' boxplot(M)
#' 
#' # Call outliers in each column
#' outliers_by_iqr(M)
#' 
outliers_by_iqr <- function(X, scale_factor = 1.5) {
  if (!requireNamespace("matrixStats", quietly = TRUE)) {
    stop("matrixStats package is required.")
  }
  
  if (is(X, "data.frame")) X <- data.matrix(X)
  iqr <- matrixStats::colIQRs(X, useNames = FALSE)
  q1 <- matrixStats::colQuantiles(X, probs = 0.25, useNames = FALSE)
  q3 <- matrixStats::colQuantiles(X, probs = 0.75, useNames = FALSE)
  lower_thresh <- q1 - (iqr * scale_factor)
  upper_thresh <- q3 + (iqr * scale_factor)
  outlier_mat <- t(ifelse(t(X) > upper_thresh | t(X) < lower_thresh, 1L, 0L))
  dimnames(outlier_mat) <- dimnames(X)

  return(outlier_mat)
}