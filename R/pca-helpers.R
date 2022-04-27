#' Calculate associations between variables and principle components
#' 
#' Calculate associations between metadata variables and PCA rotations. This 
#' function is inspired by \code{methylKit::assocComp} but is designed to work
#' with arbitrary input data. 
#' @details This function returns the p-values from testing the association 
#' between a given metadata column and all PC rotations. For numeric values the
#' p-value returned is computed using the \code{cor.test()} function. For 
#' factor variables and variables that can be converted to factor variables the 
#' function will return p-values from the \code{wilcox.test()} function or the
#' \code{kruskal.test()} function (when the metadata variable has > 2 levels).
#' @param x data.frame or matrix of values passed to \code{prcomp()}
#' @param metadata data.frame of metadata variables to test associations. rownames(metadata) must be identical to colnames(x).
#' @param N First N PCs to test associations between. Default 10.
#' @param ... Additional arguments passed to the \code{prcomp()} function.
#' @export
#' @return data.frame with rows for each metadata variable, columns for each PC, and p-values from the given test in cells.
#' @examples 
#' # Specify metadata
#' metadata <- data.frame(
#'   age = c(30, 80, 34, 30, 80, 40),
#'   treatment =  factor(c(rep("Treatment", 3), rep("Control", 3))),
#'   class = factor(c(rep("A", 2), rep("B", 2), rep("C", 2))),
#'   row.names = c(paste0("trt", 1:3), paste0("ctrl", 1:3))
#' )
#' 
#' # Create values to perform PCA on
#' df <- data.frame(replicate(6, runif(1000, 0, 100)))
#' colnames(df) <- c(paste0("trt", 1:3), paste0("ctrl", 1:3))
#' 
#' # Test for associations
#' res <- associate_components(df, metadata)
#' 
#' # Show results
#' head(res)
associate_components <- function(x, metadata, N = 10, ...) {
  stopifnot("colnames of x and rownames of metadata do no match." = all(colnames(x) == rownames(metadata)))
  stopifnot("No metadata columns present" = ncol(metadata) >= 1)
  stopifnot("NA values detected in metadata" = all(sapply(metadata, anyNA) == FALSE))
  
  # Perform PCA
  res <- prcomp(x, ...)
  
  # Test only the specified number of PCs
  keep_cols <- if (N > ncol(res$rotation)) ncol(res$rotation) else N 
  rotation <- res$rotation[, 1:keep_cols]
  
  # Test the association of each metadata column with each PC rotation
  results <- matrix(data = NA, nrow = ncol(metadata), ncol = ncol(rotation))
  for (i in seq_along(metadata)) {
    v <- metadata[, i]
    for (j in 1:ncol(rotation)) {
      if (is.factor(v) | is.character(v) | is.logical(v)) {
        v <- as.factor(v)
        d <- split(x = rotation[, j], f = v)
        if (nlevels(v) == 2) {
          results[i, j] <- wilcox.test(x = d[[1]], y = d[[2]])$p.value
        } else {
          results[i, j] <- kruskal.test(d)$p.value
        }
      } else {
        results[i, j] <- cor.test(v, rotation[, j])$p.value
      }
    }
  }
  rownames(results) <- colnames(metadata)
  colnames(results) <- colnames(rotation)
  as.data.frame(results)
}
