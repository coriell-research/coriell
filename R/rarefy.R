#' Rarefy (subsample) a matrix
#'
#' This function will randomly subsample counts from rows of the input matrix
#' such that the colSums all have even depth.
#' @param x numeric matrix or data.frame that can be converted to a numeric
#' matrix. Samples in columns features in rows.
#' @param depth desired sampling depth applied to each library.
#' @export
#' @examples
#' set.seed(123)
#' m <- matrix(
#'   sample.int(100),
#'   nrow = 10,
#'   dimnames = list(
#'     Gene = paste0("gene.", 1:10),
#'     Sample = paste0("sample.", 1:10)
#'   )
#' )
#'
#' colSums(m)
#'
#' rarefied <- rarefy(m, depth = 100)
#' colSums(rarefied)
rarefy <- function(x, depth) {
  stopifnot("NAs present in input matrix" = !any(apply(x, 2, anyNA)))
  stopifnot(
    "depth is greater than smallest library size" = !any(
      depth > min(colSums(x))
    )
  )
  subsample <- function(x, d) {
    n <- length(x)
    idx <- 1:n
    pool <- rep(idx, x)
    rarefied <- sample(pool, size = d)
    tabulate(rarefied, nbins = n)
  }
  m <- apply(x, 2, subsample, d = depth)
  dimnames(m) <- list(rownames(x), colnames(x))
  m
}
