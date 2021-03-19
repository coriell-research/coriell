#' Z-score a dataframe
#'
#' Convenience function for scaling a dataframe by rows or columns.
#' @param df data.frame or matrix of numeric columns
#' @param by a character string indicating which axis to scale by. Either "row" or "column". Default ("row")
#' @param ... Additional arguments passed to `scale`
#' @return data.frame with scaled values
#' @examples
#' # generate example data
#' df <- data.frame(
#'   a = runif(10),
#'   b = runif(10),
#'   c = runif(10)
#' )
#' rownames(df) <- paste0("row", rownames(df))
#'
#' # z-score each row of the data.frame
#' zscore_df(df, by = "row")
#' @export
zscore_df <- function(df, by = c("row", "column"), ...) {
  stopifnot("Non-numeric columns present in df" = all(unlist(lapply(df, is.numeric))) == TRUE)
  
  X <- as.matrix(df)
  d <- match.arg(by)
  res <- switch(d,
    row = t(scale(t(X), ...)),
    column = scale(X, ...)
  )
  as.data.frame(res)
}
