#' Z-score a dataframe
#' 
#' Convenience function for scaling a dataframe by rows or columns.
#' 
#' @param df data.frame or matrix of numeric columns
#' @param by a character string indicating which axis to scale by. Either "row" or "column". Default ("row")
#' @return data.frame with either
#' @examples
#' df <- data.frame(a = runif(10),
#'                  b = runif(10),
#'                  c = runif(10))
#' rownames(df) <- paste0("row", rownames(df))
#' 
#' zscore_df(df, by = "row")
#' 
#' @export
zscore_df = function(df, by = "row") {
  stopifnot("Non-numeric columns present in df" = all(unlist(lapply(df, is.numeric))) == TRUE)
  stopifnot("by must be either 'row' or 'column'" = by %in% c("row", "column"))
  
  if (by == "row") {
    res <- as.data.frame(t(apply(df, 1, scale)))
  } else {
    res <- as.data.frame(apply(df, 2, scale))
  }
  
  dimnames(res) <- dimnames(df)
  res
} 