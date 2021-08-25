#' Convert a list of vectors to a binary matrix
#'
#' Convert a list of sets into a binary matrix.
#'
#' @param sets named list of sets (vectors) to be converted to binary matrix.
#' @return binary matrix with a column for each set. rownames of the matrix represent the union of all sets.
#' A '1' indicates the inclusion of the element in the set.
#' @export
#' @examples
#' sets <- list(
#'   "set1" = letters[1:5],
#'   "set2" = letters[2:6],
#'   "set3" = letters[1:7]
#' )
#' 
#' list_to_matrix(sets)
list_to_matrix <- function(sets) {
  stopifnot("List of vectors must be supplied" = is(sets, "list"))
  union_all <- Reduce(union, sets)

  if (sum(is.na(union_all)) > 0) {
    message("NA values present in union of all sets. NA values will be dropped in final matrix")
    union_all <- union_all[!is.na(union_all)]
  }

  mat <- matrix(
    data = 0,
    nrow = length(union_all),
    ncol = length(sets)
  )
  colnames(mat) <- names(sets)
  rownames(mat) <- union_all

  for (i in seq_along(sets)) {
    mat[unique(sets[[i]][!is.na(sets[[i]])]), i] <- 1
  }
  mat
}
