#' Compute Pairwise intersection Statistics
#'
#' Given a list of sets, compute all pairwise intersections and return
#' statistics about the intersections.
#' @param sets lists of sets to perform pairwise comparisons on.
#' @return list with elements representing the statistics on intersections of the elements of sets.
#' @details The function returns a single list with several elements: intersection_size, size_A, size_B,
#' size_AB, prop_A, prop_B, prop_AB. Each of the elements contains a named vector with the given statistic.
#' Named vectors represent the comparison of set A to set B from the list of sets given as input.
#' intersection_size is the number of elements from Intersect(A, B). uniq_A is the number of elements unique to
#' set A when comparing set A to set B. uniq_B is the number of elements unique to set B when comparing set A to set B.
#' union_AB is the number of elements in the Union(A, B). prop_A is the number of elements unique to A
#' divided by the Union(A, B). prop_B is the number of elements unique to B divided by the Union(A, B).
#' prop_AB is the number of elements in Intersect(A, B) divided by the number of elements in Union(A, B)
#' @examples
#' sets <- list(
#'   "set1" = letters[1:5],
#'   "set2" = letters[1:10],
#'   "set3" = letters[5:10],
#'   "set4" = letters[2:6]
#' )
#'
#' set_stats <- pairwise_intersection_stats(sets)
#'
#' # access the intersection sizes
#' set_stats$intersection_size
#'
#' # names of the vector represent sets A : B
#' names(set_stats$intersection_size)
#' @export
pairwise_intersection_stats <- function(sets) {
  comb_names <- vector("character", length = length(names(sets)))
  a_names <- vector("character", length = length(names(sets)))
  b_names <- vector("character", length = length(names(sets)))
  lens <- vector("double", length = length(names(sets)))
  pcts <- vector("double", length = length(names(sets)))
  pcts_a <- vector("double", length = length(names(sets)))
  pcts_b <- vector("double", length = length(names(sets)))
  size_a <- vector("double", length = length(names(sets)))
  size_b <- vector("double", length = length(names(sets)))
  size_tot <- vector("double", length = length(names(sets)))

  combinations <- utils::combn(names(sets), m = 2, simplify = FALSE)
  for (i in seq_along(combinations)) {
    c <- combinations[[i]]
    a <- c[1]
    b <- c[2]

    comb_names[[i]] <- paste(a, ":", b)
    int_len <- length(base::intersect(sets[[a]], sets[[b]]))
    union_len <- length(base::union(sets[[a]], sets[[b]]))
    uniq_A <- length(base::setdiff(sets[[a]], sets[[b]]))
    uniq_B <- length(base::setdiff(sets[[b]], sets[[a]]))

    pct_tot <- int_len / union_len
    pct_a <- uniq_A / union_len
    pct_b <- uniq_B / union_len

    lens[[i]] <- int_len
    pcts[[i]] <- pct_tot
    pcts_a[[i]] <- pct_a
    pcts_b[[i]] <- pct_b
    size_a[[i]] <- length(uniq_A)
    size_b[[i]] <- length(uniq_B)
    size_tot[[i]] <- union_len
  }

  names(lens) <- comb_names
  names(pcts) <- comb_names
  names(pcts_a) <- comb_names
  names(pcts_b) <- comb_names
  names(size_a) <- comb_names
  names(size_b) <- comb_names
  names(size_tot) <- comb_names

  res <- list(
    "intersection_size" = lens,
    "uniq_A" = size_a,
    "uniq_B" = size_b,
    "union_AB" = size_tot,
    "prop_A" = pcts_a,
    "prop_B" = pcts_b,
    "prop_AB" = pcts
  )

  res
}
