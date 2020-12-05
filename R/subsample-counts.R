#' Sub-sample a count matrix to the smallest library size
#'
#' Randomly sub-sample a count matrix so that the column sums are equivalent
#' to the library size of the smallest sample (column). Returns a numeric matrix
#' of subsampled counts with the same dimensions as the input matrix.
#'
#' @param counts numeric matrix of integer counts.
#' @param seed integer. Seed for random number generator. Default (123)
#' @return numeric matrix of sub-sampled counts.
#' @export
#' @examples
#' # generate simulated count data
#' mat <- simulate_counts()$table
#' colSums(mat)
#'
#' # subsample the count matrix
#' ss <- subsample_counts(mat)
#' colSums(ss)
subsample_counts <- function(counts, seed = 123) {
  stopifnot("Non-numeric columns in matrix" = all(apply(counts, 2, is.numeric) == TRUE))
  stopifnot("colnames must be non-NULL" = length(colnames(counts)) == ncol(counts))
  stopifnot("Seed must be numeric" = is.numeric(seed))

  set.seed(seed)
  target_depth <- min(colSums(counts))
  target_cols <- which(colSums(counts) > target_depth)
  skip_cols <- which(!1:ncol(counts) %in% target_cols)
  features <- rownames(counts)
  subsamples <- vector("list", length(target_cols))

  # generate subsample dataframes
  for (i in target_cols) {
    sample_id <- colnames(counts)[i]
    pool <- rep(features, counts[, i, drop = TRUE])
    subsample <- sample(pool, size = target_depth, replace = FALSE)
    subsample_df <- as.data.frame(table(feature_id = subsample), responseName = sample_id)
    subsamples[[i]] <- subsample_df
  }

  # add skip columns as dataframes to subsamples dataframe list
  for (i in skip_cols) {
    sample_id <- colnames(counts)[i]
    d <- data.frame(feature_id = features, sample_id = counts[, i])
    names(d)[names(d) == "sample_id"] <- sample_id
    subsamples[[i]] <- d
  }

  # combine all dfs, convert columns to rownames and replace NAs with 0s
  df <- Reduce(function(df1, df2) {
    merge(df1, df2, by = "feature_id", all = TRUE)
  }, subsamples)

  # preserve original ordering
  df <- df[match(features, df$feature_id), ]

  # convert feature_id col to rownames
  rownames(df) <- df[, 1]
  df[, 1] <- NULL

  # convert to matrix and replace NAs
  mat <- as.matrix(df)
  mat[is.na(mat)] <- 0
  mat
}
