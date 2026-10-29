#' Generate all permutations of a vector
#'
#' Given the input vector, generate a matrix of permutations where each row
#' represents a permutation of the data. Permutations are returned in
#' lexicographic order of the element positions.
#' @param x vector
#' @return factorial(length(x)) x length(x) matrix
#' @examples
#' v <- letters[1:4]
#' p <- permutations(v)
#' p
#' @export
permutations <- function(x) {
  n <- length(x)
  if (n == 1) {
    return(x)
  }
  # Build index permutations of 1:k from those of 1:(k-1) by choosing the
  # first element and shifting the remaining indices up
  idx <- matrix(1L, nrow = 1, ncol = 1)
  for (k in 2:n) {
    idx <- do.call(rbind, lapply(seq_len(k), function(i) {
      cbind(i, idx + (idx >= i))
    }))
  }
  matrix(x[idx], nrow = nrow(idx), ncol = n)
}

#' Perform an exact correlation test on every row of a matrix
#'
#' Perform an exact correlation permutation test. Computes correlation value for
#' every row of matrix X against all distinct permutations of vector y.
#' @param X numeric matrix or data.frame that can be converted to a numeric matrix
#' @param y numeric vector of data to correlate
#' @param method character. One of "pearson" (default) or "spearman"
#' @param alternative character. One of "two.sided" (default), "greater" or "less"
#' @return data.frame containing the original values in X along with columns
#' containing the correlation of Xi and Y, the exact permutation p-value, and
#' the FDR corrected p-value
#' @keywords internal
#' @export
exact_cor_test <- function(X, y, method = c("pearson", "spearman"),
                           alternative = c("two.sided", "greater", "less")) {
  perm_cor_test(
    X, y,
    n_perm = Inf,
    method = match.arg(method),
    alternative = match.arg(alternative)
  )
}

#' Perform a permutation correlation test on every row of a matrix
#'
#' Compute a correlation value for every row of X against the vector y and
#' n_perm random permutations of y. If the number of distinct permutations of y
#' is less than or equal to n_perm then an exact test using all permutations is
#' performed instead. In both cases the function returns a data.frame of the
#' original data with additional columns for the test statistic, permutation
#' p-value, and FDR corrected p-value.
#'
#' Exact p-values are computed as the proportion of all distinct permutations
#' (including the observed ordering) with a correlation at least as extreme as
#' the observed value. When permutations are randomly sampled, p-values are
#' computed with \code{statmod::permp()} (if installed), which gives exact
#' p-values for permutations sampled with replacement (Phipson & Smyth 2010).
#' Otherwise the conservative estimate (b + 1) / (n_perm + 1) is used.
#'
#' Correlations are computed by standardizing the rows of X once and performing
#' matrix multiplication against chunks of permuted y. Rows of X containing
#' missing values or with zero variance return NA.
#' @param X numeric matrix or data.frame that can be converted to a numeric matrix
#' @param y numeric vector of values to correlate with rows of X
#' @param n_perm integer. The desired number of permutations to sample from. Default (10,000)
#' @param method character. One of "pearson" (default) or "spearman"
#' @param alternative character. One of "two.sided" (default), "greater" or "less"
#' @param chunk_size integer. Number of permutations to process at once. Controls
#' memory usage, nrow(X) x chunk_size correlations are held in memory at a time. Default (1000)
#' @param n_core integer. The number of cores to use for processing chunks of
#' permutations. Default (1)
#' @return data.frame containing the original values in X along with columns
#' 'cor', 'empirical.p', and 'FDR'
#' @references Phipson B, Smyth GK (2010). Permutation P-values should never be
#' zero: calculating exact P-values when permutations are randomly drawn.
#' Statistical Applications in Genetics and Molecular Biology 9, Article 39.
#' @export
#' @examples
#' # generate example data
#' X <- matrix(runif(1e3 * 10), nrow = 1e3, ncol = 10)
#' y <- 1:10
#' dimnames(X) <- list(paste("feature", 1:1e3, sep = "."), paste("sample", 1:10, sep = "."))
#'
#' # correlate each row of X with 1,000 random permutations of vector y
#' res <- permutation_correlation_test(X, y, n_perm = 1e3, method = "spearman")
#'
#' head(res)
permutation_correlation_test <- function(X, y, n_perm = 1e4,
                                         method = c("pearson", "spearman"),
                                         alternative = c("two.sided", "greater", "less"),
                                         chunk_size = 1e3, n_core = 1) {
  perm_cor_test(
    X, y,
    n_perm = n_perm,
    method = match.arg(method),
    alternative = match.arg(alternative),
    chunk_size = chunk_size,
    n_core = n_core
  )
}

# Shared implementation for exact and sampled permutation correlation tests
perm_cor_test <- function(X, y, n_perm, method, alternative,
                          chunk_size = 1e3, n_core = 1) {
  X <- if (is.data.frame(X)) as.matrix(X) else X
  if (!is.matrix(X) || !is.numeric(X)) {
    stop("X must be a numeric matrix or data.frame")
  }
  if (!is.numeric(y) || anyNA(y)) {
    stop("y must be a numeric vector without missing values")
  }
  if (ncol(X) != length(y)) {
    stop("ncol(X) must equal length(y)")
  }

  n <- length(y)
  Xs <- standardize_rows(X, method)
  ys <- standardize_rows(matrix(y, nrow = 1), method)[1, ]
  test_stat <- drop(Xs %*% ys)
  test_stat[!is.finite(test_stat)] <- NA

  # Number of distinct permutations of y (accounting for exact ties)
  yi <- match(y, unique(y))
  total_nperm <- round(exp(lfactorial(n) - sum(lfactorial(tabulate(yi)))))
  tol <- sqrt(.Machine$double.eps)

  if (total_nperm <= n_perm) {
    if (is.finite(n_perm)) {
      message(paste0(
        "The number of distinct permutations of y (", total_nperm,
        ") is less than or equal to n_perm (", n_perm,
        "). Performing exact test on all permutations."
      ))
    }
    idx <- permutations(seq_len(n))
    idx <- idx[!duplicated(matrix(yi[idx], ncol = n)), , drop = FALSE]
    Yp <- matrix(ys[t(idx)], nrow = n)
    count <- count_extreme(Xs %*% Yp, test_stat, alternative, tol)
    empirical_p <- count / total_nperm
  } else {
    # Generate all permutation indices up front so results are reproducible
    # with set.seed() regardless of n_core
    idx <- vapply(seq_len(n_perm), function(i) sample.int(n), integer(n))
    chunks <- split(seq_len(n_perm), ceiling(seq_len(n_perm) / chunk_size))
    counts <- parallel::mclapply(
      chunks,
      FUN = function(j) {
        Yp <- matrix(ys[idx[, j, drop = FALSE]], nrow = n)
        count_extreme(Xs %*% Yp, test_stat, alternative, tol)
      },
      mc.cores = n_core
    )
    count <- Reduce(`+`, counts)
    empirical_p <- sampled_pvalue(count, n_perm, total_nperm)
  }

  df <- as.data.frame(X)
  df$cor <- test_stat
  df$empirical.p <- empirical_p
  df$FDR <- p.adjust(empirical_p, method = "BH")
  df
}

# Center and scale rows to unit norm so that tcrossprod gives correlations.
# Spearman correlation is Pearson correlation of (average) ranks
standardize_rows <- function(X, method) {
  if (method == "spearman") {
    X <- if (requireNamespace("matrixStats", quietly = TRUE)) {
      matrixStats::rowRanks(X, ties.method = "average")
    } else {
      t(apply(X, 1, rank, na.last = "keep"))
    }
    X <- matrix(as.numeric(X), nrow = nrow(X))
  }
  X <- X - rowMeans(X)
  X / sqrt(rowSums(X^2))
}

# Count the permuted correlations at least as extreme as the observed values
count_extreme <- function(R, test_stat, alternative, tol) {
  switch(alternative,
    two.sided = rowSums(abs(R) >= abs(test_stat) - tol),
    greater = rowSums(R >= test_stat - tol),
    less = rowSums(R <= test_stat + tol)
  )
}

# P-values for permutations sampled with replacement
sampled_pvalue <- function(count, n_perm, total_nperm) {
  p <- rep(NA_real_, length(count))
  ok <- !is.na(count)
  u <- unique(count[ok])
  pu <- if (requireNamespace("statmod", quietly = TRUE)) {
    # permp() can build a total_nperm x length(x) matrix so evaluate in chunks
    unlist(lapply(split(u, ceiling(seq_along(u) / 100)), function(x) {
      statmod::permp(x, nperm = n_perm, total.nperm = total_nperm)
    }), use.names = FALSE)
  } else {
    (u + 1) / (n_perm + 1)
  }
  p[ok] <- pu[match(count[ok], u)]
  p
}

#' Generate a null distribution of correlation values
#'
#' Selects n random rows from a numeric matrix with replacement. For each random
#' row, permute vector y and perform correlation.
#' @param X numeric matrix or data.frame that can be converted to numeric matrix.
#' @param y numeric vector. Numeric vector of values used to correlate with each row of df
#' @param n integer. Number of random correlations to return. Default (1e4)
#' @param ... Additional arguments passed to `cor` function
#' @return numeric vector of correlation values. Vector can contain NAs if row variance == 0.
#' @examples
#' # generate example data
#' X <- matrix(runif(1e3 * 6), nrow = 1e3, ncol = 6)
#' y <- 1:6
#' dimnames(X) <- list(paste("feature", 1:1e3, sep = "."), paste("sample", 1:6, sep = "."))
#'
#' # sample random correlations from permuted data
#' null_dist <- sample_n_random_cor(X, y, n = 1e3, method = "spearman")
#'
#' head(null_dist)
#' @export
sample_n_random_cor <- function(X, y, n = 1e4, ...) {
  X <- if (is.data.frame(X)) as.matrix(X) else X
  random_rows <- sample.int(nrow(X), n, replace = TRUE)
  cor_results <- apply(X[random_rows, , drop = FALSE], MARGIN = 1, FUN = function(x) {
    cor(
      as.numeric(x),
      y = sample(y),
      ...
    )
  })
  cor_results
}
