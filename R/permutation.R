#' Generate all permutations of a vector
#' 
#' Given the input vector, generate a matrix of permutations where each row
#' represents a permutation of the data. Stolen from:
#' \url{https://stackoverflow.com/a/34287541}
#' @param x vector
#' @return factorial(x) x length(x) matrix
#' @examples
#' v <- letters[1:4]
#' p <- permutations(v)
#' p
#' @export
permutations <- function(x) {
  if (length(x) == 1) {
    return(x)
  }
  else {
    res <- matrix(nrow = 0, ncol = length(x))
    for (i in seq_along(x)) {
      res <- rbind(res, cbind(x[i], Recall(x[-i])))
    }
    return(res)
  }
}

#' Perform and exact correlation test on every row of a matrix
#' 
#' Perform an exact correlation permutation test. Computes correlation value for
#' every row of matrix X against all permutations of vector y.
#' @param X numeric matrix or data.frame that can be converted to a numeric matrix
#' @param y numeric vector of data to correlate
#' @param ... arguments passed to the `cor` function
#' @return data.frame containing the original values in X along with columns 
#' containing the correlation of Xi and Y, the empirical p-value of the 
#' permutation test, and the FDR corrected empirical p-value
#' @keywords internal
exact_cor_test <- function(X, y, ...) {
  test_stat <- cor(t(X), y, ...)
  Y_perms <- coriell::permutations(y)
  cors <- apply(Y_perms, 1, cor, x = t(X), ...)
  empirical_p <- vector("numeric", length = nrow(X))
  for (i in seq_along(test_stat)) {
    if (is.na(test_stat[[i]])) {
      empirical_p[[i]] <- NA
    } else if (test_stat[[i]] > 0) {
      empirical_p[[i]] <- mean(cors[i, ] > test_stat[[i]])
    } else if (test_stat[[i]] < 0) {
      empirical_p[[i]] <- mean(cors[i, ] < test_stat[[i]])
    } else {
      empirical_p[[i]] <- 0
    }
  }
  fdr <- p.adjust(empirical_p, method = 'BH')
  df <- as.data.frame(X)
  df$cor <- test_stat
  df$empirical.p <- empirical_p
  df$FDR <- fdr
  df
}

#' Perform a permutation correlation test on every row of a matrix
#' 
#' Compute a correlation value for every row of X against the vector y and n 
#' random permutations of y. If the number of possible permutations is less than
#' the the argument n_perm then an exact test is performed instead. In both 
#' cases the function returns a data.frame of the original data with additional
#' columns for the test statistic, empirical p-value, and FDR corrected empirical
#' p-value.
#' @param X numeric matrix or data.frame that can be converted to a numeric matrix
#' @param y numeric vector of values to correlate with rows of X
#' @param n_perm integer. The desired number of permutations to sample from. Default (10,000)
#' @param n_core integer. The number of cores to use for processing. Default (1)
#' @param ... Additional arguments to pass to `cor` function
#' @export
#' @examples
#' # generate example data
#' X <- matrix(runif(1e3 * 10), nrow = 1e3, ncol = 10)
#' y <- 1:10
#' dimnames(X) <- list(paste("feature", 1:1e3, sep = "."), paste("sample", 1:10, sep = "."))
#' 
#' # correlate each row of X with 1,000 random permutations of vector y
#' res <- permutation_correlation_test(X, y, n_perm = 1e3, n_core = 8, method = "spearman")
#' 
#' head(res)
permutation_correlation_test <- function(X, y, n_perm = 1e4, n_core = 1, ...) {
  X <- if (is.data.frame(X)) as.matrix(X) else X
  if (factorial(length(y)) < n_perm) {
    message(paste0("The number of permutations of y (", factorial(length(y)), ") is less than n_perm (", n_perm, "). Performing exact test on all permutations."))
    df <- coriell::exact_cor_test(X, y, ...)
    return(df)
  }
  
  test_stat <- cor(t(X), y, ...)
  y_perms <- replicate(n_perm, sample(y), simplify = FALSE) 
  perm_cors <- parallel::mclapply(y_perms, FUN = function(y) {cor(t(X), y, ...)}, mc.cores = n_core)
  perm_cors <- do.call(cbind, perm_cors)
  empirical_p <- vector("numeric", length = nrow(X))
  for (i in seq_along(test_stat)) {
    if (is.na(test_stat[[i]])) {
      empirical_p[[i]] <- NA
    } else if (test_stat[[i]] > 0) {
      empirical_p[[i]] <- mean(perm_cors[i, ] > test_stat[[i]])
    } else if (test_stat[[i]] < 0) {
      empirical_p[[i]] <- mean(perm_cors[i, ] < test_stat[[i]])
    } else {
      empirical_p[[i]] <- 0
    }
  }
  fdr <- p.adjust(empirical_p, method = 'BH')
  df <- as.data.frame(X)
  df$cor <- test_stat
  df$empirical.p <- empirical_p
  df$FDR <- fdr
  df
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
  cor_results <- apply(X[random_rows, ],
    MARGIN = 1,
    FUN = function(x) {
      cor(
        as.numeric(x),
        y = sample(y),
        ...
      )
    }
  )
  cor_results
}
