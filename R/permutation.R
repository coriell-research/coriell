#' Correlation Permutation Test Function
#'
#' Internal function for performing a correlation permutation test
#' @param x numeric vector
#' @param y numeric vector
#' @param n_perm numeric. The number of permutations to perform
#' @param cor_method a character string indicating which correlation coefficient (or covariance) is to be computed.
#'   One of "pearson", "kendall", or "spearman": can be abbreviated.
#' @importFrom stats cor
#' @return numeric vector of length n_perm
perm_cor <- function(x, y, n_perm, cor_method) {
  # calculate the actual value
  test_stat <- cor(x, y, method = cor_method)

  # test_stat will be NA for rows without variance
  if (is.na(test_stat)) {
    return(NA)
  }

  # Perform the permutation test
  cor_results <- vector("double", length = n_perm)
  for (i in 1:n_perm) {
    perm_vec <- sample(x)
    cor_results[[i]] <- cor(perm_vec, y, method = cor_method)
  }

  # calculate and return the empirical p-value
  if (test_stat < 0) {
    mean(cor_results <= test_stat)
  } else {
    mean(cor_results >= test_stat)
  }
}

#' Parallel Implementation of a Correlation Permutation Test
#'
#' Performs a correlation permutation test over all rows of a data.frame using multiple cores
#' @param df data.frame. Site by Samples data.frame
#' @param y numeric vector. Numeric vector of values used to correlate with each row of df
#' @param n_cores numeric. Number of cores to use for parallel processing.  Default 1.
#' @param n_perm numeric. Number of permutations to perform on each row. Default 1000.
#' @param cor_method a character string indicating which correlation coefficient (or covariance) is to be computed.
#'   One of "pearson", "kendall", or "spearman" (default): can be abbreviated.
#' @param p_adjust_method a character string indicating which of the p.adjust.methods to use for correction. Default "fdr". 
#' @return df with additional columns for correlations, empirical p-values, and fdr adjusted p-values.
#' @importFrom stats cor p.adjust
#' @importFrom foreach %dopar%
#' @importFrom parallel detectCores
#' @examples
#' \dontrun{
#' library(methylKit)
#'
#' ages <- data.frame(age = c(30, 80, 34, 30, 80, 40))
#'
#' sim_meth <- dataSim(
#'   replicates = 6,
#'   sites = 1000,
#'   treatment = c(rep(1, 3), rep(0, 3)),
#'   covariates = ages,
#'   sample.ids = c(paste0("test", 1:3), paste0("ctrl", 1:3))
#' )
#'
#' # extract the methylation as percentages and coerce to data.frame
#' perc_meth <- as.data.frame(percMethylation(sim_meth))
#'
#' # perform permutation testing using 4 cores and 1000 permutations
#' # return a dataframe (res) with the new columns
#' res <- permutation_correlation_test(perc_meth, y = ages$age, n_cores = 4, n_perm = 1000)
#' }
#' @export
permutation_correlation_test <- function(df, y, n_cores = 1, n_perm = 1000, cor_method = "spearman", p_adjust_method = "fdr") {
  # Simple error checking
  stopifnot("Number of columns is not equal to length of vector y" = length(colnames(df)) == length(y))
  stopifnot("n_cores exceeds number of cores detected" = n_cores <= parallel::detectCores())
  stopifnot("Incorrect p.adjust method specified" = p_adjust_method %in% c('holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', 'none'))

  doParallel::registerDoParallel(cores = n_cores)
  
  # compute the empirical p.values in parallel
  empirical_p <- foreach::foreach(i = 1:nrow(df), .combine = rbind, .inorder = TRUE) %dopar%
    (perm_cor(x = as.numeric(df[i, ]), y = y, n_perm = n_perm, cor_method = cor_method))
  
  p_cor <- apply(df, 1, function(x) cor(as.numeric(x), y = y, method = cor_method))
  p_res_df <- cbind(df, p_cor, empirical_p)
  p_res_df[p_adjust_method] <- p.adjust(p_res_df$empirical_p, method = p_adjust_method)
  
  # clean up column and rownames
  names(p_res_df)[names(p_res_df) == "p_cor"] <- cor_method
  rownames(p_res_df) <- rownames(df)

  p_res_df
}


#' Sample Random Correlations from a dataframe
#' 
#' Selects n random rows from a dataframe with replacement. For each random row, 
#' permute vector y and perform correlation with the given cor_method.
#' @param df data.frame. Site by Samples data.frame
#' @param y numeric vector. Numeric vector of values used to correlate with each row of df
#' @param n integer. Number of random correlations to return. Default (10000)
#' @param cor_method a character string indicating which correlation coefficient (or covariance) is to be computed.
#'   One of "pearson", "kendall", or "spearman" (default): can be abbreviated.
#' @return named numeric vector of correlation values and row indeces. Vector can contain NAs if row variance == 0.
#' @export
#' @importFrom stats cor
#' @examples
#' \dontrun{
#' library(methylKit)
#'
#' ages <- data.frame(age = c(30, 80, 34, 30, 80, 40))
#'
#' sim_meth <- dataSim(
#'   replicates = 6,
#'   sites = 1000,
#'   treatment = c(rep(1, 3), rep(0, 3)),
#'   covariates = ages,
#'   sample.ids = c(paste0("test", 1:3), paste0("ctrl", 1:3))
#' )
#'
#' # extract the methylation as percentages and coerce to data.frame
#' perc_meth <- as.data.frame(percMethylation(sim_meth))
#'
#' # Get 1_000_000 random correlations from perc_meth data.frame
#' cors <- sample_n_random_cor(perc_meth, y = ages$age)
#' 
#' # histogram of distribution -- filter NAs if present
#' hist(cors[!is.na(cors)])
#' }
sample_n_random_cor = function(df, y, n = 1e+06, cor_method = "spearman") {
  stopifnot("Number of columns is not equal to length of vector y" = length(colnames(df)) == length(y))
  random_rows <- sample(1:nrow(df), n, replace = TRUE)
  cor_results <- apply(df[random_rows, ],
                       MARGIN = 1, 
                       FUN = function(x) {
                         cor(
                           as.numeric(x), 
                           y = sample(y), 
                           method = cor_method)
                         }
                       )
  return(cor_results)
}
