#' Correlation Permutation Test Function
#'
#' Function for performing a correlation permutation test
#' @param x a numeric vector
#' @param y a numeric vector
#' @param n_perm numeric. the number of permutations to perform
#' @param method a character string indicating which correlation coefficient (or covariance) is to be computed.
#'   One of "pearson", "kendall", or "spearman" (default): can be abbreviated.
#' @importFrom stats cor
#' @return numeric vector of length n_perm
perm_cor = function(x, y, n_perm, method = "spearman") {
  # calculate the actual value
  test_stat <- cor(x, y, method = method)

  # test_stat will be NA for rows without variance
  if (is.na(test_stat)) {
    return(NA)
  }

  # Perform the permutation test
  cor_results <- vector("double", length = n_perm)
  for (i in 1:n_perm) {
    perm_vec <- sample(x)
    cor_results[[i]] <- cor(perm_vec, y, method = method)
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
#' @param method a character string indicating which correlation coefficient (or covariance) is to be computed.
#'   One of "pearson", "kendall", or "spearman" (default): can be abbreviated.
#' @return df with additional columns for correlations, empirical p-values, and fdr adjusted p-values.
#' @importFrom stats cor p.adjust
#' @importFrom foreach %dopar%
#' @examples
#' \dontrun{
#' library(methylKit)
#'
#' ages = data.frame(age = c(30, 80, 34, 30, 80, 40))
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
#' res <- parallel_permutation_test(perc_meth, y = ages$age, n_cores = 4, n_perm = 1000)
#' }
#' @export
parallel_permutation_test = function(df, y, n_cores = 1, n_perm = 1000, method = "spearman") {
  doParallel::registerDoParallel(cores = n_cores)

  # compute the empirical p.values in parallel
  empirical_p <- foreach::foreach(i = 1:nrow(df), .combine = rbind) %dopar%
    (perm_cor(x = as.numeric(df[i, ]), y = y, n_perm = n_perm, method = method))

  # strip 'result.' from rownames added by foreach call
  rownames(empirical_p) <- gsub("^result.", "", rownames(empirical_p))

  p_cor <- apply(df, 1, function(x) cor(as.numeric(x), y = y, method  = method))
  p_res_df <- cbind(df, "cor" = p_cor, empirical_p)
  p_res_df$fdr <- p.adjust(p_res_df$empirical_p, method = "BH")

  p_res_df
}
