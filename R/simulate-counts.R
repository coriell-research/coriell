#' Simulate RNASeq Count Data
#' 
#' Simulate an RNASeq count matrix. Adapted from edgeR documentation.
#' @param n_genes integer. The number of gene to simulate
#' @param n_up integer. The number of genes the count offset will be added to
#' @param n_down integer. the number of genes the count offset will be subtracted from. Negative numbers will be coerced to zeros.
#' @param n_samples integer. The number of samples to simulate
#' @param groups character vector. Character vector specifying the names of the groups to simulate
#' @param de_group character. The name of the group to be differentially expressed
#' @param mu numeric. Mean value of counts
#' @param phi numeric. Dispersion parameter
#' @param count_offset numeric. Count offset to be applied to the up/down genes.
#' @return list containing the count matrix and a vectors indicating which rows were modified.
#' @examples
#' # simulate counts using default parameters
#' sim <- simulate_counts()
#' 
#' # extract the count matrix from the simulation result
#' counts <- sim$table
#' 
#' # show rows where counts were modified
#' counts[sim$de_genes, ]
#' @export
simulate_counts = function(n_genes = 1000, 
                           n_up = 50,
                           n_down = 50,
                           n_samples = 6, 
                           groups = c("ctl", "trt"),
                           de_group = "trt", 
                           mu = 10, 
                           phi = 0.1, 
                           count_offset = 10) {
  
  stopifnot("de_group must be a member of groups" = de_group %in% groups)
  stopifnot("More groups than samples" = length(groups) <= n_samples)
  stopifnot("n_samples not a multiple of levels(groups)" = n_samples %% length(groups) == 0)
  stopifnot("Sum(n_up, n_down) must be than n_genes" = n_up + n_down <= n_genes)
  
  de_genes <- sample(n_genes, size = n_up + n_down)
  up <- sample(de_genes, size = n_up)
  down <- setdiff(de_genes, up)
  
  # generate the count table
  counts <- matrix(stats::rnbinom(n_genes * n_samples,
                           mu = mu,
                           size = 1 / phi),
                  nrow = n_genes,
                  ncol = n_samples)
  
  colnames(counts) <- rep(groups, each = n_samples / length(groups))
  
  group <- gl(length(groups), n_samples / length(groups), labels = groups)
  
  counts[up, group == de_group] <- counts[up, group == de_group] + count_offset
  counts[down, group == de_group] <- counts[down, group == de_group] - count_offset
  
  # silently convert negative numbers to zero if present
  counts <- apply(counts, 2, function(x) ifelse(x < 0, 0, x))
  
  rownames(counts) <- paste("gene", 1:n_genes, sep = ".")
  
  list("table" = counts, 
       "up_genes" = paste("gene", up, sep = "."),
       "down_genes" = paste("gene", down, sep = "."))
}
