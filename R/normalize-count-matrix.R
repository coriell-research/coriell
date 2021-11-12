#' Compute various normalization methods on a count matrix
#' 
#' Given a matrix of gene counts and a grouping factor identifying the group
#' membership for each column calculate normalized expression values using a 
#' variety of gene expression normalization techniques.
#' 
#' @details The following normalization methods are computed
#' \itemize{
#'  \item{"TMM"}{ Trimmed Mean of M-values computed by \code{edgeR::calcNormFactors}}
#'  \item{"RLE"}{ Relative log expression computed by \code{edgeR::calcNormFactors}}
#'  \item{"UQ"}{ Upper Quartile normalization computed by \code{edgeR::calcNormFactors}}
#'  \item{"VST"}{ Variance Stabalizing Transformation computed by \code{DESeq2::vst}}
#'  \item{"RLog"}{ Regularized Log transformation computed by \code{DESeq2::rlog}}
#'  \item{"QSmooth"}{ Smoothed Quantile Normalization computed by \code{qsmooth::qsmooth}}
#'  \item{"RUVg"}{ Remove Unwanted Variation using control genes by \code{RUVseq::RUVg}}
#'  \item{"Total Count"}{ Library Size Normalization only}
#' }
#' 
#' The function returns a list with the log count per million normalized values 
#' for each of the above methods. The log2 scaled original counts are also returned.
#' For all methods except \code{vst} and \code{rlog} a pseudocount of 2 is added
#' to avoid taking logs of 0.
#' 
#' @param x \code{SummarizedExperiment} object containing a counts assay with 
#' raw counts or a matrix with features as rows and samples as columns.
#' @param groups a group level continuous or categorial covariate associated 
#' with each sample or column in the object. The order of \code{groups} 
#' must match the order of the columns in object.
#' @param filter bool. Default = TRUE. Should the lowly expressed genes be removed
#' prior to normalization? Uses the \code{filterByExpr} function in the \code{edgeR}
#' package.
#' @param min_count Passed to \code{edgeR::filterByExpr} the minimum count required
#'  for at least some samples. Default 10.
#' @param min_total_count Passed to \code{edgeR::filterByExpr} the minimum total count
#' required. Default 15.
#' @param large_n Passed to \code{edgeR::filterByExpr} the number of samples per group
#' that is considered to be "large". Default 10.
#' @param min_prop Passed to \code{edgeR::filterByExpr} the minimum proportion of samples
#' in the smallest group that express the gene. Default 0.7.
#' @param ref_column Column to use as reference for method = "TMM". Can be a column 
#' number or a numeric vector of length nrow(object). Default = NULL. Used by \code{edgeR::calcNormFactors}
#' @param log_ratio_trim the fraction (0 to 0.5) of observations to be trimmed from 
#' each tail of the distribution of log-ratios (M-values) before computing the mean. 
#' Used by method="TMM" for each pair of samples. Default = 0.3. Used by \code{edgeR::calcNormFactors}
#' @param sum_trim the fraction (0 to 0.5) of observations to be trimmed from 
#' each tail of the distribution of A-values before computing the mean. 
#' Used by method="TMM" for each pair of samples. Default = 0.05. Used by \code{edgeR::calcNormFactors}
#' @param do_weighting logical, whether to use (asymptotic binomial precision) 
#' weights when computing the mean M-values. Used by method="TMM" for each pair 
#' of samples. Default = TRUE. Used by \code{edgeR::calcNormFactors}
#' @param acutoff minimum cutoff applied to A-values. Count pairs with lower 
#' A-values are ignored. Used by method="TMM" for each pair of samples. Default = -1e10.
#' Used by \code{edgeR::calcNormFactors}
#' @param p numeric value between 0 and 1 specifying which quantile of the counts 
#' should be used by method="upperquartile". Default = 0.75. Used by \code{edgeR::calcNormFactors}
#' @param n_sub the number of genes to subset. Default = 1000. Used by \code{DESeq2::vst}
#' @param fit_type either "parametric", "local", "mean", or "glmGamPoi" for the 
#' type of fitting of dispersions to the mean intensity. Used by \code{DESeq2::vst}. 
#' Default = "parametric"
#' @param batch. (Optional) batch covariate (multiple batches are not allowed). 
#' If batch covariate is provided, Combat() from sva is used prior to qsmooth 
#' normalization to remove batch effects. See Combat() for more details. Used by
#' \code{qsmooth::qsmooth}
#' @param norm_factors optional normalization scaling factors. Used by \code{qsmooth::qsmooth}
#' @param window window size for running median which is a fraction of the 
#' number of rows in object. Default is 0.05. Used by \code{qsmooth::qsmooth}
#' @param control_genes A character, logical, or numeric vector indicating the 
#' subset of genes to be used as negative controls in the estimation of the 
#' factors of unwanted variation. Used by \code{RUVseq::RUVg}. Default = NULL.
#' @param k The number of factors of unwanted variation to be estimated from the data.
#' Used by \code{RUVseq::RUVg}. Default = 2.
#' @param drop The number of singular values to drop in the estimation of the 
#' factors of unwanted variation. This number is usually zero, but might be set 
#' to one if the first singular value captures the effect of interest. It must 
#' be less than k. Used by \code{RUVseq::RUVg}. Default = 0.
#' @param center 	If TRUE, the counts are centered, for each gene, to have mean 
#' zero across samples. This is important to ensure that the first singular 
#' value does not capture the average gene expression. Used by \code{RUVseq::RUVg}.
#' Default = TRUE.
#' @param epsilon A small constant (usually no larger than one) to be added to 
#' the counts prior to the log transformation to avoid problems with log(0). Used by 
#' \code{RUVseq::RUVg}. Default = 1.
#' @param tolerance Tolerance in the selection of the number of positive singular 
#' values, i.e., a singular value must be larger than tolerance to be considered 
#' positive. Used by \code{RUVseq::RUVg}. Default = 1e-08.
#' @param round If TRUE, the normalized measures are rounded to form pseudo-counts.
#' Used by \code{RUVseq::RUVg}. Default = TRUE.
#' @export
#' @return List with items containing matrices of normalized log count per million
#'  values for each normalization method along with the log-scaled original counts.
#' @examples
#' library(coriell)
#' 
#' 
#' # Simulate count matrix with differential expression
#' sim <- simulate_counts(n_genes = 10000, n_up = 5000, n_down = 500, count_offset = 25)
#' 
#' # Extract 1000 non-DE control genes for RUVg method
#' controls <- setdiff(rownames(sim$table), union(sim$up_genes, sim$down_genes))
#' controls <- sample(controls, 1000, replace = FALSE)
#' 
#' # Define grouping factor
#' Group <- factor(rep(c("Control", "Treatment"), each = 3))
#' 
#' # Run all normalization methods
#' normed <- normalize_count_matrix(x = sim$table, groups = Group, control_genes = controls)
#' 
#' # View the normalized distributions for each
#' for (method in names(normed)) boxplot(normed[[method]], main = method)
normalize_count_matrix <- function(
  x, 
  groups,
  filter = TRUE, 
  min_count = 10,
  min_total_count = 15, 
  large_n = 10,
  min_prop = 0.7, 
  ref_column = NULL,
  log_ratio_trim = 0.3, 
  sum_trim = 0.05,
  do_weighting = TRUE, 
  acutoff = -1e10,
  p = 0.75,
  n_sub = 1000,
  fit_type = c("parametric", "local", "mean", "glmGamPoi"),
  batch = NULL,
  norm_factors = NULL,
  window = 0.05,
  control_genes = NULL,
  k = 1,
  drop = 0,
  center = TRUE,
  round = TRUE,
  epsilon = TRUE,
  tolerance = 1e-08
  ) {
  if (is(x, "SummarizedExperiment")) {
    stopifnot("Could not find 'counts' in assays of x" = "counts" %in% SummarizedExperiment::assayNames(x))
    M <- SummarizedExperiment::assay(se, "counts")
  } else if (is(x, "matrix")) {
    M <- x
  }
  stopifnot("Number of items in groups does not match number of columns in count matrix" = ncol(M) == length(groups))
  
  if (filter) {
    keep <- edgeR::filterByExpr(y = M, group = groups, min.count = min_count, 
                                min.total.count = min_total_count,
                                large.n = large_n, min.prop = min_prop)
    counts <- M[keep, ]
  } else {
    counts <- M
  }
  
  # edgeR normalizations
  dge <- edgeR::DGEList(counts = counts, group = groups)
  tmm <- edgeR::calcNormFactors(dge, method = "TMM", refColumn = ref_column, logratioTrim = log_ratio_trim, sumTrim = sum_trim, doWeighting = do_weighting, Acutoff = acutoff)
  rle <- edgeR::calcNormFactors(dge, method = "RLE")
  uq <- edgeR::calcNormFactors(dge, method = "upperquartile", p = p)
  lcpm.tmm <- edgeR::cpm(tmm, log = TRUE)
  lcpm.rle <- edgeR::cpm(rle, log = TRUE)
  lcpm.uq <- edgeR::cpm(uq, log = TRUE)
  
  # DESeq2 normalizations
  coldata <- data.frame(Group = groups)
  rownames(coldata) <- colnames(counts)
  dds <- DESeq2::DESeqDataSetFromMatrix(counts, colData = coldata, design = ~0 + Group)
  vsd <- DESeq2::vst(dds, blind = FALSE, nsub = n_sub, fitType = fit_type)
  rl <- DESeq2:::rlog(dds, blind = FALSE, fitType = fit_type)
  lcpm.vst <- SummarizedExperiment::assay(vsd)
  lcpm.rlog <- SummarizedExperiment::assay(rl)
  
  # QSmooth normalization
  qs <- qsmooth::qsmooth(object = counts, group_factor = groups, batch = batch, norm_factors = norm_factors, window = window)
  lcpm.qs <- edgeR::cpm(qsmooth::qsmoothData(qs), log = TRUE)
  
  # Total counts per million
  tc <- apply(counts, MARGIN = 2, FUN = function(x) x / sum(x) * 1e06)
  lcpm.tc <- log2(tc + 2)
  
  # RUVseq normalization
  lcpm.ruv <- NULL
  if (!is.null(control_genes)) {
    filtered_controls <- intersect(rownames(counts), control_genes)
    if (length(filtered_controls) != length(control_genes)) message("Some control genes have been dropped in the filtering step. Using the set of control genes present after filtering")
    ruvg <- RUVSeq::RUVg(x = counts, cIdx = filtered_controls, k = k, drop = drop, center = center, round = round, epsilon = epsilon, tolerance = tolerance)
    lcpm.ruv <- edgeR::cpm(ruvg$normalizedCounts, log = TRUE)
  }
  
  return(list(logcounts.tmm = lcpm.tmm,
              logcounts.rle = lcpm.rle,
              logcounts.uq = lcpm.uq,
              logcounts.vst = lcpm.vst,
              logcounts.rlog = lcpm.rlog,
              logcounts.qs = lcpm.qs,
              logcounts.ruv = lcpm.ruv,
              logcounts.tc = lcpm.tc,
              logcounts.orig = log2(counts + 2)))
}
