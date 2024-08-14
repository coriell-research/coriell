#' Helper function for converting columns to assay data
#' @noRd
#' @keywords internal
.col2assay <- function(df, rows, cols, vals) {
  m <- dcast(df, get(rows) ~ get(cols), value.var = vals, fill = NA)
  as.matrix(m, rownames = "rows")
}

#' Helper function for extracting representative results for summary stats
#' @noRd
#' @keywords internal
.getRepresentative <- function(x, m) {
  by_gene <- asplit(m, 1)
  repVal <- vector("numeric", length(by_gene))
  for (i in seq_along(by_gene)) {
    repVal[i] <- by_gene[[i]][x[i]]
  }
  repVal
}

#' Convert a list of differential expression data.frames to a SummarizedExperiment
#'
#' This function takes as input a list of data.frames containing
#' differential expression results and converts this list to a
#' \code{SummarizedExperiment} object containing assays for each of the
#' reported statistics columns. This function is intended to be used upstream
#' of the meta-analysis functions implemented in this package.
#'
#' @param x List of data.frames containing differential expression results. All
#' data.frames must have matching colnames.
#' @param feature_col Column name in the data.frames containing the gene or
#' feature ids. Default "feature_id"
#' @param import Character vector of columns from the data.frames to import. These
#' columns will be converted to assays in the final SummarizedExperiment object.
#' @param complete Use only features found across all datasets. Default FALSE,
#' i.e. fill data for missing features with NAs.
#' @return SummarizedExperiment object containing assays for each of the columns
#' in 'import'.
#' @export
#' @examples
#'
#' # data.frames containing differential expression data
#' exp1 <- data.frame(
#'   feature_id = c("geneA", "geneB", "geneC"),
#'   PValue = c(0.01, 0.5, 0.05),
#'   FDR = c(0.02, 0.5, 0.07),
#'   logFC = c(1.2, -2.5, 3.7),
#'   logCPM = c(12, 9, 0)
#' )
#'
#' exp2 <- data.frame(
#'   feature_id = c("geneA", "geneB", "geneD"),
#'   PValue = c(0.07, 0.3, 0.8),
#'   FDR = c(0.08, 0.4, 1.0),
#'   logFC = c(1.5, -2.0, 3.0),
#'   logCPM = c(14, 10, 2)
#' )
#'
#' # Combine into a single list
#' l <- list(experiment1 = exp1, experiment2 = exp2)
#'
#' # Convert the data to a SummarizedExperiment
#' se <- dfs2se(l)
#'
#' # Data is converted to assays
#' SummarizedExperiment::assays(se)
dfs2se <- function(x, feature_col = "feature_id",
                   import = c("logFC", "logCPM", "PValue", "FDR"),
                   complete = FALSE) {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("SummarizedExperiment package is required.")
  }
  same_cols <- all(sapply(x, function(l) all(colnames(l) == colnames(x[[1]]))))
  stopifnot("All column names must match across experiments" = same_cols)
  has_cols <- all(import %in% colnames(x[[1]]))
  stopifnot("One of the needed columns is missing from the data.frames" = has_cols)

  universe <- NULL
  if (isTRUE(complete)) {
    universe <- Reduce(base::intersect, lapply(x, function(df) df[[feature_col]]))
  }

  dt <- data.table::rbindlist(x, idcol = "Experiment")
  if (!is.null(universe)) {
    dt <- dt[get(feature_col) %chin% universe]
  }

  assays <- vector("list", length(import))
  names(assays) <- import
  for (i in seq_along(import)) {
    assays[[i]] <- .col2assay(dt, feature_col, "Experiment", import[[i]])
  }

  SummarizedExperiment::SummarizedExperiment(assays = assays)
}

#' Perform p-value combination for sets of differential expression tests
#'
#' This function performs p-value combination for all genes and estimates summary
#' statistics for average effect sizes for all experiments in the input
#' \code{SummarizedExperiment} object.
#'
#' @param x SummarizedExperiment object containing combined differential
#' expression results from different studies. The SE object must contain at
#' least two assays, one for the P-values to combine and the other for the
#' effect sizes to compute (e.g. logFC).
#' @param FUN One of the 'parallel' functions provided by \code{metapod}. One
#' of "parallelBerger", "parallelFisher", "parallelHolmMin", "parallelPearson",
#' "parallelSimes", "parallelStouffer", or "parallelWilkinson".
#' @param pval assay name in SE object containing the P-values to combine.
#' @param lfc assay name in the SE object containing the logFC values to combine.
#' @param impute_missing TRUE/FALSE should missing values in the logFC and P-Value 
#' assays be imputed prior p-value combination? Default TRUE, missing p-values 
#' are imputed with 1 and missing logFCs are imputed with 0.  
#' @param ... Additional arguments passed to FUN. See the \code{metapod} package
#' for details.
#' @return data.table with summary stats of the p-value combination of all
#' experiments. Please see the documentation in the \code{metapod} package for
#' more details. The returned columns, "Rep.LogFC" and "Rep.Pval" contain the
#' results of extracting the representative effect and P=value from all
#' influential tests. These are individual tests in the data that are particularly
#' important for calculating the combined effects.
#' @export
#' @examples
#' # Example taken from ?dfs2se()
#'
#' # Define two differential expression dataset data.frames
#' exp1 <- data.frame(
#'   feature_id = c("geneA", "geneB", "geneC"),
#'   PValue = c(0.01, 0.5, 0.05),
#'   FDR = c(0.02, 0.5, 0.07),
#'   logFC = c(1.2, -2.5, 3.7),
#'   logCPM = c(12, 9, 0)
#' )
#'
#' exp2 <- data.frame(
#'   feature_id = c("geneA", "geneB", "geneD"),
#'   PValue = c(0.07, 0.3, 0.8),
#'   FDR = c(0.08, 0.4, 1.0),
#'   logFC = c(1.5, -2.0, 3.0),
#'   logCPM = c(14, 10, 2)
#' )
#'
#' # Combine into a single list
#' l <- list(experiment1 = exp1, experiment2 = exp2)
#'
#' # Convert the data to a SummarizedExperiment
#' se <- dfs2se(l)
#'
#' # Perform p-value combination across experiments for each gene
#' #  using Wilkinson's method and passing additional values
#' result <- meta_de(se, metapod::parallelWilkinson, min.prop = 0.1)
#' head(result)
#' 
meta_de <- function(x, FUN, pval = "PValue", lfc = "logFC", impute_missing = TRUE, ...) {
  stopifnot("SummarizedExperiment object expected" = is(x, "SummarizedExperiment"))
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("SummarizedExperiment package is required.")
  }
  if (!requireNamespace("metapod", quietly = TRUE)) {
    stop("metapod package is required.")
  }
  if (!requireNamespace("matrixStats", quietly = TRUE)) {
    stop("matrixStats package is required.")
  }
  
  lfc_m <- SummarizedExperiment::assay(x, lfc)
  pval_m <- SummarizedExperiment::assay(x, pval)
  if (isTRUE(impute_missing)) {
    lfc_m[is.na(lfc_m)] <- 0
    pval_m[is.na(pval_m)] <- 1
  }

  # Split the matrices into lists by sample
  pval_l <- asplit(pval_m, 2)
  lfc_l <- asplit(lfc_m, 2)

  # Compute the combined p-values
  comb <- do.call(FUN, list(p.values = pval_l, ...))

  # Summarize the direction
  direction <- metapod::summarizeParallelDirection(lfc_l, influential = comb$influential)

  # Calculate summary statistics of the logFC
  median_lfc <- matrixStats::rowMedians(lfc_m, na.rm = TRUE, useNames = FALSE)
  avg_lfc <- matrixStats::rowMeans2(lfc_m, na.rm = TRUE, useNames = FALSE)
  min_lfc <- matrixStats::rowMins(lfc_m, na.rm = TRUE, useNames = FALSE)
  max_lfc <- matrixStats::rowMaxs(lfc_m, na.rm = TRUE, useNames = FALSE)

  # Combine all results into a data.table
  data.table(
    Feature = rownames(x),
    Combined.Pval = comb$p.value,
    Direction = direction,
    Rep.logFC = .getRepresentative(comb$representative, lfc_m),
    Rep.Pval = .getRepresentative(comb$representative, pval_m),
    Median.logFC = median_lfc,
    Mean.logFC = avg_lfc,
    Min.logFC = min_lfc,
    Max.logFC = max_lfc
  )
}

#' Perform jackknife resampling on all columns of a SummarizeExperiment object
#' 
#' This function provides a simple wrapper to perform jackknife resampling on
#' all columns of a SummarizedExperiment object and returns the results of each 
#' resample in a list. This function is designed to be used to assess the 
#' robustness of p-value combination techniques of the included \code{meta_de()} 
#' function but in theory any arbitrary function which operates on the columns 
#' of a SummarizedExperiment object could be used.
#' 
#' @param x SummarizedExperiment object to perform jackknife resampling of columns on
#' @param FUN Function to perform on each resample.
#' @param ... Additional arguments passed to FUN
#' @return List of jackknife resampled results
#' @export
#' @examples
#' # Define three differential expression dataset data.frames
#' exp1 <- data.frame(
#'   feature_id = c("geneA", "geneB", "geneC"),
#'   PValue = c(0.01, 0.5, 0.05),
#'   FDR = c(0.02, 0.5, 0.07),
#'   logFC = c(1.2, -2.5, 3.7),
#'   logCPM = c(12, 9, 0)
#' )
#'
#' exp2 <- data.frame(
#'   feature_id = c("geneA", "geneB", "geneD"),
#'   PValue = c(0.07, 0.3, 0.8),
#'   FDR = c(0.08, 0.4, 1.0),
#'   logFC = c(1.5, -2.0, 3.0),
#'   logCPM = c(14, 10, 2)
#' )
#' 
#' exp3 <- data.frame(
#'   feature_id = c("geneA", "geneB", "geneC", "geneD"),
#'   PValue = c(0.03, 0.3, 0.01, 0.8),
#'   FDR = c(0.08, 0.4, 0.04, 0.9),
#'   logFC = c(1.5, -2.0, 3.0, 4.1),
#'   logCPM = c(14, 10, 1, 2.1)
#' )
#'
#' # Combine into a single list
#' l <- list(experiment1 = exp1, experiment2 = exp2, experiment3 = exp3)
#'
#' # Convert the data to a SummarizedExperiment
#' se <- dfs2se(l)
#' 
#' # Perform the jackknife using meta_de on each subset of the data
#' result <- jackknifeSE(se, \(x) meta_de(x, metapod::parallelWilkinson, min.prop = 0.5))
#' 
#' # Combine the results from calling meta_de on each resample and show
#' result <- data.table::rbindlist(result, idcol = "Jackknife")
#' head(result[order(Feature)])
jackknifeSE <- function(x, FUN, ...) {
  stopifnot("SummarizedExperiment object expected" = is(x, "SummarizedExperiment"))
  idx <- 1:ncol(x)
  lapply(idx, function(i) do.call(FUN, list(x = x[, setdiff(idx, i)]), ...))
}
