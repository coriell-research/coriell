#' Combine p-values across multiple differential expression studies
#' 
#' This function provides an interface to p-value combination methods 
#' implemented by the \code{metapod} package for use in the context of 
#' differential expression. 
#' 
#' @param l List of data.frames containing differential expression results. 
#' Every data.frame in the list must have the same columns names.
#' @param p_method p-value combination method. One of c("simes", "holm-min", 
#' "berger", "fisher", "pearson", "wilkinson", "stouffer"). Default "fisher".
#' @param lfc_threshold numeric scalar defining the threshold at which an effect 
#' is "up" or "down"
#' @param pval_col character string of the column name in each data.frame that 
#' contains the p-values to combine. Default "p.value"
#' @param lfc_col character string of the column name in each data.frame that 
#' contains the log-fold change values to combine. Default "logFC"
#' @param feature_col character string of the column name in each data.frame 
#' that contains the feature names (i.e. genes) over which to combine. Default 
#' "gene_id"
#' @param ... Additional arguments passed to \code{metapod::combineParallelPValues()}
#' @export
#' @return data.frame with columns containing the combined p-values, the 
#' summarized direction of the effect, and additional summary statistics about 
#' combined logFC values for each gene.
#' @details
#' The "n.logFC" column of the resulting data.frame reports the number of non-NA
#' logFC values. The "n.pos.logFC" and the "n.neg.logFC" columns report the 
#' number of logFC values across experiments that exceed the value given by the
#' "lfc_threshold" argument. All other summary statistics on logFC values are 
#' performed after removing NA values.
#' @examples
#' df1 <- data.frame(
#'   gene_id = c("A", "B", "C", "D", "E"),
#'   logFC = c(-1, 0, 1.5, 3.2, 0.25),
#'   p.value = c(0.015, 0.1, 0.000005, 0.006, 0.9)
#' )
#' 
#' df2 <- data.frame(
#'   gene_id = c("A", "B", "C", "D", "F"),
#'   logFC = c(-1.5, 0.1, 1.6, 2.9, -2.1),
#'   p.value = c(0.0006, 0.3, 0.000058, 0.0005, 0.0087)
#' )
#' 
#' df3 <- data.frame(
#'   gene_id = c("G", "B", "C", "D", "F"),
#'   logFC = c(2.6, 0.2, 1.25, 3.3, -2.7),
#'   p.value = c(0.085, 0.5, 0.005, 0.0006, 0.00007)
#' )
#' 
#' # Note: genes do not have to be present across all studies 
#' L <- list(exp1=df1, exp2=df2, exp3=df3)
#'
#' meta_de(L)
meta_de <- function(l, p_method = "fisher", lfc_threshold = 0, 
                    pval_col = "p.value", lfc_col = "logFC", 
                    feature_col = "gene_id", ...) {
  
  if (!requireNamespace("metapod", quietly = TRUE)) {
    stop("metapod package is required.")
  }
  dt <- data.table::rbindlist(l, idcol = "Experiment")
  p_dt <- data.table::dcast(dt, get(feature_col) ~ Experiment, value.var = pval_col)
  l_dt <- data.table::dcast(dt, get(feature_col) ~ Experiment, value.var = lfc_col)

  par.combined <- metapod::combineParallelPValues(
    p_dt[, .SD, .SDcols = names(l)], 
    method = p_method,
    ...
    )
  
  par.dir <- metapod::summarizeParallelDirection(
    l_dt[, .SD, .SDcols = names(l)], 
    influential = par.combined$influential,
    threshold = lfc_threshold
    )
  
  # Calculate summary statistics for logFCs
  lfc <- as.matrix(l_dt, rownames = "feature_col")
  summaryStats <- function(x) {
    list(
      n = sum(!is.na(x)),
      n.pos = sum( (x[!is.na(x)]) > lfc_threshold ),
      n.neg = sum( (x[!is.na(x)]) < -lfc_threshold ),
      mean = mean(x, na.rm = TRUE),
      median = median(x, na.rm = TRUE),
      min = min(x, na.rm = TRUE),
      max = max(x, na.rm = TRUE)
      )
  }
  ss <- apply(lfc, 1, summaryStats, simplify = TRUE)

  # Collect results into a final data.frame
  data.frame(
    feature = rownames(lfc),
    p.value = par.combined$p.value,
    direction = par.dir,
    n.logFC = sapply(ss, function(x) x$n, simplify = TRUE),
    n.pos.logFC = sapply(ss, function(x) x$n.pos, simplify = TRUE),
    n.neg.logFC = sapply(ss, function(x) x$n.neg, simplify = TRUE),
    mean.logFC = sapply(ss, function(x) x$mean, simplify = TRUE),
    median.logFC = sapply(ss, function(x) x$median, simplify = TRUE),
    min.logFC = sapply(ss, function(x) x$min, simplify = TRUE),
    max.logFC = sapply(ss, function(x) x$max, simplify = TRUE)
  )
}
