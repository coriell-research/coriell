#' Normalize counts in a DGEList object
#'
#' This function wraps additional normalization methods (qsmooth and cyclic loess) for a given
#' DGEList object in addition to the standard normalization methods implemented by
#' \code{edgeR::normLibSizes()}.
#'
#' @details
#' If "TMM", "TMMwsp", "RLE", "upperquartile", or "none" is selected, then
#' \code{edgeR::normLibSizes(x, ...)} is used. If method = "qsmooth" then
#' \code{qsmooth::qsmooth(x, group_factor, ...)} is used. If method = "loess" then
#' \code{csaw::normOffsets(x, ...)} is used.
#'
#' @param x DGEList object
#' @param method Normalization method. One of "TMM", "TMMwsp", "RLE", "upperquartile", "none",
#' "qsmooth", or "loess".
#' @param refColumn column to use as reference for method="TMM". Can be a column number or a numeric
#' vector of length nrow(object).
#' @param logratioTrim the fraction (0 to 0.5) of observations to be trimmed from each tail of the
#' distribution of log-ratios (M-values) before computing the mean. Used by method="TMM" for each
#' pair of samples.
#' @param sumTrim fraction (0 to 0.5) of observations to be trimmed from each tail of the
#' distribution of A-values before computing the mean. Used by method="TMM" for each pair of samples.
#' @param doWeighting logical, whether to use (asymptotic binomial precision) weights when computing
#' the mean M-values. Used by method="TMM" for each pair of samples.
#' @param Acutoff minimum cutoff applied to A-values. Count pairs with lower A-values are ignored.
#' Used by method="TMM" for each pair of samples.
#' @param p numeric value between 0 and 1 specifying which quantile of the counts should be used
#' by method="upperquartile".
#' @param group_factor a group level continuous or categorial covariate associated with each sample
#' or column in the object. The order of the group_factor must match the order of the columns in
#' object. Used in qsmooth normalization.
#' @param batch (Optional) batch covariate (multiple batches are not allowed). If batch covariate
#' is provided, \code{ComBat()} from sva is used prior to qsmooth normalization to remove batch effects.
#' See \code{ComBat()} for more details. Used in qsmooth.
#' @param norm_factors optional normalization scaling factors. Used in qsmooth.
#' @param window window size for running median which is a fraction of the number of rows in
#' object. Default is 0.05. Used in qsmooth.
#' @param prior_count prior count to add to 0 counts when computing offsets. Default = 0.1. Used
#' in qsmooth.
#' @param weights numeric vector of non-negative prior weights. Missing values are treated as zero.
#' Default = NULL. Used in loess.
#' @param span positive numeric value between 0 and 1 specifying proportion of data to be used in
#' the local regression moving window. Larger numbers give smoother fits. Default = 0.3. Used in
#' loess.
#' @param iterations number of local regression fits. Values greater than 1 produce robust fits.
#' Default = 4L. Used in loess.
#' @param min.weight minimum weight. Any lower weights will be reset. Default = 1e-5. Used in loess.
#' @param max.weight maximum weight. Any higher weights will be reset. Default = 1e5. Used in loess.
#' @param equal.weights.as.null should equal weights be treated as if weights were NULL, so that
#' lowess is called? Applies even if all weights are all zero. Default = TRUE. Used in loess.
#' @param loess_method method used for weighted lowess. Possibilities are "weightedLowess",
#' "loess" or "locfit". Used in loess.
#'
#' @returns DGEList with updated norm.factors (for standard methods) or an observation-level offset matrix in x$offset (for qsmooth and loess).
#' @export
#' @examples
#' counts <- matrix(rnbinom(1e3, mu = 10, size = 20), ncol = 20)
#'
#' y <- edgeR::DGEList(
#'   counts = counts,
#'   group = gl(n = 2, k = 10, labels = c("control", "treatment"))
#' )
#'
#' # TMM normalization
#' y.tmm <- normalize_counts(y, method = "TMM")
#'
#' # Qsmooth normalization
#' y.qs <- normalize_counts(y, method = "qsmooth", group_factor = y$samples$group)
#'
#' # Cyclic loess
#' y.loess <- normalize_counts(y, method = "loess")
#'
normalize_counts <- function(
  x,
  method = c(
    "TMM",
    "TMMwsp",
    "RLE",
    "upperquartile",
    "none",
    "qsmooth",
    "loess"
  ),
  refColumn = NULL,
  logratioTrim = 0.3,
  sumTrim = 0.05,
  doWeighting = TRUE,
  Acutoff = -1e10,
  p = 0.75,
  group_factor = NULL,
  batch = NULL,
  norm_factors = NULL,
  window = 0.05,
  prior_count = 0.1,
  weights = NULL,
  span = 0.3,
  iterations = 4L,
  min.weight = 1e-5,
  max.weight = 1e5,
  equal.weights.as.null = TRUE,
  loess_method = c("weightedLowess", "loess", "locfit")
) {
  method <- match.arg(method)
  loess_method <- match.arg(loess_method)

  if (!inherits(x, "DGEList")) {
    stop(paste("Expected DGEList as input. x is", class(x)[1]))
  }

  if (method %in% c("TMM", "TMMwsp", "RLE", "upperquartile", "none")) {
    if (!requireNamespace("edgeR", quietly = TRUE)) {
      stop(paste0("edgeR package is required for method=", method))
    }
    result <- edgeR::normLibSizes(
      x,
      method = method,
      refColumn = refColumn,
      logratioTrim = logratioTrim,
      sumTrim = sumTrim,
      doWeighting = doWeighting,
      Acutoff = Acutoff,
      p = p
    )
  }

  if (method == "qsmooth") {
    if (!requireNamespace("qsmooth", quietly = TRUE)) {
      stop("qsmooth package is required for method='qsmooth'")
    }
    if (is.null(group_factor)) {
      stop("group_factor must be supplied for method='qsmooth'")
    }
    qs <- qsmooth::qsmooth(
      x$counts,
      group_factor = group_factor,
      batch = batch,
      norm_factors = norm_factors,
      window = window
    )
    qsd <- qsmooth::qsmoothData(qs)
    offset <- log(x$counts + prior_count) - log(qsd + prior_count)
    result <- edgeR::scaleOffset(x, offset)
  }

  if (method == "loess") {
    if (!requireNamespace("csaw", quietly = TRUE)) {
      stop("csaw package is required for method='loess'")
    }
    result <- csaw::normOffsets(
      x,
      se.out = TRUE,
      weights = weights,
      span = span,
      iterations = iterations,
      min.weight = min.weight,
      max.weight = max.weight,
      equal.weights.as.null = equal.weights.as.null,
      method = loess_method
    )
  }

  return(result)
}
