#' Convert EdgeR results object to a data.frame
#'
#' Create a data.frame from an \code{edgeR} results object. This function calls
#' \code{edgeR::topTags()} on the object and extracts the \code{table} data.frame
#' with all features. This function returns all rows unsorted by default
#' i.e. \code{topTags(..., n=Inf, sort.by="none")}.
#'
#' @param x \code{edgeR} results object to be converted
#' @param ... Additional arguments passed to \code{edgeR::topTags()}. \code{n} and
#' \code{sort.by} are fixed and cannot be supplied. If filtering arguments (e.g.
#' \code{p.value}) leave no rows, a 0-row data.frame is returned.
#' @export
#' @return data.frame
#' @examples
#' library(edgeR)
#' library(coriell)
#'
#' # create some fake data
#' x <- data.frame(
#'   ctl1 = rnbinom(1000, size = 0.4, prob = 1e-5),
#'   ctl2 = rnbinom(1000, size = 0.4, prob = 1e-5),
#'   trt1 = rnbinom(1000, size = 0.4, prob = 1e-5),
#'   trt2 = rnbinom(1000, size = 0.4, prob = 1e-5),
#'   row.names = paste0("gene", 1:1000)
#' )
#'
#' # run edger pipeline
#' group <- factor(c(1, 1, 2, 2))
#' y <- DGEList(counts = x, group = group)
#' y <- calcNormFactors(y)
#' design <- model.matrix(~group)
#' y <- estimateDisp(y, design)
#'
#' # To perform quasi-likelihood F-tests:
#' fit <- glmQLFit(y, design)
#' qlf <- glmQLFTest(fit, coef = 2)
#'
#' # convert the results object to a dataframe -- do not filter the results
#' res_df <- edger_to_df(qlf)
#'
#' head(res_df)
#'
edger_to_df <- function(x, ...) {
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("edgeR package is required.")
  }
  tt <- edgeR::topTags(x, n = Inf, sort.by = "none", ...)
  if (inherits(tt, "TopTags")) {
    df <- tt$table
  } else {
    # Filtering left no rows: return a 0-row data.frame with the usual columns
    args <- list(...)
    args[c("p.value", "lfc")] <- NULL
    tt <- do.call(edgeR::topTags, c(list(x, n = Inf, sort.by = "none"), args))
    df <- tt$table[0, , drop = FALSE]
  }

  # Prevent duplicated colnames when adding feature ids
  if ("feature_id" %in% colnames(df)) {
    df <- cbind(feature_id.1 = rownames(df), df)
  } else {
    df <- cbind(feature_id = rownames(df), df)
  }

  rownames(df) <- NULL

  return(df)
}
