#' edgeR result to dataframe
#'
#' Create a dataframe from an edger results object. The function can also be used to filter the results object by the FDR and/or log-fold-change values by 
#' supplying fdr and lfc arguments. By default the function returns all rows of the results object.
#'
#' @param res_obj edgeR results object to be converted
#' @param fdr numeric. FDR adjusted p-value used as a filter. Return only records with this value or less.
#' @param lfc numeric. Log-fold-change used as filter. Return only records where abs(logFC) > lfc.
#' @export
#' @return tibble
#' @examples
#' \dontrun{
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
#' }
#'
edger_to_df <- function(res_obj, fdr = 1, lfc = 0) {
  edgeR::topTags(res_obj, n = nrow(res_obj$table))$table %>%
    dplyr::as_tibble(rownames = "feature_id") %>% 
    dplyr::filter(FDR <= fdr & abs(logFC) >= lfc)
}

#' Summarize Results
#'
#' Summarize a results dataframe. Return dataframe of counts of up/down/non-DE genes based on
#' log-fold-change and significance values. 
#' @param df dataframe of results. Must have columns containing significance values and log-fold changes.
#' @param fdr_col dataframe column. Column of dataframe containing the significance level values.
#' @param lfc_col dataframe column. Column of dataframe containing the lof-fold change values.
#' @param fdr numeric. FDR or significance value below which genes are considered significant.
#' @param lfc numeric. abs(log-fold change) value above which genes are considered significant.
#' @export
#' @return tibble
#' @examples
#' \dontrun{
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
#' # convert the results object to a dataframe and summarize
#' edger_to_df(qlf) %>%
#'   summarize_dge()
#' }
summarize_dge <- function(df, fdr_col = FDR, lfc_col = logFC, fdr = 0.05, lfc = 0) {
  df %>%
    dplyr::mutate(dge = dplyr::case_when(
      ({{ fdr_col }} < fdr) & (abs({{ lfc_col }}) > lfc) & ({{ lfc_col }} < 0) ~ "down",
      ({{ fdr_col }} < fdr) & (abs({{ lfc_col }}) > lfc) & ({{ lfc_col }} > 0) ~ "up",
      TRUE ~ "non-dge"
    )) %>%
    dplyr::group_by(dge) %>%
    dplyr::summarize(n = n()) %>%
    dplyr::mutate(perc = n / sum(n) * 100)
}
