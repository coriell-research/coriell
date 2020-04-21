#' edgeR result to dataframe
#'
#' Create a dataframe from an edger results object
#'
#' @param res_obj edgeR results object to be converted
#' @export
#' @return tibble
#' @examples
#' \dontrun{
#' library(edgeR)
#' library(coriell)
#'
#' # create some fake data
#' x <- data.frame(ctl1 = rnbinom(1000, size = 0.4, prob = 1e-5),
#'                 ctl2 = rnbinom(1000, size = 0.4, prob = 1e-5),
#'                 trt1 = rnbinom(1000, size = 0.4, prob = 1e-5),
#'                 trt2 = rnbinom(1000, size = 0.4, prob = 1e-5),
#'                 row.names = paste0('gene', 1:1000))
#'
#' # run edger pipeline
#' group <- factor(c(1,1,2,2))
#' y <- DGEList(counts=x,group=group)
#' y <- calcNormFactors(y)
#' design <- model.matrix(~group)
#' y <- estimateDisp(y,design)
#'
#' # To perform quasi-likelihood F-tests:
#' fit <- glmQLFit(y,design)
#' qlf <- glmQLFTest(fit,coef=2)
#'
#' # convert the results object to a dataframe
#' res_df <- edger_to_df(qlf)
#'
#' }
#'
edger_to_df = function(res_obj) {
  edgeR::topTags(res_obj, n = nrow(res_obj$table))$table %>%
    dplyr::as_tibble(rownames = 'gene_id')
}

#' Summarize Results
#'
#' Summarize a results dataframe. Return dataframe of counts of up/down/non-DE genes based on
#' log-fold-change and significance values
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
#' x <- data.frame(ctl1 = rnbinom(1000, size = 0.4, prob = 1e-5),
#'                 ctl2 = rnbinom(1000, size = 0.4, prob = 1e-5),
#'                 trt1 = rnbinom(1000, size = 0.4, prob = 1e-5),
#'                 trt2 = rnbinom(1000, size = 0.4, prob = 1e-5),
#'                 row.names = paste0('gene', 1:1000))
#'
#' # run edger pipeline
#' group <- factor(c(1,1,2,2))
#' y <- DGEList(counts=x,group=group)
#' y <- calcNormFactors(y)
#' design <- model.matrix(~group)
#' y <- estimateDisp(y,design)
#'
#' # To perform quasi-likelihood F-tests:
#' fit <- glmQLFit(y,design)
#' qlf <- glmQLFTest(fit,coef=2)
#'
#' # convert the results object to a dataframe and summarize
#' res_df <- edger_to_df(qlf) %>%
#'   summarize_dge(fdr_col = FDR, lfc_col = logFC)
#'
#' }
summarize_dge = function(df, fdr_col, lfc_col, fdr = 0.05, lfc = 1.5) {
  df %>%
    dplyr::mutate(dge = dplyr::case_when(({{ fdr_col }} < fdr) & (abs({{ lfc_col }}) > lfc) & ({{ lfc_col }} < 0) ~ "down",
                                         ({{ fdr_col }} < fdr) & (abs({{ lfc_col }}) > lfc) & ({{ lfc_col }} > 0) ~ "up",
                                         TRUE ~ "non-dge")) %>%
    dplyr::group_by(dge) %>%
    dplyr::summarize(n = n()) %>%
    dplyr::mutate(perc = n / sum(n) * 100)
}
