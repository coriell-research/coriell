#' edgeR result to dataframe
#'
#' Create a dataframe from an edger results object
#'
#' @param res_obj edgeR results object to be converted
#' @param fdr numeric. adjusted p-value under which genes are considered significant
#' @param lfc numeric. abs(log-fold-change) value above which genes are considered significant
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
edger_to_df = function(res_obj, fdr = 1, lfc = 0) {
  edgeR::topTags(res_obj, n = nrow(res_obj$table))$table %>%
    dplyr::as_tibble(rownames = 'gene_id') %>%
    dplyr::mutate(signif = dplyr::if_else(.data$FDR < fdr & abs(.data$logFC) > lfc, 'yes', 'no'),
                  direction = dplyr::case_when(.data$FDR < fdr & lfc < 0 ~ "down",
                                               .data$FDR < fdr & lfc > 0 ~ "up",
                                               TRUE ~ "non-de"))
}
