#' Summarize RNA-seq expression results
#'
#' Summarize a results data.frame. Return data.frame of counts of up/down/non-DE 
#' genes based on log-fold-change and significance values.
#' 
#' @param df dataframe of results. Must have columns containing significance values and log-fold changes.
#' @param fdr_col dataframe column (unquoted). Column of dataframe containing the significance level values.
#' @param lfc_col dataframe column (unquoted). Column of dataframe containing the lof-fold change values.
#' @param fdr numeric. FDR or significance value below which genes are considered significant.
#' @param lfc numeric. abs(log-fold change) value above which genes are considered significant.
#' @export
#' @importFrom magrittr %>%
#' @return tibble of summarized results
#' @examples
#' library(coriell)
#' 
#' 
#' # Create a simulated results data.frame
#' # 50 will have logFC > 0
#' # 50 will have logFC < 0
#' # 25 of the up-regulated will have significant FDR values (< 0.0001)
#' # 25 of the down-regulated will have significant FDR values (< 0.0001)
#' df <- data.frame(feature_id = paste("gene", 1:100,sep = "."),
#'                  logFC = c(runif(50, min = -8, max = -0.01), 
#'                            runif(50, min = 0.01, max = 8)),
#'                  FDR = c(runif(25, min = 1e-12, max = 1e-4),
#'                          runif(25, min = 0.06, max = 1),
#'                          runif(25, min = 1e-12, max = 1e-4),
#'                          runif(25, min = 0.06, max = 1))
#' )
#' 
#' # view the results data.frame
#' head(df)
#' 
#' # summarize results
#' summarize_dge(df)
#' 
summarize_dge <- function(df, fdr_col = FDR, lfc_col = logFC, fdr = 0.05, lfc = 0) {
  df %>%
    dplyr::mutate(dge = dplyr::case_when(
      ({{ fdr_col }} < fdr) & (abs({{ lfc_col }}) > lfc) & ({{ lfc_col }} < 0) ~ "down",
      ({{ fdr_col }} < fdr) & (abs({{ lfc_col }}) > lfc) & ({{ lfc_col }} > 0) ~ "up",
      TRUE ~ "non-dge"
    )) %>%
    dplyr::mutate(dge = factor(dge, levels = c("up", "down", "non-dge"))) %>%
    dplyr::group_by(dge, .drop = FALSE) %>%
    dplyr::summarize(n = dplyr::n()) %>%
    dplyr::mutate(perc = n / sum(n) * 100)
}
