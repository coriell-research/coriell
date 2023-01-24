#' Summarize RNA-seq expression results
#'
#' Summarize a results data.frame. Return data.frame of counts of up/down/non-DE 
#' genes based on log-fold-change and significance values.
#' 
#' @param df dataframe of results. Must have columns containing significance values and log-fold changes.
#' @param fdr_col Column name of the data.frame containing the significance level values.
#' @param lfc_col Column name of data.frame containing the log-fold change values.
#' @param fdr numeric. FDR or significance value below which genes are considered significant.
#' @param lfc numeric. abs(log-fold change) value above which genes are considered significant.
#' @export
#' @return data.frame with columns describing the count and percentages of up/down/unperturbed genes
#' @examples
#' summarize_dge(GSE161650_de)
summarize_dge <- function(df, fdr_col = "FDR", lfc_col = "logFC", fdr = 0.05, lfc = 0) {
  stopifnot("Provided fdr_col column not in data.frame" = fdr_col %in% colnames(df))
  stopifnot("Provided lfc_col not in data.frame" = lfc_col %in% colnames(df))
  
  up <- sum(df[[fdr_col]] < fdr & abs(df[[lfc_col]]) > lfc & df[[lfc_col]] > 0, na.rm = TRUE)
  down <- sum(df[[fdr_col]] < fdr & abs(df[[lfc_col]]) > lfc & df[[lfc_col]] < 0, na.rm = TRUE)
  unperturbed <- nrow(df) - (up + down)
  total <- up + down + unperturbed
    
  data.frame(
    Direction = factor(c("Up", "Down", "Unperturbed"), levels = c("Up", "Down", "Unperturbed")),
    N = c(up, down, unperturbed),
    Percent = c(round(up / total * 100, 2),
                round(down / total * 100, 2),
                round(unperturbed / total * 100, 2)
                )
    )
}
