#' MD Plot Function
#'
#' Create an MD plot from the given dataframe
#' @param df dataframe containing log-fold-change, p-value, and logCPM columns.
#' @param x column in dataframe containing the logCPM data. Default (logCPM)
#' @param y column in dataframe containing the log-fold-change values. Default (logFC)
#' @param sig_col column in dataframe containing the results from significance testing. Default (FDR)
#' @param fdr numeric. Significance level cutoff for plotting. Values below the given fdr threshold are considered significant. Default (0.05)
#' @param lfc numeric. Log-fold-change cutoff for plotting. Values greater than the abs(lfc) and less than fdr are displayed as differentially expressed. Default(0)
#' @return ggplot MD plot
#' @export
#' @importFrom rlang .data
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
#' # convert the results object to a dataframe
#' res_df <- edger_to_df(qlf)
#'
#' # Create md plot
#' plot_md(res_df)
#' }
#'
plot_md <- function(df, x = logCPM, y = logFC, sig_col = FDR, fdr = 0.05, lfc = 0) {
  plot_df <- df %>%
    dplyr::mutate(
      DE = dplyr::case_when(
        {{ sig_col }} < fdr & {{ y }} < -lfc ~ "Down",
        {{ sig_col }} < fdr & {{ y }} > lfc ~ "Up",
        TRUE ~ "Non-DE"
      ),
      DE = factor(.data$DE, levels = c("Up", "Non-DE", "Down"))
    )
  
  md_plot <- ggplot2::ggplot(data = plot_df, ggplot2::aes(x = {{ x }}, y = {{ y }})) +
    ggplot2::geom_point(data = dplyr::filter(plot_df, DE == "Non-DE"), ggplot2::aes(color = .data$DE)) +
    ggplot2::geom_point(data = dplyr::filter(plot_df, DE == "Up"), ggplot2::aes(color = .data$DE)) +
    ggplot2::geom_point(data = dplyr::filter(plot_df, DE == "Down"), ggplot2::aes(color = .data$DE)) +
    ggplot2::scale_color_manual(values = c("Up" = "red", "Non-DE" = "black", "Down" = "blue")) +
    ggplot2::geom_hline(yintercept = 0, linetype = 1) +
    ggplot2::geom_hline(yintercept = lfc, linetype = 2) +
    ggplot2::geom_hline(yintercept = -lfc, linetype = 2) +
    ggplot2::labs(
      subtitle = paste0("FDR = ", fdr, "; lfc cutoff = ", lfc),
      x = "Average logCPM",
      y = "Log-fold change"
    ) +
    ggplot2::theme_classic()
  
  md_plot
}