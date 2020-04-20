#' Volcano Plot Function
#'
#' Create a volcano plot from the given dataframe
#' @param df dataframe containing columns with gene names, p-values, and log-fold changes.
#' @param x logFC column
#' @param y FDR or p-value column
#' @param lab column containing gene id or labels
#' @param fdr significance level cutoff for plotting. Values below the given fdr threshold are considered significant. Default (0.05)
#' @param lfc log-fold-change cutoff for plotting. Values greater than the abs(lfc) and less than fdr are displayed as differentially expressed. Default (0)
#' @param label_sig TRUE/FALSE. apply ggrepel::geom_text_labels to significant DE genes.
#' @return ggplot volcano plot
#' @export
#' @importFrom rlang .data
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
#' # Create volcano plot
#' plot_volcano(res_df, x = logFC, y = FDR)
#' }
#'

plot_volcano = function(df, x, y, lab = NULL, fdr = 0.05, lfc = 0, label_sig = FALSE) {
  plot_df <- df %>%
    dplyr::mutate(signif = dplyr::if_else({{ y }} < fdr & abs({{ x }}) > lfc, 'yes', 'no'))

  vplot <- ggplot2::ggplot(data = plot_df, ggplot2::aes(x = {{ x }}, y = -log10({{ y }}))) +
           ggplot2::geom_point(ggplot2::aes(color = .data$signif)) +
           ggplot2::scale_colour_manual(values = c('gray40','red2')) +
           ggplot2::geom_vline(xintercept = 0, linetype = 2) +
           ggplot2::geom_vline(xintercept = lfc, linetype = 3) +
           ggplot2::geom_vline(xintercept = -lfc, linetype = 3) +
           ggplot2::geom_hline(yintercept = -log10(fdr), linetype = 3) +
           ggplot2::labs(color = 'Significant') +
           ggplot2::labs(subtitle = paste("FDR = ", fdr, "; lfc cutoff = ", lfc),
                         x = "logFC",
                         y = "-log10(FDR") +
           ggplot2::theme_classic()

  if (label_sig == TRUE) {
    vplot + ggrepel::geom_text_repel(data = plot_df %>% dplyr::filter(signif == 'yes'),
                                     ggplot2::aes(label = {{ lab }}))
  } else {
    vplot
  }
}

#' MD Plot Function
#'
#' Create an MD plot from the given dataframe
#' @param df dataframe containing log-fold-change, p-value, and logCPM columns.
#' @param x column in dataframe containing the logCPM data
#' @param y column in dataframe containing the log-fold-change values
#' @param sig_col column in dataframe containing the results from significance testing.
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
#' # Create md plot
#' plot_md(res_df, x = logCPM, y = logFC, sig_col = FDR)
#' }

plot_md = function(df, x, y, sig_col, fdr = 0.05, lfc = 0) {
  plot_df <- df %>%
    dplyr::mutate(DE = dplyr::case_when({{ sig_col }} < fdr & {{ y }} < -lfc ~ "Down",
                                        {{ sig_col }} < fdr & {{ y }} > lfc ~ "Up",
                                        TRUE ~ "Non-DE"),
                  DE = factor(.data$DE, levels = c("Up", "Non-DE", "Down")))

  md_plot <- ggplot2::ggplot(data = plot_df, ggplot2::aes(x = {{ x }}, y = {{ y }})) +
    ggplot2::geom_point(data = dplyr::filter(plot_df, DE == "Up"), ggplot2::aes(color = .data$DE)) +
    ggplot2::geom_point(data = dplyr::filter(plot_df, DE == "Down"), ggplot2::aes(color = .data$DE)) +
    ggplot2::geom_point(data = dplyr::filter(plot_df, DE == "Non-DE"), ggplot2::aes(color = .data$DE), size = 0.5) +
    ggplot2::scale_color_manual(values = c('Up' = 'red', 'Non-DE' = 'black', 'Down' = 'blue')) +
    ggplot2::geom_hline(yintercept = 0, linetype = 1) +
    ggplot2::geom_hline(yintercept = lfc, linetype = 2) +
    ggplot2::geom_hline(yintercept = -lfc, linetype = 2) +
    ggplot2::labs(subtitle = paste("FDR = ", fdr, "; lfc cutoff = ", lfc),
                  x = "Average logCPM",
                  y = "Log-fold change") +
    ggplot2::theme_classic()

  md_plot
}
