#' Volcano Plot Function
#'
#' Create a volcano plot from an edgeR results object
#' @param result edgeR results object. Required
#' @param graph_title Title of the plot. Required.
#' @param fdr Significance level cutoff for plotting. Values below the given fdr threshold are considered significant. Default (0.05)
#' @param lfc Log-fold-change cutoff for plotting. Values greater than the abs(lfc) and less than fdr are displayed as differentially expressed. Default (0)
#' @param label_sig Plot the gene id labels on the plot with ggrepel. Default (FALSE)
#' @param pt_alpha Alpha level of the points. Default (0.25)
#' @param pt_size Size of the points. Default (2.25)
#' @return ggplot volcano plot
#' @export
#' @keywords plot, volcano, volcano plot
#' @importFrom rlang .data
#'
plot_volcano = function(result, graph_title, label_sig = FALSE, fdr = 0.05, lfc = 0, pt_alpha = 0.25, pt_size = 2.25) {
  vplot <- edgeR::topTags(result, n = nrow(result$table))$table %>%
    dplyr::as_tibble(rownames = 'gene_id') %>%
    dplyr::mutate(signif = dplyr::if_else(.data$FDR < fdr & abs(.data$logFC) > lfc, 'yes', 'no')) %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$logFC, y = -log10(.data$FDR))) +
    ggplot2::geom_point(ggplot2::aes(color = .data$signif), alpha = pt_alpha, size = pt_size) +
    ggplot2::scale_colour_manual(values = c('gray40','red2')) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::labs(color = 'Significant') +
    ggplot2::labs(title = graph_title,
                  subtitle = paste("FDR = ", fdr, "; lfc cutoff = ", lfc),
                  x = "logFC",
                  y = "-log10(FDR") +
    ggplot2::theme_classic()

  if (label_sig == FALSE) {
    vplot
  } else {
    vplot + ggrepel::geom_text_repel(data = .data %>% dplyr::filter(.data$signif == 'yes'), ggplot2::aes(label = .data$gene_id))
  }

}

#' MD plot function
#'
#' Create an MD plot from an edgeR results object
#' @param result edgeR result object. Required.
#' @param graph_title Title of the plot. Required.
#' @param fdr Significance level cutoff for plotting. Values below the given fdr threshold are considered significant. Default (0.05)
#' @param lfc Log-fold-change cutoff for plotting. Values greater than the abs(lfc) and less than fdr are displayed as differentially expressed. Default(0)
#' @param up_alpha Alpha value of points for up-regulated genes. Default (0.8)
#' @param down_alpha Alpha value of points for down-regulated genes. Default (0.8)
#' @param non_alpha Alpha value of points for non-differentially expressed genes. Default (0.1)
#' @return ggplot MD plot
#' @export
#' @importFrom rlang .data
#'
plot_MD = function(result, graph_title, fdr = 0.05, lfc = 0, up_alpha = 0.8, down_alpha = 0.8, non_alpha = 0.1) {
  edgeR::topTags(result, n = nrow(result$table))$table %>%
  dplyr::as_tibble(rownames = 'gene_id') %>%
  dplyr::mutate(DE = dplyr::case_when(.data$FDR < fdr & .data$logFC < -lfc ~ "Down",
                                      .data$FDR < fdr & .data$logFC > lfc ~ "Up",
                                      TRUE ~ "Non-DE"),
           DE = factor(.data$DE, levels = c("Up", "Non-DE", "Down"))) %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$logCPM, y = .data$logFC)) +
    ggplot2::geom_point(data = .data %>% dplyr::filter(.data$DE == "Up"), ggplot2::aes(color = .data$DE), alpha = up_alpha) +
    ggplot2::geom_point(data = .data %>% dplyr::filter(.data$DE == "Down"), ggplot2::aes(color = .data$DE), alpha = down_alpha) +
    ggplot2::geom_point(data = .data %>% dplyr::filter(.data$DE == "Non-DE"), ggplot2::aes(color = .data$DE), alpha = non_alpha, size = 0.5) +
    ggplot2::scale_color_manual(values = c('Up' = 'red', 'Non-DE' = 'black', 'Down' = 'blue')) +
    ggplot2::geom_hline(yintercept = 0, linetype = 1) +
    ggplot2::geom_hline(yintercept = lfc, linetype = 2) +
    ggplot2::geom_hline(yintercept = -lfc, linetype = 2) +
    ggplot2::labs(title = graph_title,
         subtitle = paste("FDR = ", fdr, "; lfc cutoff = ", lfc),
         x = "Average logCPM",
         y = "Log-fold change") +
    ggplot2::theme_classic()
}
