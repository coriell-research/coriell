#' Volcano Plot Function
#'
#' Create a volcano plot from the given dataframe
#' @param df dataframe containing columns with gene names, p-values, and log-fold changes.
#' @param x logFC column. Default (logFC)
#' @param y FDR or p-value column. Default (FDR)
#' @param lab column containing gene id or labels. Default (NULL)
#' @param fdr significance level cutoff for plotting. Values below the given fdr threshold are considered significant. Default (0.05)
#' @param lfc log-fold-change cutoff for plotting. Values greater than the abs(lfc) and less than fdr are displayed as differentially expressed. Default (0)
#' @param label_sig TRUE/FALSE. apply ggrepel::geom_text_labels to significant DE genes.
#' @param annotate_counts TRUE/FALSE. Annotate the plot with the summarized gene counts
#' @param xmin_label_offset numeric. Value between 0 and 1 inclusive to control the x-position of the count annotation label for the 'down' counts
#' @param xmax_label_offset numeric. Value between 0 and 1 inclusive to control the x-position of the count annotation label for the 'up' counts
#' @param ymax_label_offset numeric. Value between 0 and 1 inclusive to control the y-position of the count labels.
#' @return ggplot volcano plot
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
#' # Create volcano plot
#' plot_volcano(res_df)
#' }
#'
plot_volcano <- function(df,
                         x = logFC,
                         y = FDR,
                         lab = feature_id,
                         fdr = 0.05,
                         lfc = 0,
                         label_sig = FALSE,
                         annotate_counts = TRUE,
                         xmin_label_offset = 0.5,
                         xmax_label_offset = 0.5,
                         ymax_label_offset = 0.8) {
  stopifnot("xmin_label_offset must be between 0 and 1" = xmin_label_offset >= 0 & xmin_label_offset <= 1)
  stopifnot("xmax_label_offset must be between 0 and 1" = xmax_label_offset >= 0 & xmax_label_offset <= 1)
  stopifnot("ymax_label_offset must be between 0 and 1" = ymax_label_offset >= 0 & ymax_label_offset <= 1)

  plot_df <- df %>%
    dplyr::mutate(signif = dplyr::if_else({{ y }} < fdr & abs({{ x }}) > lfc, "yes", "no"))

  vplot <- ggplot2::ggplot(data = plot_df, ggplot2::aes(x = {{ x }}, y = -log10({{ y }}))) +
    ggplot2::geom_point(ggplot2::aes(color = .data$signif)) +
    ggplot2::scale_colour_manual(values = c("no" = "gray40", "yes" = "red2")) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::geom_vline(xintercept = lfc, linetype = 3) +
    ggplot2::geom_vline(xintercept = -lfc, linetype = 3) +
    ggplot2::geom_hline(yintercept = -log10(fdr), linetype = 3) +
    ggplot2::labs(color = "Significant") +
    ggplot2::labs(
      subtitle = paste("FDR = ", fdr, "; lfc cutoff = ", round(lfc, digits = 2)),
      x = "logFC",
      y = "-log10(FDR)"
    ) +
    ggplot2::theme_classic()

  # add text labels to significant genes
  if (label_sig) {
    vplot <- vplot +
      ggrepel::geom_text_repel(
        data = plot_df %>% dplyr::filter(signif == "yes"),
        ggplot2::aes(label = {{ lab }})
      )
  }

  if (annotate_counts) {
    d <- coriell::summarize_dge(df, fdr = fdr, lfc = lfc)
    plot_lims <- coriell::get_axis_limits(vplot)

    up_count <- d[d$dge == "up", "n", drop = TRUE]
    down_count <- d[d$dge == "down", "n", drop = TRUE]
    up_pct <- round(d[d$dge == "up", "perc", drop = TRUE], digits = 2)
    down_pct <- round(d[d$dge == "down", "perc", drop = TRUE], digits = 2)

    vplot <- vplot +
      ggplot2::annotate(
        geom = "label",
        x = xmin_label_offset * plot_lims$x_min,
        y = ymax_label_offset * plot_lims$y_max,
        label = paste0(down_count, "\n", down_pct, "%")
      ) +
      ggplot2::annotate(
        geom = "label",
        x = xmax_label_offset * plot_lims$x_max,
        y = ymax_label_offset * plot_lims$y_max,
        label = paste0(up_count, "\n", up_pct, "%")
      )
  }
  vplot
}
