#' Create an MD plot from expression data
#'
#' Create an MD plot from the given dataframe
#' @param df dataframe containing log-fold-change, p-value, and logCPM columns.
#' @param x column in dataframe containing the logCPM data. Default (logCPM)
#' @param y column in dataframe containing the log-fold-change values. Default (logFC)
#' @param sig_col column in dataframe containing the results from significance testing. Default (FDR)
#' @param fdr numeric. Significance level cutoff for plotting. Values below the given fdr threshold are considered significant. Default (0.05)
#' @param lfc numeric. Log-fold-change cutoff for plotting. Values greater than the abs(lfc) and less than fdr are displayed as differentially expressed. Default(0)
#' @param annotate_counts TRUE/FALSE. Annotate the plot with the summarized gene counts
#' @param xmax_label_offset numeric. Value between 0 and 1 inclusive. Controls the x-position of the count labels
#' @param ymax_label_offset numeric. Value between 0 and 1 inclusive. Controls the y-position of the 'up' count label
#' @param ymin_label_offset numeric. Value between 0 and 1 inclusive. Controls the y-position of the 'down' count label
#' @return ggplot MD plot
#' @export
#' @importFrom magrittr %>%
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
plot_md <- function(df,
                    x = logCPM,
                    y = logFC,
                    sig_col = FDR,
                    fdr = 0.05,
                    lfc = 0,
                    annotate_counts = TRUE,
                    xmax_label_offset = 0.8,
                    ymax_label_offset = 0.5,
                    ymin_label_offset = 0.5) {
  stopifnot("xmax_label_offset must be between 0 and 1" = xmax_label_offset >= 0 & xmax_label_offset <= 1)
  stopifnot("ymax_label_offset must be between 0 and 1" = ymax_label_offset >= 0 & ymax_label_offset <= 1)
  stopifnot("ymin_label_offset must be between 0 and 1" = ymin_label_offset >= 0 & ymin_label_offset <= 1)

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
      subtitle = paste0("FDR = ", fdr, "; lfc cutoff = ", round(lfc, digits = 2)),
      x = "Average logCPM",
      y = "Log-fold change"
    ) +
    ggplot2::theme_classic()


  if (annotate_counts) {
    d <- coriell::summarize_dge(df, fdr_col = {{ sig_col }}, lfc_col = {{ y }}, fdr = fdr, lfc = lfc)
    plot_lims <- coriell::get_axis_limits(md_plot)

    up_count <- d[d$dge == "up", "n", drop = TRUE]
    down_count <- d[d$dge == "down", "n", drop = TRUE]
    up_pct <- round(d[d$dge == "up", "perc", drop = TRUE], digits = 2)
    down_pct <- round(d[d$dge == "down", "perc", drop = TRUE], digits = 2)

    md_plot <- md_plot +
      ggplot2::annotate(
        geom = "label",
        x = xmax_label_offset * plot_lims$x_max,
        y = ymax_label_offset * plot_lims$y_max,
        label = paste0(up_count, "\n", up_pct, "%")
      ) +
      ggplot2::annotate(
        geom = "label",
        x = xmax_label_offset * plot_lims$x_max,
        y = ymin_label_offset * plot_lims$y_min,
        label = paste0(down_count, "\n", down_pct, "%")
      )
  }
  md_plot
}
