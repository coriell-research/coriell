#' Create an MD plot from expression data
#'
#' Create an MD (MA) plot from a data.frame containing differential expression results.
#'
#' @param df dataframe containing log-fold-change, p-value, and logCPM columns.
#' @param x column in dataframe containing the logCPM data. Default ("logCPM")
#' @param y column in dataframe containing the log-fold-change values. Default ("logFC")
#' @param sig_col column in dataframe containing the results from significance testing. Default ("FDR")
#' @param lab column in dataframe containing the labels to plot if label_sig = TRUE. Default NULL
#' @param fdr numeric. Significance level cutoff for plotting. Values below the given fdr threshold are considered significant. Default (0.05)
#' @param lfc numeric. Log-fold-change cutoff for plotting. Values greater than the abs(lfc) and less than fdr are displayed as differentially expressed. Default(0)
#' @param annotate_counts TRUE/FALSE. Annotate the plot with the summarized gene counts
#' @param label_sig logical. Apply \code{ggrepel::geom_text_repel()} to significant DE genes. Default FALSE
#' @param up_color Point color of the up-regulated features. Default ("red2")
#' @param down_color Point color of the down-regulated features. Default ("royalblue2")
#' @param nonde_color Point color of the unperturbed features. Default ("grey40")
#' @param up_alpha Point alpha value of the up-regulated features. Default (1)
#' @param down_alpha Point alpha value of the down-regulated features. Default (1)
#' @param nonde_alpha Point alpha value of the unperturbed features. Default (1)
#' @param up_size Point size of the up-regulated features. Default (1)
#' @param down_size Point size of the down-regulated features. Default (1)
#' @param nonde_size Point size of the unperturbed features. Default (1)
#' @param up_shape Point shape of the up-regulated features
#' @param down_shape Point shape of the down-regulated features
#' @param nonde_shape Point shape of the unperturbed features
#' @param xmax_label_offset numeric. Value between 0 and 1 inclusive. Controls the x-position of the count labels
#' @param ymax_label_offset numeric. Value between 0 and 1 inclusive. Controls the y-position of the 'up' count label
#' @param ymin_label_offset numeric. Value between 0 and 1 inclusive. Controls the y-position of the 'down' count label
#' @param lab_size numeric. If annotate_counts = TRUE specify the label size. Default = 8.
#' @param lab_digits numeric. The number of digits used when rounding percentage values when annotate_counts=TRUE. Default (2)
#' @param x_axis_limits numeric vector of axis limits supplied to ggplot2::xlim(). Default (NULL)
#' @param y_axis_limits numeric vector of axis limits supplied to ggplot2::ylim(). Default (NULL)
#' @param ... Additional arguments passed to \code{ggrepel::geom_text_repel()}
#' @return ggplot MD plot
#' @import data.table
#' @export
#' @examples
#' plot_md(GSE161650_de, fdr = 0.01, lfc = log2(2))
plot_md <- function(df, x = "logCPM", y = "logFC", sig_col = "FDR", lab = NULL,
                    fdr = 0.1, lfc = 0, annotate_counts = TRUE,
                    label_sig = FALSE, up_color = "red2",
                    down_color = "royalblue2", nonde_color = "grey40",
                    up_alpha = 1, down_alpha = 1, nonde_alpha = 1,
                    up_size = 1, down_size = 1, nonde_size = 1,
                    up_shape = 19, down_shape = 19, nonde_shape = 19,
                    xmax_label_offset = 0.8, ymax_label_offset = 0.5,
                    ymin_label_offset = 0.5, lab_size = 8, lab_digits = 2,
                    x_axis_limits = NULL, y_axis_limits = NULL, ...) {
  if (label_sig && is.null(lab)) {
    message("'label_sig = TRUE' but 'lab = NULL'. Please specifiy a column name of features in order to plot labels.")
  }

  # Add new label for Up, Down and Non-DE genes
  dt <- as.data.table(df)
  dt[, direction := fcase(
    get(sig_col) < ..fdr & abs(get(y)) > ..lfc & get(y) > 0, "Up",
    get(sig_col) < ..fdr & abs(get(y)) > ..lfc & get(y) < 0, "Down",
    default = "Unperturbed"
  )]

  # Set up the base plot object
  p <- ggplot2::ggplot(data = dt, ggplot2::aes(x = .data[[x]], y = .data[[y]])) +
    ggplot2::geom_point(data = dt[direction == "Unperturbed"], color = nonde_color, size = nonde_size, alpha = nonde_alpha, shape = nonde_shape) +
    ggplot2::geom_point(data = dt[direction == "Down"], color = down_color, size = down_size, alpha = down_alpha, shape = down_shape) +
    ggplot2::geom_point(data = dt[direction == "Up"], color = up_color, size = up_size, alpha = up_alpha, shape = up_shape) +
    ggplot2::geom_hline(yintercept = 0, linetype = 1) +
    ggplot2::geom_hline(yintercept = lfc, linetype = 2) +
    ggplot2::geom_hline(yintercept = -lfc, linetype = 2) +
    ggplot2::labs(
      caption = paste0("FDR = ", fdr, "\nlfc cutoff = ", round(lfc, digits = 2)),
      x = "Average logCPM",
      y = "Log-fold change"
    ) +
    coriell::theme_coriell()

  # Apply new limits if set
  if (!is.null(x_axis_limits)) {
    p <- p + ggplot2::xlim(x_axis_limits)
  }
  if (!is.null(y_axis_limits)) {
    p <- p + ggplot2::ylim(y_axis_limits)
  }

  # Add text labels to significant genes
  if (label_sig && !is.null(lab)) {
    if (!requireNamespace("ggrepel", quietly = TRUE)) {
      stop("ggrepel package is required.")
    }

    p <- p +
      ggrepel::geom_text_repel(
        data = dt[direction %in% c("Up", "Down")],
        ggplot2::aes(label = .data[[lab]]),
        ...
      )
  }

  # Add DE count annotations
  if (annotate_counts) {
    d <- coriell::summarize_dge(df, fdr_col = sig_col, lfc_col = y, fdr = fdr, lfc = lfc)
    plot_lims <- coriell::get_axis_limits(p)
    up_count <- d[d$Direction == "Up", "N", drop = TRUE]
    down_count <- d[d$Direction == "Down", "N", drop = TRUE]
    up_pct <- round(d[d$Direction == "Up", "Percent", drop = TRUE], digits = lab_digits)
    down_pct <- round(d[d$Direction == "Down", "Percent", drop = TRUE], digits = lab_digits)

    p <- p +
      ggplot2::annotate(
        geom = "label",
        size = lab_size,
        x = xmax_label_offset * plot_lims$x_max,
        y = ymax_label_offset * plot_lims$y_max,
        label = paste0(up_count, "\n", up_pct, "%")
      ) +
      ggplot2::annotate(
        geom = "label",
        size = lab_size,
        x = xmax_label_offset * plot_lims$x_max,
        y = ymin_label_offset * plot_lims$y_min,
        label = paste0(down_count, "\n", down_pct, "%")
      )
  }
  p
}
