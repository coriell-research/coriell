#' Create a volcano plot from expression data
#'
#' Create a volcano plot from a data.frame containing differential expression results.
#' @param df dataframe containing columns with gene names, p-values, and log-fold changes.
#' @param x logFC column. Default (logFC)
#' @param y FDR or p-value column. Default (FDR)
#' @param lab column containing gene id or labels. Default (NULL)
#' @param fdr significance level cutoff for plotting. Values below the given fdr threshold are considered significant. Default (0.05)
#' @param lfc log-fold-change cutoff for plotting. Values greater than the abs(lfc) and less than fdr are displayed as differentially expressed. Default (0)
#' @param label_sig TRUE/FALSE. apply ggrepel::geom_text_labels to significant DE genes.
#' @param annotate_counts TRUE/FALSE. Annotate the plot with the summarized gene counts
#' @param up_color Point color of the up-regulated features
#' @param down_color Point color of the down-regulated features
#' @param nonde_color Point color of the unperturbed features
#' @param up_alpha Point alpha of the up-regulated features
#' @param down_alpha Point alpha of the down-regulated features
#' @param nonde_alpha Point alpha of the unperturbed features
#' @param up_size Point size of the up-regulated features
#' @param down_size Point size of the down-regulated features
#' @param nonde_size Point size of the unperturbed features
#' @param up_shape Point shape of the up-regulated features
#' @param down_shape Point shape of the down-regulated features
#' @param nonde_shape Point shape of the unperturbed features
#' @param xmin_label_offset numeric. Value between 0 and 1 inclusive to control the x-position of the count annotation label for the 'down' counts
#' @param xmax_label_offset numeric. Value between 0 and 1 inclusive to control the x-position of the count annotation label for the 'up' counts
#' @param ymax_label_offset numeric. Value between 0 and 1 inclusive to control the y-position of the count labels.
#' @param lab_size numeric. Size of the label if annotate_counts = TRUE. Default 6.
#' @param lab_digits numeric. The number of digits used when rounding percentage values when annotate_counts=TRUE. Default (1)
#' @param x_axis_limits numeric vector of axis limits supplied to ggplot2::coord_cartesian(). Default (NULL)
#' @param y_axis_limits numeric vector of axis limits supplied to ggplot2::coord_cartesian(). Default (NULL)
#' @param ... Additional arguments passed to \code{ggrepel::geom_text_repel}
#' @return ggplot volcano plot
#' @import data.table
#' @export
#' @examples
#' plot_volcano(GSE161650_de, fdr = 0.01, lfc = log2(2))
plot_volcano <- function(df, x = "logFC", y = "FDR", lab = NULL, fdr = 0.05,
                         lfc = 0, label_sig = FALSE, annotate_counts = TRUE,
                         up_color = "red2", down_color = "royalblue2",
                         nonde_color = "grey40", up_alpha = 1, down_alpha = 1,
                         nonde_alpha = 1, up_size = 1, down_size = 1,
                         nonde_size = 1, up_shape = 16, down_shape = 16,
                         nonde_shape = ".", xmin_label_offset = 0.8,
                         xmax_label_offset = 0.8, ymax_label_offset = 0.9,
                         lab_size = 6, lab_digits = 1, x_axis_limits = NULL,
                         y_axis_limits = NULL, ...) {
  if (label_sig && is.null(lab)) {
    message("'label_sig = TRUE' but 'lab = NULL'. Please specifiy a column name of features in order to plot labels.")
  }

  dt <- as.data.table(df)

  # Modify the input df for plotting
  dt[, `:=`(
    direction = fcase(
      get(y) < ..fdr & abs(get(x)) > ..lfc & get(x) > 0, "Up",
      get(y) < ..fdr & abs(get(x)) > ..lfc & get(x) < 0, "Down",
      default = "Unperturbed"
    ),
    logPval = -log10(get(y))
  )]

  # Set up the base plot object
  p <- ggplot2::ggplot(data = dt, ggplot2::aes(x = .data[[x]], y = .data[["logPval"]])) +
    ggplot2::geom_point(data = dt[direction == "Unperturbed"], color = nonde_color, alpha = nonde_alpha, size = nonde_size, shape = nonde_shape) +
    ggplot2::geom_point(data = dt[direction == "Down"], color = down_color, alpha = down_alpha, size = down_size, shape = down_shape) +
    ggplot2::geom_point(data = dt[direction == "Up"], color = up_color, alpha = up_alpha, size = up_size, shape = up_shape) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::geom_vline(xintercept = lfc, linetype = 3) +
    ggplot2::geom_vline(xintercept = -lfc, linetype = 3) +
    ggplot2::geom_hline(yintercept = -log10(fdr), linetype = 3) +
    ggplot2::labs(
      caption = paste(y, " = ", fdr, "\nlfc cutoff = ", round(lfc, digits = 2)),
      x = x,
      y = paste0("-log10(", y, ")")
    ) +
    coriell::theme_coriell()

  # Apply new limits if set
  if (!is.null(x_axis_limits) & !(is.null(y_axis_limits))) {
    p <- suppressMessages({
      p + ggplot2::coord_cartesian(
        xlim = x_axis_limits,
        ylim = y_axis_limits
        )
      })
  } else if (!is.null(x_axis_limits)) {
    p <- suppressMessages(p + ggplot2::coord_cartesian(xlim = x_axis_limits))
  } else if (!is.null(y_axis_limits)) {
    p <- suppressMessages(p + ggplot2::coord_cartesian(ylim = y_axis_limits))
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
    d <- coriell::summarize_dge(df, fdr_col = y, lfc_col = x, fdr = fdr, lfc = lfc)
    plot_lims <- coriell::get_axis_limits(p)
    up_count <- d[d$Direction == "Up", "N", drop = TRUE]
    down_count <- d[d$Direction == "Down", "N", drop = TRUE]
    up_pct <- round(d[d$Direction == "Up", "Percent", drop = TRUE], digits = lab_digits)
    down_pct <- round(d[d$Direction == "Down", "Percent", drop = TRUE], digits = lab_digits)

    p <- p +
      ggplot2::annotate(
        geom = "label",
        size = lab_size,
        x = xmin_label_offset * plot_lims$x_min,
        y = ymax_label_offset * plot_lims$y_max,
        label = paste0(down_count, "\n", down_pct, "%")
      ) +
      ggplot2::annotate(
        geom = "label",
        size = lab_size,
        x = xmax_label_offset * plot_lims$x_max,
        y = ymax_label_offset * plot_lims$y_max,
        label = paste0(up_count, "\n", up_pct, "%")
      )
  }
  p
}
