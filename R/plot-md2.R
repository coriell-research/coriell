#' Create an MD plot from expression data using base R
#'
#' This function differs from \code{plot_md()} as it produces a plot in base R.
#' Points are colored by adjusted p-value (and logFC) significance, and the up/down
#' counts are written beneath the plot rather than on the plot area.
#'
#' @param df data.frame containing columns for average abundance, log-fold-change and adjusted p-values
#' @param x column in dataframe containing the logCPM data. Default "logCPM"
#' @param y column in dataframe containing the log-fold-change values. Default "logFC"
#' @param sig_col column in dataframe containing the adjusted p-values used to define significance. Default "FDR"
#' @param lab column in dataframe containing the labels to plot if label_sig = TRUE. Default NULL
#' @param fdr Significance threshold for features. Based on sig_col. Default 0.1
#' @param lfc logFC threshold for considering significant features. Features must have abs(logFC) greater than lfc. Default 0.
#' @param label_sig TRUE/FALSE. Label significant features with the values in \code{lab} using \code{text()}.
#' @param n_labels integer. Maximum number of features to label when label_sig = TRUE. Since base R cannot repel labels, only the \code{n_labels} significant features with the most extreme combined rank of adjusted p-value and absolute logFC are labeled. Default 10
#' @param annotate_counts TRUE/FALSE. Write the summarized feature counts beneath the plot. Default TRUE
#' @param up_color Point color of the up-regulated features
#' @param down_color Point color of the down-regulated features
#' @param nonde_color Point color of the unperturbed features
#' @param up_alpha Point alpha of the up-regulated features
#' @param down_alpha Point alpha of the down-regulated features
#' @param nonde_alpha Point alpha of the unperturbed features
#' @param up_size Point size (cex) of the up-regulated features
#' @param down_size Point size (cex) of the down-regulated features
#' @param nonde_size Point size (cex) of the unperturbed features
#' @param up_shape Point shape (pch) of the up-regulated features
#' @param down_shape Point shape (pch) of the down-regulated features
#' @param nonde_shape Point shape (pch) of the unperturbed features
#' @param x_axis_limits numeric vector of length 2 giving the x-axis limits. Default (NULL)
#' @param y_axis_limits numeric vector of length 2 giving the y-axis limits. Default (NULL)
#' @param ... Additional arguments passed to \code{plot()}. These override the defaults (e.g. \code{xlab}, \code{ylab}, \code{las}, \code{bty}, \code{main})
#' @return Invisibly, a list with the plotted data and group counts
#' @export
#' @examples
#' plot_md2(GSE161650_de, fdr = 0.01, lfc = log2(2))
plot_md2 <- function(
  df,
  x = "logCPM",
  y = "logFC",
  sig_col = "FDR",
  lab = NULL,
  fdr = 0.1,
  lfc = 0,
  label_sig = FALSE,
  n_labels = 10,
  annotate_counts = TRUE,
  up_color = "red2",
  down_color = "royalblue2",
  nonde_color = "grey40",
  up_alpha = 1,
  down_alpha = 1,
  nonde_alpha = 1,
  up_size = 1,
  down_size = 1,
  nonde_size = 1,
  up_shape = 16,
  down_shape = 16,
  nonde_shape = ".",
  x_axis_limits = NULL,
  y_axis_limits = NULL,
  ...
) {
  # Basic validation
  if (!is.data.frame(df)) {
    stop("df must be a data.frame.")
  }

  required_columns <- c(x, y, sig_col, lab)

  missing_columns <- setdiff(required_columns, colnames(df))

  if (length(missing_columns) > 0) {
    stop(
      "The following columns are missing from df: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  if (length(fdr) != 1 || is.na(fdr) || fdr <= 0 || fdr >= 1) {
    stop("fdr must be a single number between 0 and 1.")
  }

  if (length(lfc) != 1 || is.na(lfc) || lfc < 0) {
    stop("lfc must be a single non-negative number.")
  }

  if (length(n_labels) != 1 || is.na(n_labels) || n_labels < 0) {
    stop("n_labels must be a single non-negative number.")
  }

  for (lim in list(x_axis_limits, y_axis_limits)) {
    if (!is.null(lim) && (length(lim) != 2 || anyNA(lim))) {
      stop("Axis limits must be NULL or a numeric vector of length 2.")
    }
  }

  if (label_sig && is.null(lab)) {
    message(
      "'label_sig = TRUE' but 'lab = NULL'. Please specify a column name of features in order to plot labels."
    )
  }

  # Extract relevant columns
  plot_data <- data.frame(
    abundance = as.numeric(df[[x]]),
    effect = as.numeric(df[[y]]),
    p_adjusted = as.numeric(df[[sig_col]])
  )

  if (!is.null(lab)) {
    plot_data$label <- as.character(df[[lab]])
  }

  # Retain complete, valid observations
  keep <- complete.cases(plot_data[, c("abundance", "effect", "p_adjusted")]) &
    is.finite(plot_data$abundance) &
    is.finite(plot_data$effect) &
    plot_data$p_adjusted >= 0 &
    plot_data$p_adjusted <= 1

  plot_data <- plot_data[keep, , drop = FALSE]

  if (nrow(plot_data) == 0) {
    stop("No valid observations remain after filtering.")
  }

  # Significance is based on adjusted p-values and the absolute logFC threshold
  plot_data$significant <- plot_data$p_adjusted < fdr &
    abs(plot_data$effect) > lfc

  plot_data$group <- "Unperturbed"
  plot_data$group[plot_data$significant & plot_data$effect > 0] <- "Up"
  plot_data$group[plot_data$significant & plot_data$effect < 0] <- "Down"

  counts <- c(
    Up = sum(plot_data$group == "Up"),
    Down = sum(plot_data$group == "Down"),
    Unperturbed = sum(plot_data$group == "Unperturbed")
  )

  # Plot limits
  if (is.null(x_axis_limits)) {
    x_axis_limits <- range(plot_data$abundance)
    if (diff(x_axis_limits) == 0) {
      x_axis_limits <- x_axis_limits + c(-1, 1)
    }
  }

  if (is.null(y_axis_limits)) {
    y_axis_limits <- range(plot_data$effect, 0)
    if (diff(y_axis_limits) == 0) {
      y_axis_limits <- y_axis_limits + c(-1, 1)
    }
  }

  # Leave room beneath the plot for the annotations
  old_par <- par(mar = c(if (annotate_counts) 7 else 5, 4.5, 2, 1))
  on.exit(par(old_par), add = TRUE)

  # Empty plot, then draw groups in order so significant points sit on top
  plot_args <- modifyList(
    list(
      x = NA,
      xlim = x_axis_limits,
      ylim = y_axis_limits,
      xlab = "Average logCPM",
      ylab = "Log-fold change",
      las = 1,
      bty = "n"
    ),
    list(...)
  )
  do.call(plot, plot_args)

  draw_group <- function(group, color, alpha, size, shape) {
    d <- plot_data[plot_data$group == group, , drop = FALSE]
    if (nrow(d) > 0) {
      points(
        d$abundance,
        d$effect,
        pch = shape,
        cex = size,
        col = adjustcolor(color, alpha.f = alpha)
      )
    }
  }

  draw_group("Unperturbed", nonde_color, nonde_alpha, nonde_size, nonde_shape)
  draw_group("Down", down_color, down_alpha, down_size, down_shape)
  draw_group("Up", up_color, up_alpha, up_size, up_shape)

  # Reference lines
  abline(h = 0, lty = 1)

  if (lfc > 0) {
    abline(h = c(-lfc, lfc), lty = 2)
  }

  # Label significant features
  if (label_sig && !is.null(lab)) {
    sig <- plot_data[plot_data$significant, , drop = FALSE]
    if (nrow(sig) > 0 && n_labels > 0) {
      # Rank by adjusted p-value and by |logFC| (1 = most extreme); keep the best combined ranks
      score <- rank(sig$p_adjusted, ties.method = "average") +
        rank(-abs(sig$effect), ties.method = "average")
      sig <- sig[order(score)[seq_len(min(n_labels, nrow(sig)))], , drop = FALSE]
      text(sig$abundance, sig$effect, labels = sig$label, pos = 3)
    }
  }

  # Annotations beneath the plot
  if (annotate_counts) {
    mtext(
      sprintf(
        "Down-regulated (%s)",
        format(counts[["Down"]], big.mark = ",")
      ),
      side = 1,
      line = 4.5,
      adj = 0.15,
      col = down_color,
      font = 2
    )
    mtext(
      sprintf("Up-regulated (%s)", format(counts[["Up"]], big.mark = ",")),
      side = 1,
      line = 4.5,
      adj = 0.85,
      col = up_color,
      font = 2
    )
  }

  mtext(
    sprintf("%s < %.3g; lfc cutoff = %.2f", sig_col, fdr, lfc),
    side = 1,
    line = if (annotate_counts) 5.75 else 4,
    cex = 0.8
  )

  invisible(
    list(
      data = plot_data,
      counts = counts
    )
  )
}
