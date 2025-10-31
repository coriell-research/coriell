#' Create a binned histogram of methylation percentage statistics
#'
#' @param x numeric vector of methylation beta-values 
#' @param bar_color color of the histogram bars. Default = 'cornflowerblue'
#' @param border_color color of the borders for each bar. Default = 'white'
#' @param plot_title title of the plot. Default = NULL
#' @param y_label_cex size of the y-axis label. Default = 1
#' @param x_label_cex size of the x-axis label. Default = 1
#' @param label_cex size of the barplot labels. Default = 0.8
#' @param y_axis_label name of the y-axis. Default = 'Number of CpGs'
#' @param y_axis_cex size of the y-axis values. Default = 1
#' @param x_axis_cex size of the x-axis values. Default = 0.8
#'
#' @returns histogram with methylation statistics 
#' @export
#'
#' @examples
#' 
#' # Example methylation data (single column of a matrix)
#' data <- rbeta(runif(1e4), 0.6, 0.8)
#' plot_meth_hist(data, plot_title="Percentage Methylated CpGs per Binned beta-value")
#' 
#' # Plot all columns of a matrix
#' m <- replicate(6, rbeta(1e4, 0.5, 0.8))
#' colnames(m) <- paste("Sample", 1:6)
#' 
#' # Create a square grid
#' n <- ncol(m) 
#' ncols <- ceiling(sqrt(n))
#' nrows <- ceiling(n / ncols)
#' 
#' par(mfrow = c(nrows, ncols))
#' for (i in 1:n) {
#'   plot_meth_hist(as.numeric(m[, i]), plot_title = colnames(m)[i], x_label_cex = 0.8)
#' }
#' 
plot_meth_hist <- function(x, bar_color = "cornflowerblue",
                           border_color = "white", plot_title = NULL,
                           y_label_cex = 1, x_label_cex = 1, label_cex = 0.8,
                           y_axis_cex = 0.8, x_axis_cex = 0.8, 
                           y_axis_label = "Number CpGs") {
  hdata <- hist(x, plot = FALSE, breaks = 10)
  pct <- round(hdata$counts / sum(hdata$counts) * 100, 0)
  col_labels <- paste0(pct, "%")

  bin_starts <- hdata$breaks[-length(hdata$breaks)]
  bin_ends <- hdata$breaks[-1]
  bin_labels <- paste(round(bin_starts, 1), round(bin_ends, 1), sep = " - ")

  bp <- barplot(
    hdata$counts,
    names.arg = bin_labels,
    space = 0,
    col = bar_color,
    border = border_color,
    main = plot_title,
    ylab = y_axis_label,
    ylim = c(0, max(hdata$counts) * 1.2),
    cex.lab = y_label_cex, 
    cex.axis = y_axis_cex,
    cex.names = x_axis_cex,
    las = 2
  )
  text(
    x = bp,
    y = hdata$counts,
    labels = col_labels,
    pos = 3,
    cex = label_cex,
    col = "black"
  )
  axis(side = 1, at = bp, labels = FALSE, lty = 1)
  mtext("Methylation Value", side = 1, line = 4, cex = x_label_cex)
}
