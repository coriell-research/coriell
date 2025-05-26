#' Create a methylation bead plot
#'
#' @param x read x CpG matrix of 0's and 1's indicating unmethylated and methylated CpGs, respectively
#' @param m_col bead color of methylated CpGs. default "black"
#' @param un_col bead color of unmethylated CpGs. default "white"
#' @param pt_size size of bead set by cex. default 1 
#' @param x_lab x-axis label. default "CpGs"
#' @param y_lab y-axis label. default "Reads"
#' @param lab_cex x/y-axis label size. default 1.3
#' @param lab_font x/y-axis label font. default 2
#' @param lab_pos x/y-axis label position relative to the plot window. default 1.2
#' @param dash_col color of the horizontal dashed lines. default "grey80"
#'
#' @returns plot
#' @export
#'
#' @examples
#' n_reads <- 10
#' n_cpgs <- 10
#' m <- matrix(rbinom(n_reads*n_cpgs, 1, 0.25), ncol=n_cpgs, nrow=n_reads)
#' 
#' plot_beads(m, pt_size=2)
plot_beads <- function(x, m_col = "black", un_col = "white", pt_size = 1,
                       x_lab = "CpGs", y_lab = "Reads", lab_cex = 1.3,
                       lab_font = 2, lab_pos = 1.2, dash_col = "grey80") {
  plot.new()
  plot.window(xlim = c(1, ncol(x)), ylim = c(nrow(x), 1))
  abline(h = 1:nrow(x), col = dash_col, lty = "dotted", lwd = 1)
  xs <- rep(1:ncol(x), each = nrow(x))
  ys <- rep(1:nrow(x), times = ncol(x))
  cols <- as.vector(ifelse(x == 1, m_col, un_col))
  points(x = xs, y = ys, col = "black", pch = 21, bg = cols, cex = pt_size)
  title(xlab = x_lab, col.lab = "black", adj = 0, cex.lab = lab_cex, font.lab = lab_font, line = lab_pos)
  title(ylab = y_lab, col.lab = "black", adj = 1, cex.lab = lab_cex, font.lab = lab_font, line = lab_pos)
}