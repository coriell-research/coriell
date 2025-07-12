#' Plot pairwise correlations between columns of a matrix
#' 
#' Plot the pairwise correlations between pairs of samples in a matrix. 
#' Histograms, correlation values, and smooth scatters with loess fits are 
#' visualized in each panel. Inspired by `methylKit::getCorrelation(..., plot=TRUE)`
#' and examples in `pairs()`
#'
#' @param m Matrix with features in rows and Samples in columns
#' @param hist_breaks Number of breaks to use in histogram
#' @param hist_col Color of the histogram bars
#' @param cor_method Method to use for calculating correlations. One of 
#' "pearson" (default), "kendall", or "spearman": can be abbreviated
#' @param cor_use an optional character string giving a method for computing 
#' covariances in the presence of missing values. This must be 
#' (an abbreviation of) one of the strings "everything", "all.obs", 
#' "complete.obs", "na.or.complete", or "pairwise.complete.obs".
#' @param cor_digits Number of digits to report in correlations 
#' @param cor_cex Character expansion size for correlation values
#' @param scatter_colramp function accepting an integer n as an argument and 
#' returning n colors used in the smoothScatter plots
#' @param scatter_ab_col Color of the abline of the smoothScatter plots
#' @param scatter_ab_lty Line type of the abline of the smoothScatter plots
#' @param scatter_ab_lwd Line weight of the abline of the smoothScatter plots 
#' @param scatter_crv_col Color of the lowess lines of the smoothScatter plots
#' @param scatter_crv_lwd Line weight of the lowess lines of the smoothScatter plots
#' @param scatter_crv_lty Line type of the lowess lines of the smoothScatter plots
#' @param cex_labels Text panel graphics parameters passed to `pairs()`
#' @param font_labels Text panel graphics parameters passed to `pairs()`
#' @param ... Additional arguments passed to `pairs()`
#'
#' @return pairs plot
#' @export
#'
#' @examples
#' 
#' plot_cor_pairs(GSE161650_lc)
#' 
plot_cor_pairs <- function(m, hist_breaks = 30, hist_col = "grey50", 
                           cor_method = c("pearson", "kendall", "spearman"), 
                           cor_use = "everything", cor_digits = 2, 
                           cor_cex = 2, scatter_colramp = hcl.colors, 
                           scatter_ab_col = "red2", scatter_ab_lty = 2,
                           scatter_ab_lwd = 1, scatter_crv_col = "red3", 
                           scatter_crv_lwd = 2, scatter_crv_lty = 1,
                           cex_labels = 2, font_labels = 1, ...) {
  panel_hist <- function(x) {
    usr <- par("usr")
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE, breaks = hist_breaks)
    breaks <- h$breaks
    nB <- length(breaks)
    y <- h$counts
    y <- y / max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col = hist_col)
  }
  
  panel_cor <- function(x, y, ...) {
    par(usr = c(0, 1, 0, 1))
    r <- cor(x, y, method = cor_method, use = cor_use)
    txt <- format(r, nsmall = 2, digits = cor_digits)
    text(0.5, 0.5, txt, cex = cor_cex)
  }
  
  panel_smoothScatter <- function(x, y) {
    smoothScatter(x, y, nrpoints = 0, add = TRUE, colramp = hcl.colors)
    abline(0, 1, lty = scatter_ab_lty, col = scatter_ab_col, lwd = scatter_ab_lwd)
    lines(lowess(x, y), col = scatter_crv_col, lwd = scatter_crv_lwd, lty = scatter_crv_lty)
  }
  
  pairs(
    m, 
    panel = panel_smoothScatter ,
    lower.panel = panel_cor,
    diag.panel = panel_hist,
    horOdd = TRUE,
    cex.labels = cex_labels, 
    font.labels = font_labels,
    ...
  )
}
