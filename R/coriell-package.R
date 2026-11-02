#' @keywords internal
"_PACKAGE"

#' @importFrom grDevices adjustcolor hcl.colors
#' @importFrom graphics abline axis barplot bxp grid hist legend lines mtext
#'   pairs par plot.new plot.window points rect smoothScatter strheight
#'   strwidth text title
#' @importFrom methods is
#' @importFrom stats IQR coef complete.cases cor cor.test dist fisher.test
#'   kruskal.test lm lowess mad median p.adjust pnorm prcomp quantile var
#'   wilcox.test
#' @importFrom utils head modifyList
#' @importFrom ggplot2 .data
NULL

# Column names used via data.table non-standard evaluation
utils::globalVariables(c(
  ".", "..fdr", "..lfc", ".sample", "Chrom", "Location", "Me", "N",
  "P.Value", "Sample", "Start", "Un", "Var", "direction"
))
