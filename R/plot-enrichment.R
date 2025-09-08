#' GSEA enrichment plot
#'
#' Drop-in replacement for fgsea::plotEnrichment() which allows more plot
#' customization and annotation. The source code for this function was lifted
#' from: \url{https://github.com/alserglab/fgsea/blob/master/R/plot.R}.
#'
#' @param pathway Gene set to plot
#' @param stats Gene-level statistics
#' @param gseaParam  GSEA parameter
#' @param annotate Should the supplied FDR, ES, and NES values be annotated on the plot? default (FALSE)
#' @param FDR numeric. manually supplied FDR value to be annotated on the plot if annotate=TRUE
#' @param ES numeric. manually supplied ES value to be annotated on the plot if annotate=TRUE
#' @param NES numeric. manually supplied NES value to be annotated on the plot if annotate=TRUE
#' @param ticksSize size of the rank tick mark
#' @param tickAlpha alpha value of the rank tick marks
#' @param tickColor color of the rank tick marks
#' @param lineWidth width of the ES line
#' @param lineColor color = the ES line. "red2" or "blue2" if lineColor=NULL,
#' depending on the max observed ES
#' @param lineAlpha alpha value of the ES line
#' @param hlinesColor color of the min/max horizontal ES lines
#' @param hlinesType linetype of the min/max horizontal ES lines
#' @param anno_x_pos if annotate=TRUE, relative x-position of the text annotation on the plot. default (0.8)
#' @param anno_y_pos if annotate=TRUE, relative y-position of the text annotation on the plot. default (0.8)
#' @param anno_size if annotate=TRUE, size of the text annotation
#' @param fdr_digits if annotate=TRUE, number of digits to round supplied FDR value by
#' @param es_digits if annotate=TRUE, number of digits to round supplied ES value by
#' @param nes_digits if annotate=TRUE, number of digits to round supplied NES value by
#'
#' @returns ggplot2 object
#' @export
#'
#' @examples
#'
#' pathway <- fgsea::examplePathways[["5991130_Programmed_Cell_Death"]]
#' stats <- fgsea::exampleRanks
#' plot_enrichment(pathway, stats)
#'
#' # Stat annotations can be added to the plot
#' fgsea_res <- fgsea::fgsea(fgsea::examplePathways, fgsea::exampleRanks)
#' cell_death <- fgsea_res[pathway == "5991130_Programmed_Cell_Death"]
#'
#' plot_enrichment(
#'   pathway,
#'   stats,
#'   annotate = TRUE,
#'   FDR = cell_death[, padj],
#'   ES = cell_death[, ES],
#'   NES = cell_death[, NES],
#'   plotTitle = "Programmed Cell Death"
#' )
#'
plot_enrichment <- function(pathway, stats, gseaParam = 1, ticksSize = 0.5,
                            annotate = FALSE, FDR = NULL, ES = NULL, NES = NULL,
                            tickAlpha = 1, tickColor = "black",
                            lineWidth = 1, lineColor = NULL, lineAlpha = 1,
                            hlinesColor = "black", hlinesType = 2,
                            anno_x_pos = 0.8, anno_y_pos = 0.7, anno_size = 10,
                            fdr_digits = 3, es_digits = 2, nes_digits = 2,
                            plotTitle = NULL, xlab = "Rank", ylab = "Enrichment Score") {
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("fgsea package is required.")
  }

  pd <- fgsea::plotEnrichmentData(
    pathway = pathway,
    stats = stats,
    gseaParam = gseaParam
  )

  # Determine the line color as either red ('up') or blue ('down')
  if (is.null(lineColor)) {
    lineColor <- "red2"
    if (abs(pd$negES) > pd$posES) {
      lineColor <- "blue2"
    }
  }

  # Determine the base plot
  p <- ggplot2::ggplot(pd$curve, ggplot2::aes(x = rank, y = ES)) +
    ggplot2::geom_line(linewidth = lineWidth, color = lineColor, alpha = lineAlpha) +
    ggplot2::geom_segment(
      data = pd$ticks,
      mapping = ggplot2::aes(
        x = rank,
        y = -pd$spreadES * 1 / 20,
        xend = rank,
        yend = pd$spreadES * 1 / 20
      ),
      linewidth = ticksSize,
      color = tickColor,
      alpha = tickAlpha,
    ) +
    ggplot2::geom_hline(yintercept = pd$posES, colour = hlinesColor, linetype = hlinesType) +
    ggplot2::geom_hline(yintercept = pd$negES, colour = hlinesColor, linetype = hlinesType) +
    ggplot2::geom_hline(yintercept = 0, colour = "black") +
    coriell::theme_coriell() +
    ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "grey90")
    ) +
    ggplot2::labs(
      title = plotTitle,
      x = xlab,
      y = ylab
    )

  # Add annotations if requested
  if (isTRUE(annotate)) {
    l <- coriell::get_axis_limits(p)

    if (is.null(FDR)) {
      stop("annotate=TRUE but no FDR is given")
    }
    if (is.null(ES)) {
      stop("annotate=TRUE but no ES is given")
    }
    if (is.null(NES)) {
      stop("annotate=TRUE but no NES is given")
    }

    lab <- paste(
      "FDR:", round(FDR, fdr_digits),
      "\nES:", round(ES, es_digits),
      "\nNES:", round(NES, nes_digits)
    )
    p + ggplot2::annotate(
      "text",
      x = anno_x_pos * l$x_max,
      y = anno_y_pos * l$y_max,
      label = lab,
      size = anno_size
    )
  } else {
    p
  }
}
