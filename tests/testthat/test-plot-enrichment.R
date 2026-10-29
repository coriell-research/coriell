stats5 <- c(a = 5, b = 4, c = 3, d = 2, e = 1)

test_that("plot_enrichment curve matches a hand-computed GSEA walk", {
  skip_if_not_installed("fgsea")
  p <- plot_enrichment(c("a", "c"), stats5)
  expect_s3_class(p, "ggplot")
  curve <- ggplot2::layer_data(p, 1)
  es_at <- function(r) max(curve$y[curve$x == r])
  # hits weigh 5/8 and 3/8, each of the 3 misses weighs 1/3
  expect_equal(es_at(1), 5 / 8)
  expect_equal(es_at(2), 5 / 8 - 1 / 3)
  expect_equal(es_at(3), 2 / 3)
  expect_equal(max(curve$y), 2 / 3)
  expect_equal(min(curve$y), 0)
  expect_equal(unique(curve$colour), "red2")
})

test_that("plot_enrichment with gseaParam = 0 gives the unweighted walk", {
  skip_if_not_installed("fgsea")
  curve <- ggplot2::layer_data(plot_enrichment(c("a", "c"), stats5, gseaParam = 0), 1)
  expect_equal(max(curve$y[curve$x == 1]), 1 / 2)
  expect_equal(max(curve$y[curve$x == 3]), 1 / 2 - 1 / 3 + 1 / 2)
})

test_that("plot_enrichment colors down sets blue unless lineColor is given", {
  skip_if_not_installed("fgsea")
  curve <- ggplot2::layer_data(plot_enrichment(c("d", "e"), stats5), 1)
  expect_equal(min(curve$y), -1)
  expect_equal(unique(curve$colour), "blue2")
  curve <- ggplot2::layer_data(
    plot_enrichment(c("d", "e"), stats5, lineColor = "green"),
    1
  )
  expect_equal(unique(curve$colour), "green")
})

test_that("plot_enrichment ticks sit at the pathway gene ranks", {
  skip_if_not_installed("fgsea")
  set.seed(11)
  stats <- setNames(rnorm(100), paste0("g", 1:100))
  pw <- paste0("g", c(3, 17, 42, 88))
  p <- plot_enrichment(pw, stats)
  ticks <- ggplot2::layer_data(p, 2)
  curve <- ggplot2::layer_data(p, 1)
  expect_setequal(ticks$x, match(pw, names(sort(stats, decreasing = TRUE))))
  expect_equal(ticks$xend, ticks$x)
  spread <- max(curve$y) - min(curve$y)
  expect_equal(unique(ticks$yend), spread / 20)
  expect_equal(unique(ticks$y), -spread / 20)
})

test_that("plot_enrichment draws max, min and zero reference lines", {
  skip_if_not_installed("fgsea")
  set.seed(12)
  stats <- setNames(rnorm(200), paste0("g", 1:200))
  pw <- sample(names(stats), 20)
  p <- plot_enrichment(pw, stats)
  curve <- ggplot2::layer_data(p, 1)
  expect_equal(ggplot2::layer_data(p, 3)$yintercept, max(curve$y))
  expect_equal(ggplot2::layer_data(p, 4)$yintercept, min(curve$y))
  expect_equal(ggplot2::layer_data(p, 5)$yintercept, 0)
})

test_that("plot_enrichment extreme ES agrees with fgsea::calcGseaStat", {
  skip_if_not_installed("fgsea")
  examplePathways <- fgsea::examplePathways
  exampleRanks <- fgsea::exampleRanks
  pw <- examplePathways[["5991130_Programmed_Cell_Death"]]
  curve <- ggplot2::layer_data(plot_enrichment(pw, exampleRanks), 1)
  sr <- sort(exampleRanks, decreasing = TRUE)
  es <- fgsea::calcGseaStat(sr, na.omit(match(pw, names(sr))))
  ext <- if (abs(min(curve$y)) > max(curve$y)) min(curve$y) else max(curve$y)
  expect_equal(ext, es)
})

test_that("plot_enrichment annotation label and placement", {
  skip_if_not_installed("fgsea")
  p <- plot_enrichment(
    c("a", "c"),
    stats5,
    annotate = TRUE,
    FDR = 0.012345,
    ES = 0.66667,
    NES = 1.23456,
    anno_x_pos = 0.5,
    anno_y_pos = 0.5
  )
  ann <- ggplot2::layer_data(p, 6)
  expect_equal(ann$label, "FDR: 1.23e-02 \nES: 0.67 \nNES: 1.23")
  # x_max is N + 1 = 6; the y limit is the max ES of 2/3
  expect_equal(ann$x, 0.5 * 6)
  expect_equal(ann$y, 0.5 * 2 / 3)

  pd <- plot_enrichment(
    c("d", "e"),
    stats5,
    annotate = TRUE,
    FDR = 0.5,
    ES = -1,
    NES = -2,
    anno_y_pos = 0.5
  )
  expect_equal(ggplot2::layer_data(pd, 6)$y, -0.5)
  expect_error(
    plot_enrichment(c("a", "c"), stats5, annotate = TRUE, ES = 1, NES = 1),
    "FDR"
  )
})

test_that("plot_enrichment applies title and axis labels", {
  skip_if_not_installed("fgsea")
  p <- plot_enrichment(c("a", "c"), stats5, plotTitle = "PW", xlab = "R", ylab = "E")
  labs <- ggplot2::get_labs(p)
  expect_equal(labs$title, "PW")
  expect_equal(labs$x, "R")
  expect_equal(labs$y, "E")
})
