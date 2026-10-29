df <- data.frame(
  UMAP1 = c(-1, 0, 1, 2),
  UMAP2 = c(3, -2, 0, 1),
  Group = c("A", "A", "B", "B"),
  Batch = c("x", "y", "x", "y"),
  score = 1:4
)

test_that("plot_umap maps the embedding to x/y with default reference lines", {
  p <- plot_umap(df)
  expect_s3_class(p, "ggplot")
  ld <- ggplot2::layer_data(p, 1)
  expect_equal(nrow(ld), 4)
  expect_equal(ld$x, df$UMAP1)
  expect_equal(ld$y, df$UMAP2)
  expect_length(unique(ld$colour), 1)
  expect_length(unique(ld$shape), 1)
  expect_equal(ggplot2::layer_data(p, 2)$yintercept, 0)
  expect_equal(ggplot2::layer_data(p, 3)$xintercept, 0)
  labs <- ggplot2::get_labs(p)
  expect_equal(labs$x, "UMAP1")
  expect_equal(labs$y, "UMAP2")
  expect_null(labs$colour)
  expect_null(labs$shape)
})

test_that("plot_umap colour and shape partitions follow colBy and shapeBy", {
  p <- plot_umap(df, colBy = "Group", shapeBy = "Batch")
  ld <- ggplot2::layer_data(p, 1)
  expect_equal(
    as.integer(factor(ld$colour, levels = unique(ld$colour))),
    as.integer(factor(df$Group))
  )
  expect_equal(
    as.integer(factor(ld$shape, levels = unique(ld$shape))),
    as.integer(factor(df$Batch))
  )
  labs <- ggplot2::get_labs(p)
  expect_equal(labs$colour, "Group")
  expect_equal(labs$shape, "Batch")
})

test_that("plot_umap respects x/y and gives a continuous scale for numeric colBy", {
  p <- plot_umap(df, x = "UMAP2", y = "UMAP1", colBy = "score")
  ld <- ggplot2::layer_data(p, 1)
  expect_equal(ld$x, df$UMAP2)
  expect_equal(ld$y, df$UMAP1)
  expect_length(unique(ld$colour), 4)
  scale <- ggplot2::ggplot_build(p)$plot$scales$get_scales("colour")
  expect_s3_class(scale, "ScaleContinuous")
})

test_that("plot_umap applies point and reference line settings", {
  p <- plot_umap(
    df,
    pointSize = 5,
    pointAlpha = 0.4,
    hline = 0.5,
    vline = -0.5,
    hlineType = 3,
    vlineType = 1
  )
  ld <- ggplot2::layer_data(p, 1)
  expect_true(all(ld$size == 5))
  expect_true(all(ld$alpha == 0.4))
  h <- ggplot2::layer_data(p, 2)
  v <- ggplot2::layer_data(p, 3)
  expect_equal(h$yintercept, 0.5)
  expect_equal(v$xintercept, -0.5)
  expect_equal(h$linetype, 3)
  expect_equal(v$linetype, 1)
})

test_that("plot_umap works on the output of UMAP()", {
  skip_if_not_installed("umap")
  skip_if_not_installed("PCAtools")
  set.seed(123)
  md <- data.frame(
    row.names = colnames(GSE161650_lc),
    Group = rep(c("DMSO", "THZ1"), each = 3)
  )
  pc <- PCAtools::pca(GSE161650_lc, md, center = TRUE, scale = TRUE)
  u <- suppressWarnings(UMAP(pc, n_neighbors = 2))
  ld <- ggplot2::layer_data(plot_umap(u, colBy = "Group"), 1)
  expect_equal(nrow(ld), 6)
  expect_equal(ld$x, u$UMAP1)
  expect_equal(ld$y, u$UMAP2)
  expect_length(unique(ld$colour), 2)
})

test_that("plot_umap errors on unknown colBy or shapeBy columns", {
  expect_error(plot_umap(df, colBy = "nope"), "nope is not a column")
  expect_error(plot_umap(df, shapeBy = "nah"), "nah is not a column")
})
