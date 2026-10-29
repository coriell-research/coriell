make_mat <- function() {
  matrix(
    c(1, 2, 3, 4, 5, 2, 4, 6, 8, 100, 0, 0, 1, 1, 2, 10, 9, 8, 7, 6),
    nrow = 5,
    ncol = 4,
    dimnames = list(paste0("g", 1:5), c("s3", "s1", "s4", "s2"))
  )
}

# Mock the graphics primitives used by plot_boxplot2 and capture their arguments
local_capture_bxp <- function(cap, env = parent.frame()) {
  testthat::local_mocked_bindings(
    bxp = function(z, ...) {
      cap$z <- z
      cap$dots <- list(...)
      invisible(NULL)
    },
    legend = function(...) {
      cap$legend <- list(...)
      invisible(NULL)
    },
    abline = function(...) {
      cap$abline <- list(...)
      invisible(NULL)
    },
    .package = "graphics",
    .env = env
  )
}

# plot_boxplot ----------------------------------------------------------------

test_that("plot_boxplot box statistics match grDevices::boxplot.stats per sample", {
  x <- make_mat()
  p <- plot_boxplot(x)
  expect_s3_class(p, "ggplot")
  b <- ggplot2::ggplot_build(p)$data[[1]]
  lv <- levels(p$data$.sample)
  for (i in seq_along(lv)) {
    ref <- grDevices::boxplot.stats(unname(x[, lv[i]]))
    expect_equal(unname(c(b$ymin[i], b$lower[i], b$middle[i], b$upper[i], b$ymax[i])), ref$stats)
    expect_equal(as.numeric(unlist(b$outliers[i])), ref$out)
  }
})

test_that("plot_boxplot long data is a lossless reshape of the matrix", {
  x <- make_mat()
  d <- plot_boxplot(x)$data
  expect_equal(nrow(d), length(x))
  expect_equal(d$.value, unname(x[cbind(d$.feature, as.character(d$.sample))]))
})

test_that("plot_boxplot reference line is the global median", {
  x <- make_mat()
  x[2, 3] <- NA
  p <- suppressWarnings(plot_boxplot(x))
  b <- suppressWarnings(ggplot2::ggplot_build(p))
  expect_equal(b$data[[2]]$yintercept, median(x, na.rm = TRUE))
})

test_that("plot_boxplot rle = TRUE subtracts NA-aware row medians", {
  x <- make_mat()
  x[1, 2] <- NA
  p <- suppressWarnings(plot_boxplot(x, rle = TRUE))
  d <- p$data
  rm <- apply(x, 1, median, na.rm = TRUE)
  expect_equal(d$.value, unname((x - rm)[cbind(d$.feature, as.character(d$.sample))]))
  expect_equal(unname(c(tapply(d$.value, d$.feature, median, na.rm = TRUE))), rep(0, nrow(x)))
  expect_equal(suppressWarnings(ggplot2::ggplot_build(p))$data[[2]]$yintercept, 0)
})

test_that("plot_boxplot groups samples by fillBy and keeps within-group input order", {
  x <- make_mat()
  colnames(x) <- c("z", "y", "x", "w")
  meta <- data.frame(G = c("a", "b", "a", "b"), row.names = colnames(x))
  p <- plot_boxplot(x, meta, fillBy = "G")
  expect_equal(levels(p$data$.sample), c("z", "x", "y", "w"))
  per_sample <- unique(as.data.frame(p$data)[, c(".sample", "G")])
  expect_equal(per_sample$G[match(colnames(x), per_sample$.sample)], meta$G)
  b <- ggplot2::ggplot_build(p)$data[[1]]
  expect_equal(length(unique(b$fill)), 2)
  expect_equal(b$fill[1], b$fill[2])
  expect_false(b$fill[2] == b$fill[3])
  expect_equal(p$labels$fill, "G")
})

test_that("plot_boxplot default fill is one colour per sample labelled 'Sample'", {
  p <- plot_boxplot(make_mat())
  expect_equal(p$labels$fill, "Sample")
  expect_equal(length(unique(ggplot2::ggplot_build(p)$data[[1]]$fill)), 4)
})

test_that("plot_boxplot violin = TRUE uses geom_violin and forwards ...", {
  p <- plot_boxplot(make_mat(), violin = TRUE, alpha = 0.5)
  expect_s3_class(p$layers[[1]]$geom, "GeomViolin")
  expect_equal(p$layers[[1]]$aes_params$alpha, 0.5)
  expect_s3_class(plot_boxplot(make_mat())$layers[[1]]$geom, "GeomBoxplot")
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("plot_boxplot rejects metadata whose rows are not in column order", {
  x <- make_mat()
  meta <- data.frame(G = c("a", "b", "a", "b"), row.names = rev(colnames(x)))
  expect_error(plot_boxplot(x, meta, fillBy = "G"), "do not match")
})

test_that("plot_boxplot data.frame and SummarizedExperiment methods equal the matrix method", {
  skip_if_not_installed("SummarizedExperiment")
  x <- make_mat()
  meta <- data.frame(G = c("a", "b", "a", "b"), row.names = colnames(x))
  ref <- plot_boxplot(x, meta, fillBy = "G", rle = TRUE)
  p_df <- plot_boxplot(as.data.frame(x), meta, fillBy = "G", rle = TRUE)
  se <- SummarizedExperiment::SummarizedExperiment(list(counts = x), colData = meta)
  p_se <- plot_boxplot(se, fillBy = "G", rle = TRUE)
  expect_equal(as.data.frame(p_df$data), as.data.frame(ref$data))
  expect_equal(as.data.frame(p_se$data), as.data.frame(ref$data))
  expect_equal(
    ggplot2::ggplot_build(p_se)$data[[1]]$middle,
    ggplot2::ggplot_build(ref)$data[[1]]$middle
  )
})

test_that("plot_boxplot works for matrices without rownames", {
  x <- matrix(as.numeric(1:20), 5, 4, dimnames = list(NULL, paste0("s", 1:4)))
  p <- plot_boxplot(x)
  expect_s3_class(p, "ggplot")
  d <- p$data
  expect_setequal(d$.feature, as.character(1:5))
  expect_equal(d$.value, x[cbind(as.integer(d$.feature), match(d$.sample, colnames(x)))])
  b <- ggplot2::ggplot_build(p)$data[[1]]
  expect_equal(b$middle, unname(apply(x, 2, median)))
})

# plot_boxplot2 ---------------------------------------------------------------

test_that("plot_boxplot2 passes statistics identical to graphics::boxplot", {
  set.seed(42)
  x <- matrix(rnorm(200), 50, 4, dimnames = list(NULL, c("s3", "s1", "s4", "s2")))
  x[c(3, 7), 2] <- NA
  x[1, 3] <- 10
  x[2, 3] <- -9
  x[5, 1] <- 8
  cap <- new.env()
  local_capture_bxp(cap)
  plot_boxplot2(x)
  ref <- graphics::boxplot(x, plot = FALSE)
  expect_equal(unname(cap$z$stats), ref$stats)
  expect_equal(unname(cap$z$n), ref$n)
  expect_equal(unname(as.vector(cap$z$out)), ref$out)
  expect_equal(as.numeric(cap$z$group), ref$group)
  expect_equal(cap$z$names, colnames(x))
  expect_equal(cap$dots$boxfill, "lightgrey")
  expect_null(cap$legend)
  expect_null(cap$abline)
})

test_that("plot_boxplot2 outlier groups are correct with equal outlier counts", {
  x <- cbind(a = c(1:9, 100, -100), b = c(1:9, 50, -50))
  cap <- new.env()
  local_capture_bxp(cap)
  plot_boxplot2(x)
  ref <- graphics::boxplot(x, plot = FALSE)
  expect_equal(as.vector(cap$z$out), ref$out)
  expect_equal(as.numeric(cap$z$group), ref$group)
})

test_that("plot_boxplot2 rle, fill colours, legend and zero line", {
  set.seed(42)
  x <- matrix(rnorm(200), 50, 4, dimnames = list(NULL, c("s3", "s1", "s4", "s2")))
  meta <- data.frame(G = c("b", "a", "b", "c"), row.names = colnames(x))
  cap <- new.env()
  local_capture_bxp(cap)
  plot_boxplot2(x, meta, fill_by = "G", rle = TRUE)
  ref <- graphics::boxplot(x - apply(x, 1, median), plot = FALSE)
  expect_equal(unname(cap$z$stats), ref$stats)
  pal <- grDevices::hcl.colors(3, "Dark 3")
  expect_equal(unname(cap$dots$boxfill), pal[c(2, 1, 2, 3)])
  expect_equal(cap$legend$legend, c("a", "b", "c"))
  expect_equal(unname(cap$legend$fill), pal)
  expect_equal(cap$abline$h, 0)
})

test_that("plot_boxplot2 show_legend = FALSE suppresses the legend", {
  x <- matrix(as.numeric(1:20), 5, 4, dimnames = list(NULL, letters[1:4]))
  meta <- data.frame(G = c("a", "a", "b", "b"), row.names = colnames(x))
  cap <- new.env()
  local_capture_bxp(cap)
  plot_boxplot2(x, meta, fill_by = "G", show_legend = FALSE)
  expect_null(cap$legend)
})

test_that("plot_boxplot2 rle works for DelayedMatrix input", {
  skip_if_not_installed("DelayedArray")
  skip_if_not_installed("DelayedMatrixStats")
  set.seed(1)
  m <- matrix(rnorm(60), 20, 3, dimnames = list(NULL, c("a", "b", "c")))
  cap <- new.env()
  local_capture_bxp(cap)
  plot_boxplot2(DelayedArray::DelayedArray(m), rle = TRUE)
  ref <- graphics::boxplot(m - apply(m, 1, median), plot = FALSE)
  expect_equal(unname(cap$z$stats), ref$stats)
})

test_that("plot_boxplot2 renders on a real device", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot_boxplot2(GSE161650_lc))
})
