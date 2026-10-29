make_dens_mat <- function() {
  set.seed(1)
  x <- matrix(
    rnorm(400), 100, 4,
    dimnames = list(paste0("g", 1:100), c("s3", "s1", "s4", "s2"))
  )
  x[, 2] <- x[, 2] + 3
  x
}

# Mock graphics::lines and graphics::legend and capture their arguments
local_capture_lines <- function(cap, env = parent.frame()) {
  cap$lines <- list()
  testthat::local_mocked_bindings(
    lines = function(x, y, ...) {
      cap$lines[[length(cap$lines) + 1]] <- c(list(x = x, y = y), list(...))
      invisible(NULL)
    },
    legend = function(...) {
      cap$legend <- list(...)
      invisible(NULL)
    },
    .package = "graphics",
    .env = env
  )
}

# plot_density ----------------------------------------------------------------

test_that("plot_density curves equal stats::density over the pooled range", {
  x <- make_dens_mat()
  p <- plot_density(x)
  expect_s3_class(p, "ggplot")
  d <- ggplot2::ggplot_build(p)$data[[1]]
  rng <- range(x)
  smp <- sort(colnames(x))
  for (g in seq_along(smp)) {
    ref <- stats::density(x[, smp[g]], from = rng[1], to = rng[2], n = 512)
    dd <- d[d$group == g, ]
    expect_equal(dd$x, ref$x)
    expect_equal(dd$y, ref$y)
  }
  expect_equal(p$labels$colour, "Sample")
  expect_equal(p$labels$y, "Density")
})

test_that("plot_density colBy draws one curve per sample, coloured by group", {
  x <- make_dens_mat()
  meta <- data.frame(G = c("b", "a", "b", "a"), row.names = colnames(x))
  p <- plot_density(x, meta, colBy = "G")
  d <- ggplot2::ggplot_build(p)$data[[1]]
  expect_equal(length(unique(d$group)), 4)
  cg <- tapply(d$colour, d$group, unique)
  expect_equal(unname(cg[1]), unname(cg[2]))
  expect_equal(unname(cg[3]), unname(cg[4]))
  expect_false(cg[1] == cg[3])
  expect_equal(p$labels$colour, "G")
})

test_that("plot_density data.frame and SummarizedExperiment methods equal the matrix method", {
  skip_if_not_installed("SummarizedExperiment")
  set.seed(1)
  x <- matrix(rnorm(200), 50, 4, dimnames = list(paste0("g", 1:50), paste0("s", 1:4)))
  meta <- data.frame(G = c("a", "a", "b", "b"), row.names = colnames(x))
  ref <- ggplot2::ggplot_build(plot_density(x, meta, colBy = "G"))$data[[1]]
  se <- SummarizedExperiment::SummarizedExperiment(list(counts = x), colData = meta)
  expect_equal(ggplot2::ggplot_build(plot_density(se, colBy = "G"))$data[[1]]$y, ref$y)
  expect_equal(
    ggplot2::ggplot_build(plot_density(as.data.frame(x), meta, colBy = "G"))$data[[1]]$y,
    ref$y
  )
})

test_that("plot_density works for matrices without rownames", {
  x <- unname(make_dens_mat())
  colnames(x) <- paste0("s", 1:4)
  p <- plot_density(x)
  expect_s3_class(p, "ggplot")
  d <- ggplot2::ggplot_build(p)$data[[1]]
  rng <- range(x)
  for (g in 1:4) {
    ref <- stats::density(x[, g], from = rng[1], to = rng[2], n = 512)
    expect_equal(d$y[d$group == g], ref$y)
  }
})

# plot_density2 ---------------------------------------------------------------

test_that("plot_density2 curves equal stats::density on the common range", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  set.seed(2)
  x <- matrix(rnorm(300), 100, 3, dimnames = list(NULL, c("s3", "s1", "s2")))
  x[, 2] <- x[, 2] * 2 + 1
  x[4, 1] <- NA
  cap <- new.env()
  local_capture_lines(cap)
  plot_density2(x)
  expect_length(cap$lines, 3)
  rng <- range(x, na.rm = TRUE)
  for (i in 1:3) {
    ref <- stats::density(x[!is.na(x[, i]), i], from = rng[1], to = rng[2])
    expect_equal(cap$lines[[i]]$x, ref$x)
    expect_equal(cap$lines[[i]]$y, ref$y)
    expect_equal(cap$lines[[i]]$col, "black")
  }
  expect_null(cap$legend)
})

test_that("plot_density2 colours lines by col_by and passes line parameters", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  set.seed(2)
  x <- matrix(rnorm(300), 100, 3, dimnames = list(NULL, c("s3", "s1", "s2")))
  meta <- data.frame(G = c("b", "a", "b"), row.names = colnames(x))
  cap <- new.env()
  local_capture_lines(cap)
  plot_density2(x, meta, col_by = "G", line_weight = 3, line_type = 2)
  pal <- grDevices::hcl.colors(2, "Dark 3")
  expect_equal(unname(vapply(cap$lines, function(l) unname(l$col), "")), pal[c(2, 1, 2)])
  expect_true(all(vapply(cap$lines, function(l) l$lwd == 3 && l$lty == 2, logical(1))))
  expect_equal(cap$legend$legend, c("a", "b"))
  expect_equal(unname(cap$legend$col), pal)
})

test_that("plot_density2 show_legend = FALSE suppresses the legend", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  set.seed(3)
  x <- matrix(rnorm(60), 20, 3, dimnames = list(NULL, c("a", "b", "c")))
  meta <- data.frame(G = c("x", "y", "x"), row.names = colnames(x))
  cap <- new.env()
  local_capture_lines(cap)
  plot_density2(x, meta, col_by = "G", show_legend = FALSE)
  expect_length(cap$lines, 3)
  expect_null(cap$legend)
})

test_that("plot_density2 renders on a real device", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot_density2(GSE161650_lc))
})
