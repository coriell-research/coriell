par_mat <- function() {
  set.seed(3)
  sds <- c(5, 0.1, 3, 0.2, 1, 0.3, 4, 0.5, 2, 0.4)
  matrix(
    rnorm(40, sd = rep(sds, 4)), 10, 4,
    dimnames = list(paste0("g", 1:10), c("z", "y", "x", "w"))
  )
}

par2_mat <- function() {
  set.seed(3)
  sds <- c(5, 0.1, 3, 0.2, 1, 0.3, 4, 0.5, 2, 0.4)
  matrix(
    rnorm(60, sd = rep(sds, 6)), 10, 6,
    dimnames = list(paste0("g", 1:10), paste0("s", 1:6))
  )
}

# Mock the graphics primitives used by plot_parallel2 and capture their arguments
local_capture_par2 <- function(cap, env = parent.frame()) {
  cap$lines <- list()
  cap$axis <- list()
  testthat::local_mocked_bindings(
    lines = function(x, y, ...) {
      cap$lines[[length(cap$lines) + 1]] <- c(list(x = x, y = y), list(...))
      invisible(NULL)
    },
    axis = function(side, at = NULL, labels = TRUE, ...) {
      cap$axis[[side]] <- list(at = at, labels = labels)
      invisible(NULL)
    },
    title = function(main = NULL, sub = NULL, xlab = NULL, ylab = NULL, ...) {
      cap$title <- list(main = main, xlab = xlab, ylab = ylab)
      invisible(NULL)
    },
    .package = "graphics",
    .env = env
  )
}

# plot_parallel ---------------------------------------------------------------

test_that("plot_parallel draws one line per feature through the exact values", {
  x <- par_mat()
  p <- plot_parallel(x)
  expect_s3_class(p, "ggplot")
  d <- p$data
  expect_equal(nrow(d), length(x))
  expect_equal(d$.value, unname(x[cbind(d$.feature, as.character(d$.sample))]))
  b <- ggplot2::ggplot_build(p)$data[[1]]
  expect_equal(length(unique(b$group)), nrow(x))
  expect_equal(as.vector(table(b$group)), rep(ncol(x), nrow(x)))
  expect_equal(p$labels$colour, "Sample")
})

test_that("plot_parallel colBy orders the x-axis by group and colours by group", {
  x <- par_mat()
  meta <- data.frame(G = c("a", "b", "a", "b"), row.names = colnames(x))
  p <- plot_parallel(x, meta, colBy = "G")
  expect_equal(levels(p$data$.sample), c("z", "x", "y", "w"))
  b <- ggplot2::ggplot_build(p)$data[[1]]
  expect_equal(as.vector(table(b$colour)), c(20L, 20L))
  expect_equal(p$labels$colour, "G")
})

test_that("plot_parallel removeVar keeps the highest-variance features", {
  x <- par_mat()
  p <- suppressMessages(plot_parallel(x, removeVar = 0.5))
  v <- apply(x, 1, var)
  expect_setequal(unique(p$data$.feature), names(sort(v, decreasing = TRUE))[1:5])
})

test_that("plot_parallel forwards ... and the data.frame method matches", {
  x <- par_mat()
  expect_equal(plot_parallel(x, alpha = 0.3)$layers[[1]]$aes_params$alpha, 0.3)
  expect_equal(
    as.data.frame(plot_parallel(as.data.frame(x))$data),
    as.data.frame(plot_parallel(x)$data)
  )
})

test_that("plot_parallel SummarizedExperiment method uses colData as metadata", {
  skip_if_not_installed("SummarizedExperiment")
  x <- par_mat()
  meta <- data.frame(G = c("a", "b", "a", "b"), row.names = colnames(x))
  se <- SummarizedExperiment::SummarizedExperiment(list(counts = x), colData = meta)
  p <- plot_parallel(se, colBy = "G")
  d <- as.data.frame(p$data)
  expect_equal(d$G, meta[as.character(d$.sample), "G"])
  expect_equal(levels(p$data$.sample), c("z", "x", "y", "w"))
  expect_equal(d$.value, unname(x[cbind(d$.feature, as.character(d$.sample))]))
  expect_equal(p$labels$colour, "G")
  expect_s3_class(plot_parallel(se), "ggplot")
})

test_that("plot_parallel works for matrices without rownames", {
  x <- unname(par_mat())
  colnames(x) <- paste0("s", 1:4)
  p <- plot_parallel(x)
  expect_s3_class(p, "ggplot")
  d <- p$data
  expect_setequal(d$.feature, as.character(1:10))
  expect_equal(d$.value, x[cbind(as.integer(d$.feature), match(d$.sample, colnames(x)))])
  b <- ggplot2::ggplot_build(p)$data[[1]]
  expect_equal(length(unique(b$group)), nrow(x))
})

# plot_parallel2 --------------------------------------------------------------

test_that("plot_parallel2 draws z-scored top-variance rows", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  x <- par2_mat()
  cap <- new.env()
  local_capture_par2(cap)
  suppressMessages(plot_parallel2(
    x, remove_var = 0.5, color = "red", alpha = 0.25,
    plot_title = "T", x_label = "X", y_label = "Y"
  ))
  keep <- names(sort(apply(x, 1, var), decreasing = TRUE))[1:5]
  z <- (x[keep, ] - rowMeans(x[keep, ])) / apply(x[keep, ], 1, sd)
  expect_length(cap$lines, 5)
  for (i in 1:5) {
    expect_equal(cap$lines[[i]]$x, 1:6)
    expect_equal(cap$lines[[i]]$y, unname(z[i, ]))
    expect_equal(mean(cap$lines[[i]]$y), 0, tolerance = 1e-12)
    expect_equal(sd(cap$lines[[i]]$y), 1, tolerance = 1e-12)
    expect_equal(cap$lines[[i]]$col, grDevices::adjustcolor("red", alpha.f = 0.25))
  }
  expect_equal(cap$axis[[1]]$labels, colnames(x))
  expect_equal(cap$title, list(main = "T", xlab = "X", ylab = "Y"))
})

test_that("plot_parallel2 scale = FALSE draws raw values", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  x <- par2_mat()
  cap <- new.env()
  local_capture_par2(cap)
  suppressMessages(plot_parallel2(x, remove_var = 0.5, scale = FALSE))
  keep <- names(sort(apply(x, 1, var), decreasing = TRUE))[1:5]
  expect_length(cap$lines, 5)
  for (i in 1:5) expect_equal(cap$lines[[i]]$y, unname(x[keep[i], ]))
})

test_that("plot_parallel2 renders on a real device", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(suppressMessages(plot_parallel2(GSE161650_lc, remove_var = 0.5)))
})
