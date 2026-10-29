# Record the base-graphics display list as (name, args) pairs
record_base_plot <- function(expr) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  grDevices::dev.control("enable")
  force(expr)
  dl <- grDevices::recordPlot()[[1]]
  lapply(dl, function(e) {
    f <- e[[2]][[1]]
    list(name = if (is.list(f)) f$name else NA_character_, args = e[[2]][-1])
  })
}
calls_named <- function(calls, nm) Filter(function(e) identical(e$name, nm), calls)
# Correlation labels: all text() labels except the diagonal column names
cor_labels <- function(calls, cn) {
  labs <- unlist(lapply(calls_named(calls, "C_text"), function(e) e$args[[2]]))
  trimws(labs[!labs %in% cn])
}

test_that("plot_cor_pairs shows the Pearson correlation for every pair", {
  skip_if_not_installed("KernSmooth")
  v <- 1:20
  m <- cbind(A = v, B = -v, C = v^3)
  labs <- cor_labels(record_base_plot(plot_cor_pairs(m)), colnames(m))
  r_ac <- cor(v, v^3)
  expected <- sprintf("%.2f", c(-1, r_ac, -r_ac))
  expect_equal(sort(labs), sort(expected))
})

test_that("plot_cor_pairs respects cor_method", {
  skip_if_not_installed("KernSmooth")
  v <- 1:20
  m <- cbind(A = v, B = -v, C = v^3)
  for (meth in c("spearman", "kendall")) {
    labs <- cor_labels(
      record_base_plot(plot_cor_pairs(m, cor_method = meth)),
      colnames(m)
    )
    expect_equal(sort(labs), c("-1.00", "-1.00", "1.00"), info = meth)
  }
  labs <- cor_labels(record_base_plot(plot_cor_pairs(m)), colnames(m))
  expect_false(all(labs %in% c("-1.00", "1.00")))
})

test_that("plot_cor_pairs passes cor_use through to cor()", {
  skip_if_not_installed("KernSmooth")
  set.seed(1)
  m <- cbind(A = rnorm(50), B = rnorm(50), C = rnorm(50))
  m[c(3, 7), "A"] <- NA
  labs <- cor_labels(
    record_base_plot(plot_cor_pairs(m, cor_use = "pairwise.complete.obs")),
    colnames(m)
  )
  r <- cor(m, use = "pairwise.complete.obs")
  expected <- sprintf("%.2f", c(r["A", "B"], r["A", "C"], r["B", "C"]))
  expect_equal(sort(labs), sort(expected))

  labs <- cor_labels(record_base_plot(plot_cor_pairs(m)), colnames(m))
  expect_equal(sum(labs == "NA"), 2)
})

test_that("plot_cor_pairs formats cor_digits as decimal places", {
  skip_if_not_installed("KernSmooth")
  set.seed(4)
  m <- cbind(A = rnorm(100), B = rnorm(100), C = rnorm(100))
  r <- cor(m)
  pairs_r <- c(r["A", "B"], r["A", "C"], r["B", "C"])
  for (d in c(1, 3)) {
    labs <- cor_labels(
      record_base_plot(plot_cor_pairs(m, cor_digits = d)),
      colnames(m)
    )
    expect_equal(sort(labs), sort(sprintf(paste0("%.", d, "f"), pairs_r)), info = d)
    expect_true(all(nchar(sub("^-?[0-9]+\\.", "", labs)) == d))
  }
})

test_that("plot_cor_pairs diagonal histograms are scaled hist() counts", {
  skip_if_not_installed("KernSmooth")
  set.seed(42)
  m <- cbind(A = rnorm(200), B = rexp(200), C = runif(200))
  rects <- calls_named(
    record_base_plot(plot_cor_pairs(m, hist_breaks = 15, hist_col = "orange")),
    "C_rect"
  )
  expect_length(rects, 3)
  for (j in 1:3) {
    h <- hist(m[, j], breaks = 15, plot = FALSE)
    a <- rects[[j]]$args
    expect_equal(a[[1]], head(h$breaks, -1))
    expect_equal(a[[3]], h$breaks[-1])
    expect_equal(a[[4]], h$counts / max(h$counts))
    expect_equal(a$col, "orange")
  }
})

test_that("plot_cor_pairs upper panels have a smoothScatter, identity line and lowess fit", {
  skip_if_not_installed("KernSmooth")
  set.seed(2)
  m <- cbind(A = rnorm(100), B = rnorm(100), C = rnorm(100))
  calls <- record_base_plot(plot_cor_pairs(m, scatter_ab_col = "blue"))
  ab <- calls_named(calls, "C_abline")
  expect_length(ab, 3)
  for (e in ab) {
    expect_equal(e$args[[1]], 0)
    expect_equal(e$args[[2]], 1)
    expect_equal(e$args[[6]], "blue")
  }
  expect_length(calls_named(calls, "C_image"), 3)

  xy <- lapply(calls_named(calls, "C_plotXY"), function(e) e$args[[1]])
  lw <- lowess(m[, "B"], m[, "A"])
  matches <- vapply(
    xy,
    function(p) isTRUE(all.equal(p$x, lw$x)) && isTRUE(all.equal(p$y, lw$y)),
    logical(1)
  )
  expect_true(any(matches))
})

test_that("plot_cor_pairs uses the supplied scatter_colramp", {
  skip_if_not_installed("KernSmooth")
  called <- FALSE
  ramp <- function(n) {
    called <<- TRUE
    grDevices::grey.colors(n)
  }
  set.seed(3)
  m <- cbind(A = rnorm(50), B = rnorm(50))
  calls <- record_base_plot(plot_cor_pairs(m, scatter_colramp = ramp))
  expect_true(called)
  img_cols <- calls_named(calls, "C_image")[[1]]$args[[4]]
  expect_true(all(img_cols %in% grDevices::grey.colors(256)))
})
