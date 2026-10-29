x3 <- matrix(
  c(0, 0, 3, 4, 6, 8),
  nrow = 2,
  dimnames = list(c("g1", "g2"), c("s1", "s2", "s3"))
)

# Replace pheatmap::pheatmap with a stub that records its arguments
capture_pheatmap <- function(env = parent.frame()) {
  cap <- new.env()
  local_mocked_bindings(
    pheatmap = function(...) cap$args <- list(...),
    .package = "pheatmap",
    .env = env
  )
  cap
}

test_that("plot_dist passes the hand-computed euclidean distances to pheatmap", {
  skip_if_not_installed("pheatmap")
  skip_if_not_installed("viridisLite")
  cap <- capture_pheatmap()
  plot_dist(x3)
  ids <- c("s1", "s2", "s3")
  expected <- matrix(
    c(0, 5, 10, 5, 0, 5, 10, 5, 0),
    3,
    dimnames = list(ids, ids)
  )
  expect_equal(cap$args$mat, expected)
})

test_that("plot_dist respects method and matches stats::dist", {
  skip_if_not_installed("pheatmap")
  skip_if_not_installed("viridisLite")
  cap <- capture_pheatmap()
  plot_dist(x3, method = "manhattan")
  expect_equal(unname(cap$args$mat[1, ]), c(0, 7, 14))
  plot_dist(x3, method = "maximum")
  expect_equal(unname(cap$args$mat[1, ]), c(0, 4, 8))

  set.seed(5)
  m <- matrix(rnorm(40), 8, dimnames = list(NULL, letters[1:5]))
  for (meth in c("euclidean", "maximum", "manhattan", "canberra", "minkowski")) {
    plot_dist(m, method = meth)
    expect_equal(cap$args$mat, as.matrix(dist(t(m), method = meth)), info = meth)
  }
})

test_that("plot_dist supports method = 'binary'", {
  skip_if_not_installed("pheatmap")
  skip_if_not_installed("viridisLite")
  cap <- capture_pheatmap()
  set.seed(5)
  mb <- matrix(rbinom(40, 1, 0.5), 8, dimnames = list(NULL, letters[1:5]))
  plot_dist(mb, method = "binary")
  expect_equal(cap$args$mat, as.matrix(dist(t(mb), method = "binary")))
})

test_that("plot_dist canberra matches stats::dist when there are zeros", {
  skip_if_not_installed("pheatmap")
  skip_if_not_installed("viridisLite")
  cap <- capture_pheatmap()
  m <- cbind(a = c(0, 1, 2, 0), b = c(0, 3, 0, 1), c = c(1, 1, 1, 1))
  plot_dist(m, method = "canberra")
  # a vs b: the first row is 0/0 and is dropped, so 3 terms rescaled by 4/3
  expect_equal(cap$args$mat["a", "b"], (2 / 4 + 2 / 2 + 1 / 1) * 4 / 3)
  expect_equal(cap$args$mat, as.matrix(dist(t(m), method = "canberra")))
})

test_that("plot_dist handles NAs like stats::dist", {
  skip_if_not_installed("pheatmap")
  skip_if_not_installed("viridisLite")
  cap <- capture_pheatmap()
  xn <- x3
  xn[1, 2] <- NA
  plot_dist(xn)
  # s1 vs s2 uses only g2: |0 - 4| scaled by sqrt(2 / 1)
  expect_equal(cap$args$mat["s1", "s2"], 4 * sqrt(2))
  expect_equal(cap$args$mat, as.matrix(dist(t(xn))))
})

test_that("plot_dist returns a symmetric matrix in the original sample order", {
  skip_if_not_installed("pheatmap")
  skip_if_not_installed("viridisLite")
  cap <- capture_pheatmap()
  set.seed(7)
  m <- matrix(rnorm(60), 10, dimnames = list(NULL, c("z", "a", "m", "b", "q", "c")))
  plot_dist(as.data.frame(m))
  d <- cap$args$mat
  expect_equal(d, t(d))
  expect_equal(unname(diag(d)), rep(0, 6))
  expect_true(all(d >= 0))
  expect_equal(rownames(d), colnames(m))
  expect_equal(colnames(d), colnames(m))
})

test_that("plot_dist applies defaults and lets ... override them", {
  skip_if_not_installed("pheatmap")
  skip_if_not_installed("viridisLite")
  cap <- capture_pheatmap()
  plot_dist(x3)
  expect_false(cap$args$cluster_rows)
  expect_false(cap$args$cluster_cols)
  expect_true(cap$args$display_numbers)
  expect_equal(cap$args$scale, "none")

  plot_dist(x3, cluster_rows = TRUE, number_format = "%.3f", main = "M")
  expect_true(cap$args$cluster_rows)
  expect_equal(cap$args$number_format, "%.3f")
  expect_equal(cap$args$main, "M")
  expect_equal(sum(names(cap$args) == "cluster_rows"), 1)
})

test_that("plot_dist attaches metadata as row and column annotations", {
  skip_if_not_installed("pheatmap")
  skip_if_not_installed("viridisLite")
  cap <- capture_pheatmap()
  md <- data.frame(row.names = c("s1", "s2", "s3"), Group = c("a", "a", "b"))
  plot_dist(x3, metadata = md)
  expect_identical(cap$args$annotation_row, md)
  expect_identical(cap$args$annotation_col, md)

  plot_dist(x3)
  expect_null(cap$args$annotation_row)
  expect_error(plot_dist(x3, metadata = md[c(2, 1, 3), , drop = FALSE]), "rownames")
})

test_that("plot_dist clustering tree matches hclust on stats::dist", {
  skip_if_not_installed("pheatmap")
  skip_if_not_installed("viridisLite")
  pdf(NULL)
  on.exit(dev.off())
  h <- plot_dist(GSE161650_lc, cluster_cols = TRUE, silent = TRUE)
  expect_s3_class(h, "pheatmap")
  # pheatmap clusters the columns of the distance matrix by euclidean distance
  d <- as.matrix(dist(t(GSE161650_lc)))
  expected <- hclust(dist(t(d)), method = "complete")
  expect_equal(h$tree_col$merge, expected$merge)
  expect_equal(h$tree_col$height, expected$height)
  expect_equal(h$tree_col$labels, colnames(GSE161650_lc))
})
