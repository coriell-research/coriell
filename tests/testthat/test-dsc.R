dsc_oracle <- function(x, g) {
  g <- factor(g)
  N <- ncol(x)
  mu <- rowMeans(x)
  Sb <- 0
  Sw <- 0
  for (l in levels(g)) {
    xj <- x[, g == l, drop = FALSE]
    p <- ncol(xj) / N
    mj <- rowMeans(xj)
    Sb <- Sb + p * tcrossprod(mj - mu)
    Sw <- Sw + p * tcrossprod(xj - mj) / ncol(xj)
  }
  sqrt(sum(diag(Sb))) / sqrt(sum(diag(Sw)))
}

test_that("dsc hand-computed values", {
  # global mean 6, batch means 1 and 11 -> tr(Sb) = 25, tr(Sw) = 4/4 = 1
  expect_equal(dsc(matrix(c(0, 2, 10, 12), nrow = 1), c("a", "a", "b", "b")), 5)
  # feature 1: batch means 2 and 10, global 4 -> Sb = 12, Sw = 8/4 = 2
  x2 <- rbind(c(0, 2, 4, 10), c(1, 1, 1, 1))
  expect_equal(dsc(x2, c("a", "a", "a", "b")), sqrt(12) / sqrt(2))
})

test_that("dsc matches the scatter-matrix definition and its invariances", {
  set.seed(1)
  x <- matrix(rnorm(20 * 12), 20)
  g <- rep(c("A", "B", "C"), times = c(3, 4, 5))
  x[, g == "C"] <- x[, g == "C"] + 1
  d <- dsc(x, g)
  expect_equal(d, dsc_oracle(x, g))
  expect_equal(dsc(x * 7.3, g), d)
  expect_equal(dsc(x + rnorm(20), g), d)
  p <- sample(12)
  expect_equal(dsc(x[, p], g[p]), d)
  expect_equal(dsc(x, c(A = "z", B = "y", C = "x")[g]), d)
  expect_equal(dsc(as.data.frame(x), factor(g)), d)
  expect_equal(dsc(x, factor(g, levels = c("A", "B", "C", "D"))), d)
  # total sum of squares decomposes as N * (Db^2 + Dw^2)
  Db2 <- 0
  for (l in unique(g)) {
    Db2 <- Db2 + mean(g == l) * sum((rowMeans(x[, g == l]) - rowMeans(x))^2)
  }
  expect_equal(ncol(x) * (Db2 + Db2 / d^2), sum((x - rowMeans(x))^2))
})

test_that("dsc limiting cases", {
  x <- cbind(c(1, 5), c(3, 7), c(1, 5), c(3, 7))
  expect_equal(dsc(x, c(1, 1, 2, 2)), 0)
  expect_equal(dsc(x, c(1, 2, 1, 2)), Inf)
  expect_equal(dsc(x, rep(1, 4)), 0)
  set.seed(9)
  y <- matrix(rnorm(500 * 100), 500)
  b <- rep(1:2, each = 50)
  y2 <- y
  y2[, 51:100] <- y2[, 51:100] + 50
  expect_lt(dsc(y, b), 0.2)
  expect_gt(dsc(y2, b), 10)
})

test_that("dsc drops samples with an NA batch label, with a warning", {
  x <- matrix(c(0, 2, 10, 12, 100), nrow = 1)
  expect_warning(d <- dsc(x, c("a", "a", "b", "b", NA)), "NA batch label")
  expect_equal(d, 5)
  expect_warning(d2 <- dsc(x, factor(c("a", "a", "b", "b", NA))), "NA batch label")
  expect_equal(d2, 5)
})
