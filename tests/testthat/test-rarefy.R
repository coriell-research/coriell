test_that("rarefy gives exact depth, is bounded by the input and keeps dimnames", {
  set.seed(123)
  m <- matrix(
    sample.int(100),
    nrow = 10,
    dimnames = list(Gene = paste0("gene.", 1:10), Sample = paste0("sample.", 1:10))
  )
  r <- rarefy(m, depth = 100)
  expect_equal(dim(r), dim(m))
  expect_equal(colSums(r), setNames(rep(100L, 10), colnames(m)))
  expect_true(all(r <= m))
  expect_true(all(r >= 0))
  expect_type(r, "integer")
  expect_identical(rownames(r), rownames(m))
  expect_identical(colnames(r), colnames(m))
})

test_that("rarefy deterministic edge cases", {
  m <- matrix(c(0, 5, 0, 3, 2, 0, 7, 0, 4, 4, 4, 4), nrow = 4, dimnames = list(letters[1:4], c("x", "y", "z")))
  set.seed(1)
  r <- rarefy(m, 8)
  expect_true(all(r[m == 0] == 0))
  # depth equal to the library size returns the column unchanged
  expect_equal(unname(r[, "x"]), c(0, 5, 0, 3))
  expect_equal(unname(rarefy(m, 0)), matrix(0L, 4, 3))
  m2 <- matrix(c(0, 0, 10, 3, 4, 5), nrow = 3)
  set.seed(2)
  expect_equal(rarefy(m2, 6)[, 1], c(0, 0, 6))
})

test_that("rarefy matches hypergeometric moments", {
  m <- matrix(c(60, 30, 10, 0), ncol = 1)
  set.seed(99)
  draws <- replicate(4000, rarefy(m, 20)[, 1])
  expect_equal(rowMeans(draws), c(12, 6, 2, 0), tolerance = 0.03)
  expect_equal(var(draws[1, ]), 20 * 0.6 * 0.4 * (100 - 20) / (100 - 1), tolerance = 0.1)
  expect_true(all(draws[1, ] <= 60))
})

test_that("rarefy is reproducible and accepts a data.frame", {
  m <- matrix(c(10, 20, 30, 40), 2, dimnames = list(c("a", "b"), c("s1", "s2")))
  set.seed(5)
  a <- rarefy(m, 15)
  set.seed(5)
  b <- rarefy(m, 15)
  expect_identical(a, b)
  set.seed(5)
  d <- rarefy(as.data.frame(m), 15)
  expect_identical(d, a)
})

test_that("rarefy guards", {
  m <- matrix(c(1, 2, 3, 40), 2)
  expect_error(rarefy(m, 4), "depth is greater")
  expect_silent(rarefy(m, 3))
  expect_error(rarefy(matrix(c(1, NA, 3, 4), 2), 1), "NAs present")
})

test_that("rarefy keeps the matrix shape for a single-row input", {
  m <- matrix(c(5, 6), nrow = 1, dimnames = list("g", c("s1", "s2")))
  expect_equal(rarefy(m, 3), matrix(3L, 1, 2, dimnames = dimnames(m)))
})

test_that("rarefy assigns the read to the right feature when a pool has one read", {
  m <- matrix(c(0, 0, 1, 0, 0, 5), nrow = 3)
  for (s in 1:20) {
    set.seed(s)
    expect_equal(rarefy(m, 1)[, 1], c(0L, 0L, 1L))
  }
})
