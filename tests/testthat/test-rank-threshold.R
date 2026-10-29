rosin_oracle <- function(x) {
  v <- sort(x)
  n <- length(v)
  r <- seq_len(n)
  dx <- n - 1
  dy <- v[n] - v[1]
  d <- abs(dx * (v - v[1]) - dy * (r - 1)) / sqrt(dx^2 + dy^2)
  k <- which.max(d)
  c(x = r[k], y = v[k])
}

test_that("rank_threshold hand examples", {
  # line from (1, 0) to (5, 10): distance of (r, 0) grows with r, max at r = 4
  expect_equal(rank_threshold(c(0, 0, 0, 0, 10)), c(x = 4, y = 0))
  expect_equal(rank_threshold(c(10, 0, 0, 0, 0)), c(x = 4, y = 0))
  expect_equal(rank_threshold(c(1, 2, 3, 4, 20)), c(x = 4, y = 4))
})

test_that("rank_threshold agrees with a geometric oracle and is scale/order invariant", {
  set.seed(42)
  x <- rexp(200)^3
  res <- rank_threshold(x)
  expect_equal(res, rosin_oracle(x))
  expect_equal(rank_threshold(sample(x)), res)
  expect_equal(rank_threshold(x * 1000), c(x = res[["x"]], y = res[["y"]] * 1000))
  xx <- (1:100)^exp(3)
  expect_equal(rank_threshold(xx), rosin_oracle(xx))
})

test_that("rank_threshold show = TRUE gives the same result", {
  pdf(NULL)
  on.exit(dev.off())
  expect_equal(rank_threshold(c(0, 0, 0, 0, 10), show = TRUE), c(x = 4, y = 0))
})

test_that("rank_threshold line ends at the highest rank when the maximum is tied", {
  # line from (1, 1) to (10, 10) is y = x; (6, 10) is furthest from it
  expect_equal(rank_threshold(c(1:5, rep(10, 5))), c(x = 6, y = 10))
  expect_equal(rank_threshold(c(1, 10, 10)), c(x = 2, y = 10))
  set.seed(7)
  z <- c(rexp(50), rep(20, 4))
  expect_equal(rank_threshold(z), rosin_oracle(z))
})
