test_that("permutations returns all permutations in lexicographic order", {
  p <- permutations(letters[1:4])
  expect_equal(dim(p), c(24, 4))
  expect_equal(anyDuplicated(p), 0)
  expect_equal(p[1, ], letters[1:4])
  expect_equal(p[24, ], letters[4:1])
  expect_equal(permutations("a"), "a")
})

test_that("Exact test counts the observed ordering", {
  X <- rbind(a = 1:6, b = c(3, 1, 4, 1, 5, 9))
  y <- 1:6
  res <- exact_cor_test(X, y)
  # observed and reversed orderings both give |r| == 1
  expect_equal(res["a", "empirical.p"], 2 / factorial(6))
  expect_true(all(res$empirical.p > 0))

  greater <- exact_cor_test(X, y, alternative = "greater")
  expect_equal(greater["a", "empirical.p"], 1 / factorial(6))
  less <- exact_cor_test(X, y, alternative = "less")
  expect_equal(less["a", "empirical.p"], 1)
})

test_that("Exact p-values match brute force with cor()", {
  set.seed(1)
  X <- matrix(rnorm(20 * 6), 20)
  y <- c(1, 2, 2, 3, 5, 8)
  for (method in c("pearson", "spearman")) {
    perms <- unique(permutations(y))
    null <- cor(t(X), t(perms), method = method)
    obs <- cor(t(X), y, method = method)[, 1]
    expected <- rowMeans(abs(null) >= abs(obs) - 1e-8)
    res <- exact_cor_test(X, y, method = method)
    expect_equal(res$cor, obs)
    expect_equal(res$empirical.p, expected)
  }
})

test_that("Zero correlation gives p-value of 1", {
  X <- rbind(c(1, -1, -1, 1))
  y <- 1:4
  res <- exact_cor_test(X, y)
  expect_equal(res$cor, 0)
  expect_equal(res$empirical.p, 1)
})

test_that("Single row and zero variance rows are handled", {
  set.seed(1)
  X <- matrix(rnorm(10), 1)
  expect_equal(nrow(exact_cor_test(X[, 1:5, drop = FALSE], 1:5)), 1)
  expect_equal(nrow(permutation_correlation_test(X, 1:10, n_perm = 100)), 1)

  X <- rbind(rnorm(10), rep(1, 10))
  res <- permutation_correlation_test(X, 1:10, n_perm = 100)
  expect_true(is.na(res$cor[2]))
  expect_true(is.na(res$empirical.p[2]))
  expect_false(is.na(res$empirical.p[1]))
})

test_that("Exact test is used when n_perm exceeds distinct permutations", {
  X <- rbind(1:5)
  expect_message(
    res <- permutation_correlation_test(X, 1:5, n_perm = 1e3),
    "exact test"
  )
  expect_equal(res$empirical.p, 2 / factorial(5))
})

test_that("Sampled p-values are never zero and are calibrated", {
  set.seed(1)
  X <- matrix(rnorm(2000 * 10), 2000)
  y <- rnorm(10)
  res <- permutation_correlation_test(X, y, n_perm = 1000, chunk_size = 300)
  expect_true(all(res$empirical.p > 0))
  expect_equal(res$cor, cor(t(X), y)[, 1])
  expect_gt(mean(res$empirical.p <= 0.05), 0.03)
  expect_lt(mean(res$empirical.p <= 0.05), 0.07)
})

test_that("Sampled results are reproducible across chunk sizes", {
  X <- matrix(rnorm(50 * 9), 50)
  y <- rnorm(9)
  set.seed(42)
  a <- permutation_correlation_test(X, y, n_perm = 500, chunk_size = 1000)
  set.seed(42)
  b <- permutation_correlation_test(X, y, n_perm = 500, chunk_size = 7)
  expect_equal(a, b)
})
