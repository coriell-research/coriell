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

test_that("permutations matches an independent lexicographic enumeration", {
  for (n in 2:5) {
    g <- as.matrix(expand.grid(rep(list(seq_len(n)), n)))
    g <- g[apply(g, 1, function(r) length(unique(r)) == n), , drop = FALSE]
    g <- g[do.call(order, as.data.frame(g)), , drop = FALSE]
    dimnames(g) <- NULL
    p <- permutations(seq_len(n))
    expect_equal(nrow(p), factorial(n))
    expect_identical(p, g)
    # each value appears (n-1)! times in every column
    expect_true(all(apply(p, 2, function(col) all(table(col) == factorial(n - 1)))))
  }
})

test_that("permutations permutes positions, not values", {
  x <- c(10, 10, 20)
  p <- permutations(x)
  expect_equal(dim(p), c(6, 3))
  expect_equal(nrow(unique(p)), 3)
  expect_type(permutations(c(TRUE, FALSE)), "logical")
  expect_identical(permutations(c("b", "a")), matrix(c("b", "a", "a", "b"), 2))
})

test_that("Exact Spearman with ties in X and y matches brute force", {
  set.seed(7)
  X <- matrix(sample(1:3, 30 * 6, replace = TRUE), 30)
  X <- X[apply(X, 1, var) > 0, ]
  y <- c(1, 1, 2, 3, 3, 4)
  perms <- unique(permutations(y))
  expect_equal(nrow(perms), factorial(6) / (2 * 2))
  obs <- cor(t(X), y, method = "spearman")[, 1]
  null <- cor(t(X), t(perms), method = "spearman")
  tol <- 1e-8
  res <- exact_cor_test(X, y, method = "spearman")
  expect_equal(res$cor, obs)
  expect_equal(res$empirical.p, rowMeans(abs(null) >= abs(obs) - tol))
  g <- exact_cor_test(X, y, method = "spearman", alternative = "greater")
  expect_equal(g$empirical.p, rowMeans(null >= obs - tol))
  l <- exact_cor_test(X, y, method = "spearman", alternative = "less")
  expect_equal(l$empirical.p, rowMeans(null <= obs + tol))
  expect_equal(res$FDR, p.adjust(res$empirical.p, "BH"))
})

test_that("One-sided exact p-values are complementary and lie on the 1/N lattice", {
  set.seed(3)
  X <- matrix(rnorm(25 * 6), 25)
  y <- rnorm(6)
  g <- exact_cor_test(X, y, alternative = "greater")$empirical.p
  l <- exact_cor_test(X, y, alternative = "less")$empirical.p
  t2 <- exact_cor_test(X, y)$empirical.p
  # observed ordering is counted in both tails
  expect_equal(g + l, rep(1 + 1 / factorial(6), 25))
  expect_true(all(t2 >= pmin(g, l) - 1e-12))
  expect_equal(t2 * factorial(6), round(t2 * factorial(6)))
})

test_that("Exact test is invariant to affine transforms; negating y flips tails", {
  set.seed(4)
  X <- matrix(rnorm(15 * 7), 15)
  y <- rnorm(7)
  a <- exact_cor_test(X, y)
  b <- exact_cor_test(3 * X + 5, 2 * y - 1)
  expect_equal(a$cor, b$cor)
  expect_equal(a$empirical.p, b$empirical.p)
  c <- exact_cor_test(X, -y)
  expect_equal(c$cor, -a$cor)
  expect_equal(c$empirical.p, a$empirical.p)
  expect_equal(
    exact_cor_test(X, y, alternative = "greater")$empirical.p,
    exact_cor_test(X, -y, alternative = "less")$empirical.p
  )
})

test_that("data.frame input keeps dimnames and values", {
  df <- data.frame(
    s1 = c(1, 4, 2), s2 = c(2, 3, 2), s3 = c(3, 2, 5), s4 = c(4, 1, 1),
    row.names = c("g1", "g2", "g3")
  )
  res <- exact_cor_test(df, 1:4)
  expect_identical(rownames(res), c("g1", "g2", "g3"))
  expect_identical(colnames(res), c("s1", "s2", "s3", "s4", "cor", "empirical.p", "FDR"))
  expect_equal(as.matrix(res[, 1:4]), as.matrix(df))
  expect_equal(res$cor, c(1, -1, cor(c(2, 2, 5, 1), 1:4)))
  expect_equal(res["g1", "empirical.p"], 2 / 24)
})

test_that("Near-tied y values do not inflate exact p-values", {
  # 0.1 + 0.2 != 0.3 exactly, so all 24 orderings are distinct
  y <- c(0.1 + 0.2, 0.3, 1, 2)
  X <- rbind(c(1, -1, -1, 1), c(1, 2, 3, 4))
  res <- exact_cor_test(X, y)
  P <- unique(permutations(y))
  expect_equal(nrow(P), 24)
  obs <- cor(t(X), y)[, 1]
  expect_true(all(res$empirical.p <= 1))
  expect_equal(res$empirical.p, rowMeans(abs(cor(t(X), t(P))) >= abs(obs) - 1e-8))
})

test_that("Sampled p-values equal statmod::permp on brute-force counts", {
  skip_if_not_installed("statmod")
  set.seed(1)
  X <- matrix(rnorm(40 * 9), 40)
  y <- c(rnorm(8), 0)
  n_perm <- 300
  set.seed(11)
  res <- permutation_correlation_test(X, y, n_perm = n_perm, chunk_size = 64)
  set.seed(11)
  idx <- vapply(seq_len(n_perm), function(i) sample.int(9), integer(9))
  null <- cor(t(X), matrix(y[idx], 9))
  obs <- cor(t(X), y)[, 1]
  b <- rowSums(abs(null) >= abs(obs) - sqrt(.Machine$double.eps))
  expect_equal(res$cor, obs)
  expect_equal(
    res$empirical.p,
    statmod::permp(b, nperm = n_perm, total.nperm = factorial(9))
  )
  expect_true(all(res$empirical.p > 0 & res$empirical.p <= 1))
})

test_that("Sampled Spearman one-sided test with tied y matches brute force", {
  skip_if_not_installed("statmod")
  set.seed(1)
  X <- matrix(rnorm(20 * 8), 20)
  y <- c(1, 1, 2, 3, 4, 5, 6, 7)
  set.seed(5)
  res <- permutation_correlation_test(
    X, y, n_perm = 200, method = "spearman", alternative = "greater"
  )
  set.seed(5)
  idx <- vapply(seq_len(200), function(i) sample.int(8), integer(8))
  null <- cor(t(X), matrix(y[idx], 8), method = "spearman")
  obs <- cor(t(X), y, method = "spearman")[, 1]
  b <- rowSums(null >= obs - sqrt(.Machine$double.eps))
  expect_equal(
    res$empirical.p,
    statmod::permp(b, nperm = 200, total.nperm = factorial(8) / 2)
  )
})

test_that("n_core > 1 gives identical results", {
  skip_on_os("windows")
  set.seed(1)
  X <- matrix(rnorm(30 * 8), 30)
  y <- rnorm(8)
  set.seed(9)
  a <- permutation_correlation_test(X, y, n_perm = 400, chunk_size = 50)
  set.seed(9)
  b <- permutation_correlation_test(X, y, n_perm = 400, chunk_size = 50, n_core = 2)
  expect_equal(a, b)
})

test_that("NA rows propagate NA only to themselves", {
  X <- rbind(c(1, NA, 3, 4, 5, 6, 7, 8), c(1, 3, 2, 5, 4, 7, 6, 8))
  set.seed(1)
  res <- permutation_correlation_test(X, 1:8, n_perm = 100)
  expect_true(all(is.na(unlist(res[1, c("cor", "empirical.p", "FDR")]))))
  expect_equal(res$cor[2], cor(X[2, ], 1:8))
  expect_false(is.na(res$empirical.p[2]))
})

test_that("sampled_pvalue keeps NAs and agrees with permp", {
  skip_if_not_installed("statmod")
  cnt <- c(0, 5, NA, 5, 100)
  p <- sampled_pvalue(cnt, 100, 1e6)
  expect_true(is.na(p[3]))
  expect_equal(p[-3], statmod::permp(cnt[-3], nperm = 100, total.nperm = 1e6))
  expect_equal(p[2], p[4])
  expect_true(all(diff(p[c(1, 2, 5)]) > 0))
})

test_that("count_extreme gives hand-computed tail counts", {
  R <- rbind(c(-0.9, 0.2, 0.5, 0.9), c(0.1, 0.1, -0.1, 0.3))
  ts <- c(0.5, 0.1)
  tol <- 1e-8
  expect_equal(count_extreme(R, ts, "two.sided", tol), c(3, 4))
  expect_equal(count_extreme(R, ts, "greater", tol), c(2, 3))
  expect_equal(count_extreme(R, ts, "less", tol), c(3, 3))
})

test_that("standardize_rows gives correlations via tcrossprod", {
  set.seed(2)
  X <- matrix(rnorm(5 * 6), 5)
  for (m in c("pearson", "spearman")) {
    S <- standardize_rows(X, m)
    expect_equal(rowSums(S), rep(0, 5))
    expect_equal(rowSums(S^2), rep(1, 5))
    expect_equal(tcrossprod(S), cor(t(X), method = m))
  }
})

test_that("sample_n_random_cor draws permuted correlations of random rows", {
  set.seed(1)
  X <- matrix(rnorm(50 * 5), 50)
  y <- 1:5
  set.seed(2)
  r <- sample_n_random_cor(X, y, n = 200)
  set.seed(2)
  rows <- sample.int(50, 200, replace = TRUE)
  expected <- vapply(rows, function(i) cor(X[i, ], sample(y)), numeric(1))
  expect_equal(r, expected)
  # Spearman correlations for n = 5 without ties lie on a 1/20 lattice
  set.seed(2)
  rs <- sample_n_random_cor(X, y, n = 200, method = "spearman")
  expect_true(all(abs(rs * 20 - round(rs * 20)) < 1e-8))
  Xc <- rbind(rep(1, 5), rep(2, 5))
  set.seed(3)
  expect_true(all(is.na(suppressWarnings(sample_n_random_cor(Xc, y, n = 10)))))
})

test_that("sample_n_random_cor with n = 1 returns one value", {
  set.seed(1)
  X <- matrix(rnorm(30), 10)
  r <- sample_n_random_cor(X, 1:3, n = 1)
  expect_length(r, 1)
  expect_true(abs(r) <= 1)
})
