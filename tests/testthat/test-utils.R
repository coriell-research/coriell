test_that("minmax normalizes to [0,1]", {
  x <- c(0, 1, 10, 100)
  expect_equal(minmax(x), c(0, 0.01, 0.1, 1))
  expect_equal(minmax(x, min = 0, max = 1000), x / 1000)
  expect_equal(range(minmax(c(-3, 7, 2.5))), c(0, 1))
  expect_equal(minmax(c(1, NA, 3)), c(0, NA, 1))
  # invariant to positive affine transforms
  expect_equal(minmax(c(2, 4, 6) * 5 + 3), minmax(c(2, 4, 6)))
  # constant vector is 0/0
  expect_true(all(is.nan(minmax(c(5, 5, 5)))))
})

test_that("lerp interpolates and inverts minmax", {
  expect_equal(lerp(c(0, 0.25, 0.5, 1), 10, 20), c(10, 12.5, 15, 20))
  expect_equal(lerp(c(0, 1), 5, -5), c(5, -5))
  x <- c(3, 7, 11)
  expect_equal(lerp(minmax(x), min(x), max(x)), x)
  expect_error(lerp(1.5, 0, 1), "between 0-1")
})

test_that("map_value maps between ranges", {
  expect_equal(map_value(c(0, 1, 2, 3), 0, 10, 0, 1), c(0, 0.1, 0.2, 0.3))
  # Fahrenheit to Celsius
  expect_equal(map_value(c(32, 212, 98.6), 32, 212, 0, 100), c(0, 100, 37))
  # inverted target range
  expect_equal(map_value(c(2, 4), 0, 10, 1, 0), c(0.8, 0.6))
  # round trip
  x <- c(1, 4, 9)
  y <- map_value(x, 0, 10, -1, 1)
  expect_equal(y, c(-0.8, -0.2, 0.8))
  expect_equal(map_value(y, -1, 1, 0, 10), x)
})

test_that("clamp limits values to range", {
  expect_equal(clamp(c(-1, 121, 10, 15), min = 0, max = 100), c(0, 100, 10, 15))
  expect_equal(clamp(c(-5, 0, 0.5, 1, 7), 0, 1), c(0, 0, 0.5, 1, 1))
  expect_equal(
    clamp(c(5, 5, 5), min = c(0, 6, 0), max = c(10, 10, 4)),
    c(5, 6, 4)
  )
  expect_true(is.na(clamp(NA_real_, 0, 1)))
})

test_that("clamp on empty input returns numeric(0)", {
  # regression: Vectorize() returned list() for zero-length input
  expect_identical(clamp(numeric(0), 0, 1), numeric(0))
})

test_that("geometric_mean is correct", {
  expect_equal(geometric_mean(c(2, 8)), 4)
  expect_equal(geometric_mean(c(1, 3, 9)), 3)
  x <- c(2, 5, 95, 5)
  expect_equal(geometric_mean(x), prod(x)^(1 / length(x)))
  expect_equal(geometric_mean(x * 10), 10 * geometric_mean(x))
  # zeros ignored by default
  expect_equal(geometric_mean(c(2, 8, 0, 0)), 4)
  # zeros count in n but not in the log-sum
  expect_equal(geometric_mean(c(2, 8, 0, 0), ignore_zero = FALSE), 16^(1 / 4))
  expect_equal(geometric_mean(c(0, 1, 2), zero_propagate = TRUE), 0)
  expect_identical(geometric_mean(c(NA, 1, 2)), NA)
  expect_true(is.nan(geometric_mean(c(-1, 2, 3))))
  # AM >= GM >= HM
  expect_lte(geometric_mean(x), mean(x))
  expect_gte(geometric_mean(x), length(x) / sum(1 / x))
  expect_equal(geometric_mean(rep(7, 5)), 7)
})

test_that("clr satisfies CLR invariants", {
  expect_equal(clr(c(1, 2, 4, 8)), c(-1.5, -0.5, 0.5, 1.5))
  set.seed(42)
  y <- rexp(20)
  expect_equal(sum(clr(y)), 0)
  expect_equal(clr(y * 13), clr(y))
  expect_equal(clr(y, base = exp(1)), log(y) - mean(log(y)))
  expect_equal(clr(y, base = 10), clr(y, base = 2) * log(2) / log(10))
  # zeros give -Inf which is replaced by 0; gm ignores zeros (gm(2, 8) = 4)
  expect_equal(clr(c(0, 2, 8)), c(0, -1, 1))
})

test_that("colmean averages columns within contiguous groups", {
  m <- matrix(1:12, nrow = 2, dimnames = list(c("r1", "r2"), paste0("s", 1:6)))
  g <- gl(2, 3, labels = c("A", "B"))
  res <- colmean(m, g)
  expect_equal(res, cbind(A = c(r1 = 3, r2 = 4), B = c(r1 = 9, r2 = 10)))
  expect_equal(
    colmean(GSE161650_lc, gl(2, 3, labels = c("DMSO", "THZ1")))[, "THZ1"],
    rowMeans(GSE161650_lc[, 4:6])
  )
})

test_that("colmean is correct when group order differs from level order", {
  # regression: divisor was in level order, sums in first-appearance order
  m <- matrix(1:6, nrow = 2)
  g <- factor(c("b", "b", "a"), levels = c("a", "b"))
  res <- colmean(m, g)
  expect_equal(res[, "a"], c(5, 6))
  expect_equal(res[, "b"], c(2, 3))

  # interleaved groups
  m2 <- matrix(c(1, 10, 2, 20, 3, 30, 4, 40), nrow = 2)
  g2 <- factor(c("y", "x", "y", "y"), levels = c("x", "y"))
  res2 <- colmean(m2, g2)
  expect_equal(res2[, "x"], c(2, 20))
  expect_equal(res2[, "y"], c(8 / 3, 80 / 3))
})

test_that("colmean ignores unused factor levels", {
  # regression: table() included unused levels, misaligning the divisor
  m <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2)
  g <- factor(c("a", "a", "b"), levels = c("a", "b", "c"))
  res <- colmean(m, g)
  expect_equal(sort(colnames(res)), c("a", "b"))
  expect_equal(res[, "a"], c(2, 3))
  expect_equal(res[, "b"], c(5, 6))
  expect_true(all(is.finite(res)))
})

test_that("colmean na.rm gives true means", {
  # regression: NA-removed sums were divided by the full group size
  m <- matrix(c(1, NA, 3, 4, 5, 6), 2)
  g <- factor(c("a", "a", "b"))
  res <- colmean(m, g, na.rm = TRUE)
  expect_equal(unname(res[, "a"]), c(2, 4))
  expect_equal(unname(res[, "b"]), c(5, 6))
  # without na.rm the NA propagates
  expect_true(is.na(colmean(m, g)[2, "a"]))

  # matches a per-group rowMeans(na.rm = TRUE) on random data
  set.seed(3)
  x <- matrix(rnorm(60), 6)
  x[sample(length(x), 12)] <- NA
  grp <- factor(sample(c("p", "q", "r"), 10, replace = TRUE))
  expected <- sapply(levels(grp), function(l) {
    rowMeans(x[, grp == l, drop = FALSE], na.rm = TRUE)
  })
  expect_equal(colmean(x, grp, na.rm = TRUE)[, levels(grp)], expected)
})

test_that("impute fills NAs row-wise", {
  X <- matrix(
    c(1, NA, 3, 10, 20, NA, 4, 4, 4),
    3,
    byrow = TRUE,
    dimnames = list(paste0("r", 1:3), paste0("s", 1:3))
  )
  res <- impute(X)
  expect_equal(
    res,
    rbind(r1 = c(s1 = 1, s2 = 2, s3 = 3), r2 = c(10, 20, 15), r3 = c(4, 4, 4))
  )
  expect_equal(res[!is.na(X)], X[!is.na(X)])
  expect_equal(impute(X, fun = function(v) -1)[is.na(X)], c(-1, -1))
  expect_equal(impute(as.data.frame(X)), res)
  expect_true(all(is.na(impute(rbind(X, r4 = NA))["r4", ])))
  Y <- matrix(1:6, 2)
  expect_equal(impute(Y), Y)
})

test_that("remove_var keeps highest-variance rows", {
  set.seed(1)
  m <- matrix(
    rnorm(40),
    10,
    dimnames = list(paste0("g", 1:10), paste0("s", 1:4))
  ) * (1:10)
  v <- apply(m, 1, var)
  r <- suppressMessages(remove_var(m, p = 0.5))
  expect_equal(rownames(r), names(sort(v, decreasing = TRUE))[1:5])
  expect_equal(r, m[rownames(r), ])
  expect_message(remove_var(m, p = 0.5), "Removing 50% lowest variance")
  expect_equal(suppressMessages(remove_var(as.data.frame(m), 0.5)), r)
  expect_error(remove_var(m, p = 1), "between 0 and 1")
})

test_that("remove_var SummarizedExperiment method agrees with matrix method", {
  skip_if_not_installed("SummarizedExperiment")
  set.seed(1)
  m <- matrix(
    rnorm(40),
    10,
    dimnames = list(paste0("g", 1:10), paste0("s", 1:4))
  ) * (1:10)
  r <- suppressMessages(remove_var(m, p = 0.5))
  se <- SummarizedExperiment::SummarizedExperiment(list(counts = m))
  expect_equal(suppressMessages(remove_var(se, 0.5)), r)
})

test_that("remove_var DelayedArray method agrees with matrix method", {
  skip_if_not_installed("DelayedArray")
  skip_if_not_installed("DelayedMatrixStats")
  set.seed(1)
  m <- matrix(
    rnorm(40),
    10,
    dimnames = list(paste0("g", 1:10), paste0("s", 1:4))
  ) * (1:10)
  r <- suppressMessages(remove_var(m, p = 0.5))
  expect_equal(as.matrix(remove_var(DelayedArray::DelayedArray(m), 0.5)), r)
  # regression: floating point truncation in row count
  expect_equal(nrow(remove_var(DelayedArray::DelayedArray(m), 0.8)), 2)
})

test_that("remove_var ignores NAs in variance", {
  m <- rbind(
    a = c(1, 2, NA, 4),
    b = c(0, 0, 0, 1),
    c = c(10, -10, 10, -10),
    d = c(5, 5, 5, 5)
  )
  expect_equal(rownames(suppressMessages(remove_var(m, 0.5))), c("c", "a"))
})

test_that("remove_var row count is not truncated by floating point", {
  # regression: 10 * (1 - 0.8) = 1.9999... was truncated to 1 by head()
  set.seed(1)
  m <- matrix(rnorm(40), 10)
  expect_equal(nrow(suppressMessages(remove_var(m, p = 0.8))), 2)
  r <- suppressMessages(remove_var(m, p = 0.9))
  expect_true(is.matrix(r))
  expect_equal(dim(r), c(1L, 4L))
})

test_that("remove_var errors informatively on unsupported classes", {
  # regression: default method lacked ..., giving 'unused argument'
  expect_error(remove_var(list(1), p = 0.5), "not supported")
})

test_that("pairwise_intersections computes overlaps and Fisher p-values", {
  l <- list(A = c("a", "b", "c"), B = c("b", "c", "d"), C = c("x", "y"))
  r <- pairwise_intersections(l)
  expect_s3_class(r, "data.table")
  expect_equal(nrow(r), choose(3, 2))
  expect_equal(paste(r$Set1, r$Set2), c("A B", "A C", "B C"))
  expect_equal(r$Set1Size, c(3, 3, 3))
  expect_equal(r$Set2Size, c(3, 2, 2))
  expect_equal(r$IntersectionSize, c(2, 0, 0))
  expect_equal(r$UnionSize, c(4, 5, 5))
  expect_equal(r$Jaccard, c(0.5, 0, 0))
  expect_equal(r$Intersection[[1]], c("b", "c"))
  expect_false("P.Value" %in% names(r))
  expect_equal(
    pairwise_intersections(unname(l))$Set1,
    c("Set1", "Set1", "Set2")
  )

  r2 <- pairwise_intersections(
    list(A = 1:10, B = 6:20, C = 50:60),
    universe_size = 100
  )
  # 5 overlap, |A| = 10, |B| = 15, N = 100
  expect_equal(r2$P.Value[1], phyper(5 - 1, 10, 90, 15, lower.tail = FALSE))
  tab <- table(
    factor(1:100 %in% 1:10, c(FALSE, TRUE)),
    factor(1:100 %in% 6:20, c(FALSE, TRUE))
  )
  expect_equal(r2$P.Value[1], fisher.test(tab, alternative = "greater")$p.value)
  expect_equal(r2$P.Value[2], 1)
  expect_equal(
    nrow(pairwise_intersections(lapply(1:5, function(i) letters[i:(i + 3)]))),
    10
  )
})

test_that("pairwise_intersections handles a single set", {
  # regression: 1:(n - 1) loop gave subscript out of bounds
  r <- pairwise_intersections(list(A = 1:3))
  expect_equal(nrow(r), 0)
  r2 <- pairwise_intersections(list(A = 1:3), universe_size = 10)
  expect_equal(nrow(r2), 0)
})

test_that("pairwise_intersections treats inputs as sets", {
  # regression: duplicates were counted in set sizes but not in union
  l <- list(A = c("a", "a", "b"), B = c("b", "c", "c", "c"))
  r <- pairwise_intersections(l, universe_size = 5)
  expect_equal(r$Set1Size, 2)
  expect_equal(r$Set2Size, 2)
  expect_equal(r$IntersectionSize, 1)
  expect_equal(r$UnionSize, 3)
  expect_equal(r$IntersectionSize + r$UnionSize, r$Set1Size + r$Set2Size)
  # universe {a,b,c,d,e}, A = {a,b}, B = {b,c}
  expect_equal(r$P.Value, phyper(1 - 1, 2, 3, 2, lower.tail = FALSE))
})

test_that("strip_ens removes trailing version only", {
  expect_equal(
    strip_ens(c("ENSG0000001.12", "ENSMUSG00021.3", "ENST00000556.2")),
    c("ENSG0000001", "ENSMUSG00021", "ENST00000556")
  )
  expect_equal(
    strip_ens(c("ENSG0000001", NA, "ENSG2.1")),
    c("ENSG0000001", NA, "ENSG2")
  )
  expect_equal(strip_ens("ENSG.1.2"), "ENSG.1")
  expect_equal(strip_ens("ENSG000.1_PAR_Y"), "ENSG000.1_PAR_Y")
  expect_warning(
    res <- strip_ens(c("ENSG1.1", "TP53.2")),
    "Not all entries start with 'ENS'"
  )
  expect_equal(res, c("ENSG1", "TP53"))
  expect_no_warning(strip_ens(c("ENSG1.1", NA)))
})

test_that("tpm matches the TPM definition", {
  x <- matrix(
    c(10, 20, 30, 5, 0, 15),
    3,
    dimnames = list(paste0("g", 1:3), c("s1", "s2"))
  )
  l <- c(1, 2, 3)
  res <- tpm(x, l)
  expect_equal(
    res,
    matrix(c(1e6 / 3, 1e6 / 3, 1e6 / 3, 5e5, 0, 5e5), 3, dimnames = dimnames(x))
  )
  expect_equal(unname(colSums(res)), c(1e6, 1e6))

  set.seed(7)
  y <- matrix(rpois(60, 50), 20)
  ly <- sample(500:5000, 20)
  # invariant to length units and per-sample depth
  expect_equal(tpm(y, ly * 1000), tpm(y, ly))
  expect_equal(tpm(y %*% diag(c(1, 10, 100)), ly), tpm(y, ly))
  # unit lengths reduce to CPM
  expect_equal(tpm(y, rep(1, 20)), sweep(y, 2, colSums(y), "/") * 1e6)

  expect_equal(tpm(as.data.frame(x), l), res)
  expect_error(tpm(x, 1:2), "number of supplied lengths")
})

test_that("tpm SummarizedExperiment method agrees with matrix method", {
  skip_if_not_installed("SummarizedExperiment")
  x <- matrix(
    c(10, 20, 30, 5, 0, 15),
    3,
    dimnames = list(paste0("g", 1:3), c("s1", "s2"))
  )
  se <- SummarizedExperiment::SummarizedExperiment(list(counts = x))
  expect_equal(tpm(se, c(1, 2, 3)), tpm(x, c(1, 2, 3)))
})

test_that("tpm DelayedArray method agrees with matrix method", {
  skip_if_not_installed("DelayedArray")
  x <- matrix(
    c(10, 20, 30, 5, 0, 15),
    3,
    dimnames = list(paste0("g", 1:3), c("s1", "s2"))
  )
  expect_equal(
    as.matrix(tpm(DelayedArray::DelayedArray(x), c(1, 2, 3))),
    tpm(x, c(1, 2, 3))
  )
})
