mk_dge <- function() {
  set.seed(42)
  counts <- matrix(
    rnbinom(2000, mu = 50, size = 10),
    ncol = 8,
    dimnames = list(paste0("g", 1:250), paste0("s", 1:8))
  )
  edgeR::DGEList(counts, group = gl(2, 4))
}

test_that("method = 'none' gives unit factors and keeps counts", {
  skip_if_not_installed("edgeR")
  y <- mk_dge()
  r <- normalize_counts(y, method = "none")
  expect_equal(unname(r$samples$norm.factors), rep(1, 8))
  expect_identical(r$counts, y$counts)
  expect_null(r$offset)
})

test_that("upperquartile matches hand computation", {
  skip_if_not_installed("edgeR")
  y <- mk_dge()
  x <- y$counts
  ls <- colSums(x)
  for (p in c(0.75, 0.5)) {
    f <- apply(sweep(x, 2, ls, "/"), 2, quantile, probs = p)
    f <- f / exp(mean(log(f)))
    r <- normalize_counts(y, method = "upperquartile", p = p)
    expect_equal(unname(r$samples$norm.factors), unname(f))
  }
})

test_that("RLE matches hand computation", {
  skip_if_not_installed("edgeR")
  y <- mk_dge()
  x <- y$counts
  ls <- colSums(x)
  gm <- exp(rowMeans(log(x)))
  keep <- gm > 0
  f <- apply(x[keep, ], 2, function(u) median(u / gm[keep])) / ls
  f <- f / exp(mean(log(f)))
  r <- normalize_counts(y, method = "RLE")
  expect_equal(unname(r$samples$norm.factors), unname(f))
})

test_that("standard methods agree with edgeR::normLibSizes including passed args", {
  skip_if_not_installed("edgeR")
  y <- mk_dge()
  for (m in c("TMM", "TMMwsp", "RLE", "upperquartile")) {
    expect_equal(
      normalize_counts(y, method = m)$samples$norm.factors,
      edgeR::normLibSizes(y, method = m)$samples$norm.factors
    )
  }
  expect_equal(
    normalize_counts(
      y,
      method = "TMM",
      refColumn = 3,
      logratioTrim = 0.2,
      sumTrim = 0.1,
      doWeighting = FALSE
    )$samples$norm.factors,
    edgeR::normLibSizes(
      y,
      method = "TMM",
      refColumn = 3,
      logratioTrim = 0.2,
      sumTrim = 0.1,
      doWeighting = FALSE
    )$samples$norm.factors
  )
  expect_false(isTRUE(all.equal(
    normalize_counts(
      y,
      method = "TMM",
      doWeighting = FALSE
    )$samples$norm.factors,
    normalize_counts(y, method = "TMM")$samples$norm.factors
  )))
})

test_that("TMM: geometric mean 1, scale invariance, composition bias removed", {
  skip_if_not_installed("edgeR")
  y <- mk_dge()
  expect_equal(prod(normalize_counts(y)$samples$norm.factors), 1)

  # proportions are unchanged by scaling a column (unweighted TMM)
  r <- normalize_counts(y, doWeighting = FALSE)
  y2 <- y
  y2$counts[, 2] <- y2$counts[, 2] * 3L
  y2$samples$lib.size <- colSums(y2$counts)
  expect_equal(
    normalize_counts(y2, doWeighting = FALSE)$samples$norm.factors,
    r$samples$norm.factors
  )

  set.seed(1)
  base <- rnbinom(500, mu = 100, size = 5)
  yi <- edgeR::DGEList(matrix(rep(base, 4), ncol = 4))
  expect_equal(unname(normalize_counts(yi)$samples$norm.factors), rep(1, 4))

  # 10x up-regulation of 10% of genes in samples 3-4 is trimmed: effective
  # library sizes become equal
  m <- matrix(rep(base, 4), ncol = 4)
  m[1:50, 3:4] <- m[1:50, 3:4] * 10L
  nf <- normalize_counts(edgeR::DGEList(m))$samples$norm.factors
  eff <- colSums(m) * nf
  expect_equal(nf[1], nf[2])
  expect_equal(nf[3], nf[4])
  expect_gt(nf[1], nf[3])
  expect_equal(unname(eff[3] / eff[1]), 1, tolerance = 1e-6)
})

test_that("qsmooth offsets equal the centred log-ratio to qsmoothData", {
  skip_if_not_installed("edgeR")
  skip_if_not_installed("qsmooth")
  y <- mk_dge()
  g <- y$samples$group
  r <- normalize_counts(y, method = "qsmooth", group_factor = g)
  expect_s4_class(r, "DGEList")
  expect_identical(r$counts, y$counts)
  expect_equal(dim(r$offset), dim(y$counts))
  qsd <- qsmooth::qsmoothData(qsmooth::qsmooth(y$counts, group_factor = g))
  L <- log(y$counts + 0.5) - log(qsd + 0.5)
  expected <- L - rowMeans(L) + mean(log(y$samples$lib.size))
  expect_equal(unname(as.matrix(r$offset)), unname(expected))

  r2 <- normalize_counts(
    y,
    method = "qsmooth",
    group_factor = g,
    prior_count = 1,
    window = 0.1
  )
  qsd2 <- qsmooth::qsmoothData(qsmooth::qsmooth(
    y$counts,
    group_factor = g,
    window = 0.1
  ))
  L2 <- log(y$counts + 1) - log(qsd2 + 1)
  expect_equal(
    unname(as.matrix(r2$offset)),
    unname(L2 - rowMeans(L2) + mean(log(y$samples$lib.size)))
  )
})

test_that("qsmooth on identical libraries gives constant offsets", {
  skip_if_not_installed("edgeR")
  skip_if_not_installed("qsmooth")
  set.seed(1)
  base <- rnbinom(300, mu = 100, size = 5)
  y <- edgeR::DGEList(matrix(rep(base, 4), ncol = 4))
  r <- normalize_counts(y, method = "qsmooth", group_factor = gl(2, 2))
  expect_equal(unname(as.matrix(r$offset)), matrix(log(sum(base)), 300, 4))
})

test_that("qsmooth requires group_factor", {
  skip_if_not_installed("edgeR")
  skip_if_not_installed("qsmooth")
  expect_error(normalize_counts(mk_dge(), method = "qsmooth"), "group_factor")
})

test_that("loess offsets match an independent loessFit recomputation", {
  skip_if_not_installed("edgeR")
  skip_if_not_installed("csaw")
  skip_if_not_installed("limma")
  y <- mk_dge()
  r <- normalize_counts(y, method = "loess")
  ls <- y$samples$lib.size
  expect_equal(unname(rowMeans(r$offset)), rep(mean(log(ls)), 250))
  expect_equal(r$samples$norm.factors, y$samples$norm.factors)

  cc <- 0.5 * ls / mean(ls)
  ab <- edgeR::aveLogCPM(y$counts, lib.size = ls, prior.count = 0.5) /
    log2(exp(1))
  fits <- sapply(1:8, function(j) {
    limma::loessFit(
      log(y$counts[, j] + cc[j]),
      ab,
      span = 0.3,
      iterations = 4L
    )$fitted
  })
  expect_equal(unname(r$offset), unname(fits - rowMeans(fits) + mean(log(ls))))

  r2 <- normalize_counts(y, method = "loess", span = 0.6, iterations = 1L)
  expect_equal(
    r2$offset,
    csaw::normOffsets(y, span = 0.6, iterations = 1L)$offset
  )
  expect_false(isTRUE(all.equal(r$offset, r2$offset)))
})

test_that("loess on identical libraries gives constant offsets", {
  skip_if_not_installed("edgeR")
  skip_if_not_installed("csaw")
  set.seed(1)
  base <- rnbinom(300, mu = 100, size = 5)
  y <- edgeR::DGEList(matrix(rep(base, 4), ncol = 4))
  expect_equal(
    unname(normalize_counts(y, method = "loess")$offset),
    matrix(log(sum(base)), 300, 4)
  )
})
