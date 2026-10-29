exp1 <- data.frame(
  feature_id = c("geneA", "geneB", "geneC"),
  PValue = c(0.01, 0.5, 0.05),
  logFC = c(1.2, -2.5, 3.7),
  SE = c(.5, .4, 1)
)
exp2 <- data.frame(
  feature_id = c("geneA", "geneB", "geneD"),
  PValue = c(0.07, 0.3, 0.8),
  logFC = c(1.5, -2, 3),
  SE = c(.6, .45, 1.2)
)

mk_se <- function() {
  p <- matrix(
    c(0.01, 0.2, 0.5, 0.02, 0.3, 0.6, 0.04, 0.1, 0.7),
    3,
    dimnames = list(c("g1", "g2", "g3"), c("s1", "s2", "s3"))
  )
  l <- matrix(c(1, -1, 1, 2, 1, -1, 1.5, -2, 0.5), 3, dimnames = dimnames(p))
  s <- matrix(c(.5, .5, 1, .25, 1, 1, 1, .5, 2), 3, dimnames = dimnames(p))
  SummarizedExperiment::SummarizedExperiment(
    assays = list(PValue = p, logFC = l, SE = s)
  )
}

test_that("dfs2se places each value at its feature/experiment cell", {
  skip_if_not_installed("SummarizedExperiment")
  se <- dfs2se(list(e1 = exp1, e2 = exp2), import = c("PValue", "logFC", "SE"))
  expect_s4_class(se, "SummarizedExperiment")
  expect_identical(SummarizedExperiment::assayNames(se), c("PValue", "logFC", "SE"))
  expect_identical(rownames(se), c("geneA", "geneB", "geneC", "geneD"))
  expect_identical(colnames(se), c("e1", "e2"))
  for (a in c("PValue", "logFC", "SE")) {
    m <- SummarizedExperiment::assay(se, a)
    expect_equal(unname(m[exp1$feature_id, "e1"]), exp1[[a]])
    expect_equal(unname(m[exp2$feature_id, "e2"]), exp2[[a]])
    expect_true(is.na(m["geneD", "e1"]))
    expect_true(is.na(m["geneC", "e2"]))
  }
})

test_that("dfs2se is keyed by feature id, not row position", {
  skip_if_not_installed("SummarizedExperiment")
  a <- dfs2se(list(e1 = exp1, e2 = exp2), import = "logFC")
  b <- dfs2se(list(e1 = exp1[3:1, ], e2 = exp2[c(2, 3, 1), ]), import = "logFC")
  expect_identical(SummarizedExperiment::assay(a), SummarizedExperiment::assay(b))
})

test_that("dfs2se complete = TRUE keeps only shared features", {
  skip_if_not_installed("SummarizedExperiment")
  se <- dfs2se(list(e1 = exp1, e2 = exp2), import = c("PValue", "logFC"), complete = TRUE)
  expect_identical(rownames(se), intersect(exp1$feature_id, exp2$feature_id))
  expect_false(anyNA(SummarizedExperiment::assay(se, "logFC")))
  expect_equal(unname(SummarizedExperiment::assay(se, "logFC")[, "e2"]), c(1.5, -2))
})

test_that("dfs2se honours a custom feature_col", {
  skip_if_not_installed("SummarizedExperiment")
  g <- exp1
  names(g)[1] <- "gene"
  se <- dfs2se(list(only = g), feature_col = "gene", import = "logFC")
  expect_equal(dim(se), c(3L, 1L))
  expect_equal(SummarizedExperiment::assay(se)[, 1], setNames(exp1$logFC, exp1$feature_id))
})

test_that("dfs2se stops on mismatched or missing columns", {
  skip_if_not_installed("SummarizedExperiment")
  expect_error(dfs2se(list(a = exp1, b = exp2[, c(1, 3, 2, 4)])), "column names must match")
  expect_error(dfs2se(list(a = exp1, b = exp2), import = "FDR"), "missing")
})

test_that("dfs2se works when input columns are named 'rows' or 'cols'", {
  skip_if_not_installed("SummarizedExperiment")
  f1 <- data.frame(gene = c("a", "b"), cols = c(1, 2), logFC = c(1, 2))
  f2 <- f1
  f2$logFC <- c(3, 4)
  expected <- matrix(c(1, 2, 3, 4), 2, dimnames = list(c("a", "b"), c("x", "y")))
  se <- dfs2se(list(x = f1, y = f2), feature_col = "gene", import = "logFC")
  expect_equal(SummarizedExperiment::assay(se), expected)

  names(f1)[1] <- names(f2)[1] <- "rows"
  se2 <- dfs2se(list(x = f1, y = f2), feature_col = "rows", import = "logFC")
  expect_equal(SummarizedExperiment::assay(se2), expected)
})

test_that("dfs2se errors on duplicate feature ids within an experiment", {
  skip_if_not_installed("SummarizedExperiment")
  d1 <- data.frame(feature_id = c("a", "a", "b"), logFC = c(5, 7, 2))
  d2 <- data.frame(feature_id = c("a", "b"), logFC = c(1, 2))
  expect_error(dfs2se(list(x = d1, y = d2), import = "logFC"), "Duplicate feature ids.*x")
})

test_that("meta_de combined p-values match closed-form combinations", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("metapod")
  x <- mk_se()
  p <- SummarizedExperiment::assay(x, "PValue")

  fisher <- apply(p, 1, function(v) {
    pchisq(-2 * sum(log(v)), 2 * length(v), lower.tail = FALSE)
  })
  r <- meta_de(x, metapod::parallelFisher)
  expect_equal(r$Combined.Pval, unname(fisher))
  expect_identical(r$Feature, rownames(x))

  stouffer <- apply(p, 1, function(v) {
    pnorm(sum(qnorm(v, lower.tail = FALSE)) / sqrt(length(v)), lower.tail = FALSE)
  })
  expect_equal(meta_de(x, metapod::parallelStouffer)$Combined.Pval, unname(stouffer))

  simes <- apply(p, 1, function(v) min(sort(v) * length(v) / seq_along(v)))
  expect_equal(meta_de(x, metapod::parallelSimes)$Combined.Pval, unname(simes))
})

test_that("meta_de passes ... through to FUN", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("metapod")
  x <- mk_se()
  # Holm-min with min.prop = 1/3 is the smallest Holm-adjusted p = 3 * min(p)
  r1 <- meta_de(x, metapod::parallelHolmMin, min.prop = 1 / 3)
  r2 <- meta_de(x, metapod::parallelHolmMin, min.prop = 1)
  expect_equal(r1$Combined.Pval, c(0.03, 0.3, 1))
  expect_false(isTRUE(all.equal(r1$Combined.Pval, r2$Combined.Pval)))
})

test_that("meta_de logFC summaries match base R", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("metapod")
  x <- mk_se()
  p <- SummarizedExperiment::assay(x, "PValue")
  l <- SummarizedExperiment::assay(x, "logFC")
  comb <- metapod::parallelFisher(asplit(p, 2))
  r <- meta_de(x, metapod::parallelFisher)
  idx <- cbind(seq_len(nrow(p)), comb$representative)
  expect_equal(r$Rep.Pval, unname(p[idx]))
  expect_equal(r$Rep.logFC, unname(l[idx]))
  expect_equal(r$Median.logFC, unname(apply(l, 1, median)))
  expect_equal(r$Mean.logFC, unname(rowMeans(l)))
  expect_equal(r$Min.logFC, unname(apply(l, 1, min)))
  expect_equal(r$Max.logFC, unname(apply(l, 1, max)))
})

test_that("meta_de direction follows the signs of influential tests", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("metapod")
  p <- matrix(0.01, 3, 2, dimnames = list(c("u", "d", "m"), c("a", "b")))
  l <- rbind(u = c(1, 2), d = c(-1, -3), m = c(1, -1))
  colnames(l) <- c("a", "b")
  x <- SummarizedExperiment::SummarizedExperiment(assays = list(PValue = p, logFC = l))
  expect_identical(meta_de(x, metapod::parallelFisher)$Direction, c("up", "down", "mixed"))
})

test_that("meta_de fixed-effect estimate is the inverse-variance weighted mean", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("metapod")
  x <- mk_se()
  r <- meta_de(x, metapod::parallelFisher, se = "SE")
  # g1: lfc = (1, 2, 1.5), SE = (.5, .25, 1) -> w = (4, 16, 1)
  b1 <- (4 * 1 + 16 * 2 + 1 * 1.5) / 21
  z1 <- b1 * sqrt(21)
  expect_equal(r$Meta.logFC[1], b1)
  expect_equal(r$Meta.z[1], z1)
  expect_equal(r$Meta.Pval[1], 2 * pnorm(-abs(z1)))

  # equal SEs reduce to the arithmetic mean
  SummarizedExperiment::assay(x, "SE")[] <- 0.7
  r_eq <- meta_de(x, metapod::parallelFisher, se = "SE")
  expect_equal(r_eq$Meta.logFC, unname(rowMeans(SummarizedExperiment::assay(x, "logFC"))))

  # a single study reduces to a Wald z
  r1 <- meta_de(mk_se()[, 1], metapod::parallelFisher, se = "SE")
  expect_equal(r1$Combined.Pval, c(0.01, 0.2, 0.5))
  expect_equal(r1$Meta.z, c(1, -1, 1) / c(.5, .5, 1))

  # no SE -> meta columns are NA
  r0 <- meta_de(mk_se(), metapod::parallelFisher)
  expect_true(all(is.na(r0$Meta.logFC) & is.na(r0$Meta.Pval) & is.na(r0$Meta.z)))
})

test_that("meta_de handles missing values with and without imputation", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("metapod")
  p <- matrix(c(0.01, NA, 0.02, 0.3, 0.04, 0.1), 2, dimnames = list(c("g1", "g2"), c("s1", "s2", "s3")))
  l <- matrix(c(1, NA, 2, 1, 1.5, -2), 2, dimnames = dimnames(p))
  s <- matrix(c(.5, NA, .25, 1, 1, .5), 2, dimnames = dimnames(p))
  x <- SummarizedExperiment::SummarizedExperiment(assays = list(PValue = p, logFC = l, SE = s))
  ri <- meta_de(x, metapod::parallelFisher, se = "SE")
  rn <- meta_de(x, metapod::parallelFisher, se = "SE", impute_missing = FALSE)
  expect_equal(ri$Combined.Pval[2], pchisq(-2 * sum(log(c(1, .3, .1))), 6, lower.tail = FALSE))
  expect_equal(rn$Combined.Pval[2], pchisq(-2 * sum(log(c(.3, .1))), 4, lower.tail = FALSE))
  w <- 1 / c(1, .5)^2
  expect_equal(ri$Meta.logFC[2], sum(w * c(1, -2)) / sum(w))
  expect_equal(rn$Meta.logFC[2], sum(w * c(1, -2)) / sum(w))
  expect_equal(rn$Mean.logFC[2], -0.5)
  expect_equal(rn$Median.logFC[2], -0.5)
})

test_that("meta_de drops studies with a missing logFC but present SE from the fixed-effect estimate", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("metapod")
  p <- matrix(c(.01, .02, .04, .5, .5, .5), 2, byrow = TRUE, dimnames = list(c("g1", "g2"), c("s1", "s2", "s3")))
  l <- matrix(c(NA, 2, 1.5, NA, NA, NA), 2, byrow = TRUE, dimnames = dimnames(p))
  s <- matrix(c(.5, .25, 1, 1, 1, 1), 2, byrow = TRUE, dimnames = dimnames(p))
  x <- SummarizedExperiment::SummarizedExperiment(assays = list(PValue = p, logFC = l, SE = s))

  w <- 1 / c(.25, 1)^2
  b <- sum(w * c(2, 1.5)) / sum(w)
  z <- b * sqrt(sum(w))
  for (imp in c(TRUE, FALSE)) {
    r <- meta_de(x, metapod::parallelFisher, se = "SE", impute_missing = imp)
    expect_equal(r$Meta.logFC[1], b)
    expect_equal(r$Meta.z[1], z)
    expect_equal(r$Meta.Pval[1], 2 * pnorm(-abs(z)))
    # no observed logFC at all -> NA, not NaN/0
    expect_true(is.na(r$Meta.logFC[2]) && !is.nan(r$Meta.logFC[2]))
    expect_true(is.na(r$Meta.z[2]) && !is.nan(r$Meta.z[2]))
    expect_true(is.na(r$Meta.Pval[2]) && !is.nan(r$Meta.Pval[2]))
  }

  skip_if_not_installed("DelayedArray")
  skip_if_not_installed("DelayedMatrixStats")
  for (a in c("PValue", "logFC", "SE")) {
    SummarizedExperiment::assay(x, a, withDimnames = FALSE) <-
      DelayedArray::DelayedArray(SummarizedExperiment::assay(x, a, withDimnames = FALSE))
  }
  rd <- meta_de(x, metapod::parallelFisher, se = "SE")
  expect_equal(rd$Meta.logFC[1], b)
  expect_true(is.na(rd$Meta.logFC[2]) && !is.nan(rd$Meta.logFC[2]))
})

test_that("meta_de gives the same result for DelayedArray assays", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("metapod")
  skip_if_not_installed("DelayedArray")
  skip_if_not_installed("DelayedMatrixStats")
  x <- mk_se()
  r <- meta_de(x, metapod::parallelFisher, se = "SE")
  for (a in c("PValue", "logFC", "SE")) {
    SummarizedExperiment::assay(x, a, withDimnames = FALSE) <-
      DelayedArray::DelayedArray(SummarizedExperiment::assay(x, a, withDimnames = FALSE))
  }
  expect_s4_class(SummarizedExperiment::assay(x, "logFC"), "DelayedArray")
  expect_equal(as.data.frame(meta_de(x, metapod::parallelFisher, se = "SE")), as.data.frame(r))
})

test_that("jackknifeSE drops each column once and passes ... to FUN", {
  skip_if_not_installed("SummarizedExperiment")
  x <- mk_se()
  res <- jackknifeSE(x, function(y, k) list(cols = colnames(y), k = k, n = nrow(y)), k = 42)
  expect_length(res, ncol(x))
  for (i in seq_len(ncol(x))) {
    expect_identical(res[[i]]$cols, colnames(x)[-i])
    expect_identical(res[[i]]$k, 42)
    expect_identical(res[[i]]$n, nrow(x))
  }
})

test_that("jackknifeSE of the column mean gives the classic jackknife SE", {
  skip_if_not_installed("SummarizedExperiment")
  set.seed(1)
  v <- rnorm(7)
  x <- SummarizedExperiment::SummarizedExperiment(assays = list(a = matrix(v, 1)))
  theta <- unlist(jackknifeSE(x, function(y) mean(SummarizedExperiment::assay(y))))
  n <- length(v)
  jk_se <- sqrt((n - 1) / n * sum((theta - mean(theta))^2))
  expect_equal(jk_se, sd(v) / sqrt(n))
})
