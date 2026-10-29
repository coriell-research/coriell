context("Test edgeR to data.frame function")
library(coriell)
library(edgeR)


# setup results objects ---------------------------------------------------
x <- data.frame(
  sample_A = rnbinom(1000, size = 0.4, prob = 1e-5),
  sample_B = rnbinom(1000, size = 0.4, prob = 1e-5),
  sample_C = rnbinom(1000, size = 0.4, prob = 1e-5),
  sample_D = rnbinom(1000, size = 0.5, prob = 1e-5),
  sample_E = rnbinom(1000, size = 0.5, prob = 1e-5),
  sample_F = rnbinom(1000, size = 0.5, prob = 1e-5),
  sample_G = rnbinom(1000, size = 0.6, prob = 1e-5),
  sample_H = rnbinom(1000, size = 0.6, prob = 1e-5),
  sample_I = rnbinom(1000, size = 0.6, prob = 1e-5),
  row.names = paste0("gene", 1:1000)
)

# edgeR pipeline
group <- factor(rep(c("A", "B", "C"), each = 3))
design <- model.matrix(~ 0 + group)
y <- DGEList(counts = x, group = group)
y <- calcNormFactors(y)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)

# create results objects for all test types
glmqlf_single_contrast <- glmQLFTest(fit, coef = 2)
glmqlf_anova <- glmQLFTest(fit, coef = 1:3)
glmlrt_single_contrast <- glmLRT(fit, coef = 2)
glmlrt_anova <- glmLRT(fit, coef = 1:3)
glmtreat_single_contrast <- glmTreat(
  fit,
  contrast = c(-1, 1, 0),
  lfc = log2(1.2)
)
et_single_contrast <- exactTest(y)


# tests -------------------------------------------------------------------

test_that("Convert results to data.frame", {
  expect_is(
    edger_to_df(glmqlf_single_contrast),
    "data.frame",
    label = "glmQLF test with single coef"
  )
  expect_is(
    edger_to_df(glmqlf_anova),
    "data.frame",
    label = "glmQLF with multiple coefs"
  )
  expect_is(
    edger_to_df(glmlrt_single_contrast),
    "data.frame",
    label = "glmLRT with single coef"
  )
  expect_is(
    edger_to_df(glmlrt_anova),
    "data.frame",
    label = "glmLRT with multiple coefs"
  )
  expect_is(
    edger_to_df(glmtreat_single_contrast),
    "data.frame",
    label = "glmTreat with single coef"
  )
  expect_is(edger_to_df(et_single_contrast), "data.frame", label = "exactTest")
})

make_fit <- function() {
  set.seed(1)
  mu <- matrix(100, 200, 6)
  mu[1:20, 4:6] <- 800
  counts <- matrix(
    rnbinom(200 * 6, mu = mu, size = 20), 200,
    dimnames = list(paste0("g", 1:200), paste0("s", 1:6))
  )
  group <- factor(rep(c("A", "B"), each = 3))
  d <- edgeR::calcNormFactors(edgeR::DGEList(counts, group = group))
  design <- model.matrix(~group)
  d <- edgeR::estimateDisp(d, design)
  list(d = d, design = design, qlf = edgeR::glmQLFTest(edgeR::glmQLFit(d, design), coef = 2))
}

test_that("edger_to_df returns all rows in input order with BH FDR", {
  skip_if_not_installed("edgeR")
  o <- make_fit()
  d <- o$d
  qlf <- o$qlf
  res <- edger_to_df(qlf)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), nrow(d))
  expect_identical(colnames(res)[1], "feature_id")
  expect_identical(res$feature_id, rownames(d))
  expect_identical(rownames(res), as.character(seq_len(nrow(d))))
  expect_equal(res$logFC, qlf$table$logFC)
  expect_equal(res$PValue, qlf$table$PValue)
  expect_equal(res$FDR, p.adjust(qlf$table$PValue, "BH"))
  # simulated 8-fold up genes are recovered
  expect_true(all(res$logFC[1:20] > 2))
  expect_true(all(res$FDR[1:20] < 0.05))
})

test_that("... is forwarded to topTags", {
  skip_if_not_installed("edgeR")
  o <- make_fit()
  qlf <- o$qlf
  res <- edger_to_df(qlf, adjust.method = "bonferroni")
  expect_true("FWER" %in% colnames(res))
  expect_equal(res$FWER, p.adjust(qlf$table$PValue, "bonferroni"))
  keep <- p.adjust(qlf$table$PValue, "BH") <= 0.05
  r2 <- edger_to_df(qlf, p.value = 0.05)
  expect_gt(sum(keep), 0)
  expect_identical(r2$feature_id, rownames(o$d)[keep])
  expect_true(all(r2$FDR <= 0.05))
})

test_that("Filtering that removes every row gives a 0-row data.frame", {
  skip_if_not_installed("edgeR")
  o <- make_fit()
  qlf <- o$qlf
  full <- edger_to_df(qlf)
  res <- edger_to_df(qlf, p.value = 1e-300)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 0)
  expect_identical(colnames(res), colnames(full))
  res2 <- edger_to_df(qlf, p.value = 1e-300, adjust.method = "bonferroni")
  expect_equal(nrow(res2), 0)
  expect_true("FWER" %in% colnames(res2))
})

test_that("Other edgeR result types give FDR = BH(PValue) in feature order", {
  skip_if_not_installed("edgeR")
  o <- make_fit()
  d <- o$d
  et <- edgeR::exactTest(d)
  res <- edger_to_df(et)
  expect_identical(res$feature_id, rownames(d))
  expect_equal(res$FDR, p.adjust(et$table$PValue, "BH"))
  design3 <- model.matrix(~ 0 + factor(rep(1:3, each = 2)))
  f3 <- edgeR::glmQLFit(edgeR::estimateDisp(d, design3), design3)
  an <- edgeR::glmQLFTest(f3, contrast = cbind(c(-1, 1, 0), c(-1, 0, 1)))
  r3 <- edger_to_df(an)
  expect_equal(sum(grepl("^logFC", colnames(r3))), 2)
  expect_equal(r3$FDR, p.adjust(an$table$PValue, "BH"))
  lrt <- edgeR::glmLRT(edgeR::glmFit(d, o$design), coef = 2)
  expect_equal(edger_to_df(lrt)$FDR, p.adjust(lrt$table$PValue, "BH"))
  tr <- edgeR::glmTreat(edgeR::glmQLFit(d, o$design), coef = 2, lfc = 1)
  expect_equal(edger_to_df(tr)$FDR, p.adjust(tr$table$PValue, "BH"))
})

test_that("An existing feature_id annotation column is preserved", {
  skip_if_not_installed("edgeR")
  o <- make_fit()
  d2 <- o$d
  d2$genes <- data.frame(
    feature_id = paste0("ENSG", seq_len(nrow(d2))), row.names = rownames(d2)
  )
  q2 <- edgeR::glmQLFTest(edgeR::glmQLFit(d2, o$design), coef = 2)
  res <- edger_to_df(q2)
  expect_equal(anyDuplicated(colnames(res)), 0)
  expect_identical(res$feature_id.1, rownames(o$d))
  expect_identical(res$feature_id, paste0("ENSG", seq_len(nrow(o$d))))
  expect_equal(res$FDR, p.adjust(o$qlf$table$PValue, "BH"))
})
