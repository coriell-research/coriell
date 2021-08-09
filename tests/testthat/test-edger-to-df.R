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
design <- model.matrix(~0 + group)
y <- DGEList(counts = x, group = group)
y <- calcNormFactors(y)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)

# create results objects for all test types
glmqlf_single_contrast <- glmQLFTest(fit, coef = 2)
glmqlf_anova <- glmQLFTest(fit, coef = 1:3)
glmlrt_single_contrast <- glmLRT(fit, coef = 2)
glmlrt_anova <- glmLRT(fit, coef = 1:3)
glmtreat_single_contrast <- glmTreat(fit, contrast = c(-1, 1, 0), lfc = log2(1.2))
et_single_contrast <- exactTest(y)


# tests -------------------------------------------------------------------

test_that("Convert results to data.frame", {
  expect_is(edger_to_df(glmqlf_single_contrast), "data.frame", label = "glmQLF test with single coef")
  expect_is(edger_to_df(glmqlf_anova), "data.frame", label = "glmQLF with multiple coefs")
  expect_is(edger_to_df(glmlrt_single_contrast), "data.frame", label = "glmLRT with single coef")
  expect_is(edger_to_df(glmlrt_anova), "data.frame", label = "glmLRT with multiple coefs")
  expect_is(edger_to_df(glmtreat_single_contrast), "data.frame", label = "glmTreat with single coef")
  expect_is(edger_to_df(et_single_contrast), "data.frame", label = "exactTest")
  }
)
