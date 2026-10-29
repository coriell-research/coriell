context("Test summarize_dge function")
library(coriell)


# Create 100 genes
# 50 will have logFC > 0
# 50 will have logFC < 0
# 25 of the up-regulated will have significant FDR values (< 0.0001)
# 25 of the down-regulated will have significant FDR values (< 0.0001)
df <- data.frame(
  feature_id = paste("gene", 1:100, sep = "."),
  logFC = c(
    runif(50, min = -8, max = -0.01),
    runif(50, min = 0.01, max = 8)
  ),
  FDR = c(
    runif(25, min = 1e-12, max = 1e-4),
    runif(25, min = 0.06, max = 1),
    runif(25, min = 1e-12, max = 1e-4),
    runif(25, min = 0.06, max = 1)
  )
)
res <- summarize_dge(df)

test_that("Check accuracy of summary results", {
  expect_equal(
    res[res$Direction == "Up", "N", drop = TRUE],
    25,
    label = "Up-regulated genes equal 25"
  )
  expect_equal(
    res[res$Direction == "Down", "N", drop = TRUE],
    25,
    label = "Down-regulated genes equal 25"
  )
  expect_equal(
    res[res$Direction == "Unperturbed", "N", drop = TRUE],
    50,
    label = "Unperturbed genes equal 50"
  )
  expect_equal(
    identical(as.character(res$Direction), c("Up", "Down", "Unperturbed")),
    TRUE,
    "All DE factors present"
  )
})

# change colnames of df and check function using new colnames
colnames(df) <- c("gene_id", "logFoldChange", "p.adjust")
res2 <- summarize_dge(df, fdr_col = "p.adjust", lfc_col = "logFoldChange")

test_that("Check accuracy of summary results with different colnames", {
  expect_equal(
    res2[res2$Direction == "Up", "N", drop = TRUE],
    25,
    label = "Up-regulated genes equal 25"
  )
  expect_equal(
    res2[res2$Direction == "Down", "N", drop = TRUE],
    25,
    label = "Down-regulated genes equal 25"
  )
  expect_equal(
    res2[res2$Direction == "Unperturbed", "N", drop = TRUE],
    50,
    label = "Unperturbed genes equal 50"
  )
  expect_equal(
    identical(as.character(res2$Direction), c("Up", "Down", "Unperturbed")),
    TRUE,
    "All DE factors present"
  )
})

test_that("Hand-computed counts respect strict thresholds", {
  df <- data.frame(
    FDR   = c(0.01, 0.01, 0.05, 0.01, 0.01, 0.20, 0.01, 0.01),
    logFC = c(2,    -3,   4,    0,    1,    5,    -1,   0.5)
  )
  res <- summarize_dge(df, lfc = 0)
  expect_identical(levels(res$Direction), c("Up", "Down", "Unperturbed"))
  # FDR == 0.05 and logFC == 0 are not significant
  expect_equal(res$N, c(3, 2, 3))
  expect_equal(res$Percent, c(37.5, 25, 37.5))
  # |logFC| == 1 is not significant
  expect_equal(summarize_dge(df, lfc = 1)$N, c(1, 1, 6))
  expect_equal(summarize_dge(df, fdr = 0.1, lfc = 1)$N, c(2, 1, 5))
})

test_that("NA values are counted as unperturbed", {
  df <- data.frame(FDR = c(NA, 0.01, 0.01, 0.01), logFC = c(1, NA, 2, -2))
  expect_equal(summarize_dge(df)$N, c(1, 1, 2))
})

test_that("summarize_dge agrees with base R on package data", {
  df <- GSE161650_de
  for (th in list(c(0.05, 0), c(0.01, 1), c(1, 0))) {
    res <- summarize_dge(df, fdr = th[1], lfc = th[2])
    sig <- !is.na(df$FDR) & df$FDR < th[1] & !is.na(df$logFC) & abs(df$logFC) > th[2]
    expect_equal(res$N[1], sum(sig & df$logFC > 0))
    expect_equal(res$N[2], sum(sig & df$logFC < 0))
    expect_equal(sum(res$N), nrow(df))
    expect_equal(res$Percent, round(100 * res$N / nrow(df), 2))
  }
  a <- summarize_dge(df, fdr = 0.05, lfc = 0)
  b <- summarize_dge(df, fdr = 0.01, lfc = 1)
  expect_true(all(b$N[1:2] <= a$N[1:2]))
})

test_that("Custom column names match renaming to the defaults", {
  df <- data.frame(
    id = letters[1:5], lfc = c(1, -1, 2, -2, 0),
    padj = c(0.01, 0.01, 0.5, 0.001, 0.01)
  )
  r1 <- summarize_dge(df, fdr_col = "padj", lfc_col = "lfc")
  r2 <- summarize_dge(setNames(df, c("id", "logFC", "FDR")))
  expect_identical(r1, r2)
  expect_equal(r1$N, c(1, 2, 2))
  expect_error(summarize_dge(df), "fdr_col")
})
