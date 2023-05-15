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
  expect_equal(res[res$dge == "up", "n", drop = TRUE], 25, label = "Up-regulated genes equal 25")
  expect_equal(res[res$dge == "down", "n", drop = TRUE], 25, label = "Down-regulated genes equal 25")
  expect_equal(res[res$dge == "non-dge", "n", drop = TRUE], 50, label = "Non-DE genes equal 50")
  expect_equal(all(unique(res$dge) == c("up", "down", "non-dge")), TRUE, "All DE factors present")
})

# change colnames of df and check function using new colnames
colnames(df) <- c("gene_id", "logFoldChange", "p.adjust")
res2 <- summarize_dge(df, fdr_col = p.adjust, lfc_col = logFoldChange)

test_that("Check accuracy of summary results with different colnames", {
  expect_equal(res2[res2$dge == "up", "n", drop = TRUE], 25, label = "Up-regulated genes equal 25")
  expect_equal(res2[res2$dge == "down", "n", drop = TRUE], 25, label = "Down-regulated genes equal 25")
  expect_equal(res2[res2$dge == "non-dge", "n", drop = TRUE], 50, label = "Non-DE genes equal 50")
  expect_equal(all(unique(res2$dge) == c("up", "down", "non-dge")), TRUE, "All DE factors present")
})
