test_that("UMAP dispatch works properly", {
  col <- 20
  row <- 1000
  mat1 <- matrix(
    rexp(col * row, rate = 0.1),
    ncol = col
  )
  rownames(mat1) <- paste0('gene', 1:nrow(mat1))
  colnames(mat1) <- paste0('sample', 1:ncol(mat1))

  mat2 <- matrix(
    rexp(col * row, rate = 0.1),
    ncol = col
  )
  rownames(mat2) <- paste0('gene', 1:nrow(mat2))
  colnames(mat2) <- paste0('sample', (ncol(mat1) + 1):(ncol(mat1) + ncol(mat2)))

  mat <- cbind(mat1, mat2)

  metadata <- data.frame(row.names = colnames(mat))
  metadata$Group <- rep(NA, ncol(mat))
  metadata$Group[seq(1, 40, 2)] <- 'A'
  metadata$Group[seq(2, 40, 2)] <- 'B'
  metadata$CRP <- sample.int(100, size = ncol(mat), replace = TRUE)
  metadata$ESR <- sample.int(100, size = ncol(mat), replace = TRUE)

  p <- PCAtools::pca(mat, metadata = metadata, center = TRUE, scale = TRUE)
  pr <- prcomp(t(mat), center = TRUE, scale = TRUE)

  # on PCA object
  u <- UMAP(p)
  expect_is(u, "data.frame")
  expect_equal(nrow(u), ncol(mat))

  # on prcomp object
  u2 <- UMAP(pr, metadata = metadata)
  expect_is(u2, "data.frame")
  expect_equal(nrow(u2), ncol(mat))

  # on matrix object
  u3 <- UMAP(t(mat), metadata = metadata)
  expect_is(u3, "data.frame")
  expect_equal(nrow(u3), ncol(mat))

  # on data.frame
  df <- data.frame(t(mat))
  u4 <- UMAP(df, metadata = metadata)
  expect_is(u4, "data.frame")
  expect_equal(nrow(u4), ncol(mat))

  # on distance matrix
  d <- dist(t(mat))
  u5 <- UMAP(d, metadata = metadata)
  expect_is(u5, "data.frame")
  expect_equal(nrow(u5), ncol(mat))
})
