test_that("quickmap clusters with the documented defaults", {
  skip_if_not_installed("pheatmap")
  set.seed(1)
  x <- matrix(rnorm(40 * 6), 40, dimnames = list(paste0("g", 1:40), paste0("s", 1:6)))
  x[, 4:6] <- x[, 4:6] + 2
  p <- quickmap(x, silent = TRUE)
  expect_s3_class(p, "pheatmap")
  expect_identical(p$tree_col$method, "ward.D2")
  expect_identical(p$tree_row$method, "ward.D2")
  expect_identical(p$tree_col$dist.method, "euclidean")
  # row-scaled data clustered independently
  xs <- t(scale(t(x)))
  hc <- hclust(dist(t(xs)), method = "ward.D2")
  expect_equal(p$tree_col$merge, hc$merge)
  expect_equal(p$tree_col$height, hc$height)
  hr <- hclust(as.dist(1 - cor(t(xs))), method = "ward.D2")
  expect_equal(p$tree_row$merge, hr$merge)
  expect_equal(p$tree_row$height, hr$height)
  k <- cutree(p$tree_col, 2)
  expect_true(all(k[1:3] == k[[1]]))
  expect_true(all(k[4:6] != k[[1]]))
})

test_that("quickmap user args override defaults", {
  skip_if_not_installed("pheatmap")
  set.seed(1)
  x <- matrix(rnorm(40 * 6), 40, dimnames = list(paste0("g", 1:40), paste0("s", 1:6)))
  p <- quickmap(x, cluster_rows = FALSE, cluster_cols = FALSE, silent = TRUE)
  expect_true(is.na(p$tree_row)[1])
  expect_true(is.na(p$tree_col)[1])
  p2 <- quickmap(x, clustering_method = "complete", clustering_distance_cols = "manhattan", silent = TRUE)
  expect_identical(p2$tree_col$method, "complete")
  expect_identical(p2$tree_col$dist.method, "manhattan")
})

test_that("quickmap remove_var keeps the highest-variance rows", {
  skip_if_not_installed("pheatmap")
  set.seed(2)
  x <- matrix(rnorm(100 * 6), 100, dimnames = list(paste0("g", 1:100), paste0("s", 1:6)))
  x <- x * (1:100)
  p <- suppressMessages(quickmap(x, remove_var = 0.25, silent = TRUE))
  v <- apply(x, 1, var)
  expect_setequal(p$tree_row$labels, names(sort(v, decreasing = TRUE))[1:75])
})

test_that("quickmap adds a column annotation from metadata", {
  skip_if_not_installed("pheatmap")
  set.seed(1)
  x <- matrix(rnorm(40 * 6), 40, dimnames = list(paste0("g", 1:40), paste0("s", 1:6)))
  md <- data.frame(group = rep(c("a", "b"), each = 3), row.names = colnames(x))
  expect_true("col_annotation" %in% quickmap(x, metadata = md, silent = TRUE)$gtable$layout$name)
  expect_false("col_annotation" %in% quickmap(x, silent = TRUE)$gtable$layout$name)
  expect_error(
    quickmap(x, metadata = md[6:1, , drop = FALSE], silent = TRUE),
    "rownames\\(metadata\\) != colnames\\(x\\)"
  )
})
