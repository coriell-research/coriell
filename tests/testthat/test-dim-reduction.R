# Independent reconstruction from prcomp keeping the given PCs
recon_oracle <- function(X, keep, ...) {
  p <- prcomp(X, ...)
  r <- p$x[, keep, drop = FALSE] %*% t(p$rotation[, keep, drop = FALSE])
  if (!isFALSE(p$scale)) r <- sweep(r, 2, p$scale, "*")
  if (!isFALSE(p$center)) r <- sweep(r, 2, p$center, "+")
  r
}

test_that("remove_components matches manual reconstruction", {
  X <- as.matrix(USArrests)
  r <- remove_components(
    USArrests,
    components = 1:2,
    center = TRUE,
    scale. = TRUE
  )
  expect_s3_class(r, "data.frame")
  expect_identical(dimnames(r), dimnames(USArrests))
  expect_equal(
    as.matrix(r),
    recon_oracle(X, 3:4, center = TRUE, scale. = TRUE),
    check.attributes = FALSE
  )
})

test_that("remove_components honours prcomp args", {
  X <- as.matrix(USArrests)
  expect_equal(
    as.matrix(remove_components(USArrests)),
    recon_oracle(X, 2:4),
    check.attributes = FALSE
  )
  expect_equal(
    as.matrix(remove_components(USArrests, 1, center = FALSE)),
    recon_oracle(X, 2:4, center = FALSE),
    check.attributes = FALSE
  )
  expect_equal(
    as.matrix(remove_components(USArrests, 1:3, scale. = TRUE)),
    recon_oracle(X, 4, scale. = TRUE),
    check.attributes = FALSE
  )
  expect_equal(
    as.matrix(remove_components(USArrests, 1, scale. = c(1, 2, 3, 4))),
    recon_oracle(X, 2:4, scale. = c(1, 2, 3, 4)),
    check.attributes = FALSE
  )
})

test_that("remove_components satisfies linear-algebra invariants", {
  p <- prcomp(USArrests, scale. = TRUE)
  Z <- scale(
    as.matrix(remove_components(USArrests, 1, scale. = TRUE)),
    center = p$center,
    scale = p$scale
  )
  # result is orthogonal to the removed PC and has rank p - 1
  expect_equal(unname(drop(Z %*% p$rotation[, 1])), rep(0, 50), tolerance = 1e-10)
  expect_equal(qr(Z)$rank, 3L)
  Z2 <- scale(
    as.matrix(remove_components(USArrests, 2:3, scale. = TRUE)),
    center = p$center,
    scale = p$scale
  )
  expect_equal(qr(Z2)$rank, 2L)
  # column means kept
  expect_equal(
    unname(colMeans(remove_components(USArrests, 1))),
    unname(colMeans(USArrests))
  )
  # complementary removals add back to the original data
  a <- as.matrix(remove_components(USArrests, 1:2))
  b <- as.matrix(remove_components(USArrests, 3:4))
  expect_equal(
    a + b - matrix(colMeans(USArrests), 50, 4, byrow = TRUE),
    as.matrix(USArrests),
    check.attributes = FALSE,
    tolerance = 1e-8
  )
})

test_that("removing all PCs gives column means; removing none gives the data", {
  # regression: the guard refused to remove the last PC
  r <- remove_components(USArrests, 1:4)
  expect_equal(
    as.matrix(r),
    matrix(colMeans(USArrests), 50, 4, byrow = TRUE),
    check.attributes = FALSE
  )
  expect_no_error(remove_components(USArrests, 4))
  expect_equal(
    as.matrix(remove_components(USArrests, integer(0), scale. = TRUE)),
    as.matrix(USArrests),
    check.attributes = FALSE
  )
})

test_that("remove_components keeps rownames when one PC remains", {
  # regression: res$x[, keep] dropped to a vector
  r <- remove_components(USArrests, 1:3)
  expect_identical(rownames(r), rownames(USArrests))
  expect_identical(colnames(r), colnames(USArrests))
})

test_that("remove_components errors on nonexistent PCs", {
  # regression: wide data has fewer PCs than columns; PC 6 silently ignored
  set.seed(1)
  W <- matrix(rnorm(50), 5, 10)
  expect_equal(ncol(prcomp(W)$x), 5)
  expect_error(remove_components(W, 6), "number of principal components")
  expect_error(remove_components(USArrests, 5), "number of principal components")
})

test_that("associate_components matches base tests", {
  set.seed(42)
  x <- matrix(rnorm(200 * 8), 200, 8, dimnames = list(NULL, paste0("s", 1:8)))
  md <- data.frame(
    age = c(30, 80, 34, 30, 80, 40, 55, 62),
    trt = factor(rep(c("T", "C"), each = 4)),
    cls = rep(c("A", "B", "C", "D"), 2),
    flag = rep(c(TRUE, FALSE), 4),
    row.names = colnames(x)
  )
  res <- associate_components(x, md, N = 3)
  R <- prcomp(x)$rotation
  expected <- rbind(
    sapply(1:3, function(j) cor.test(md$age, R[, j])$p.value),
    sapply(1:3, function(j) {
      wilcox.test(R[md$trt == "C", j], R[md$trt == "T", j])$p.value
    }),
    sapply(1:3, function(j) kruskal.test(R[, j], factor(md$cls))$p.value),
    sapply(1:3, function(j) wilcox.test(R[!md$flag, j], R[md$flag, j])$p.value)
  )
  dimnames(expected) <- list(c("age", "trt", "cls", "flag"), paste0("PC", 1:3))
  expect_s3_class(res, "data.frame")
  expect_equal(as.matrix(res), expected)

  # N larger than number of PCs is capped
  expect_equal(dim(associate_components(x, md)), c(4L, 8L))
  # invariant to sign and sample permutation
  expect_equal(associate_components(-x, md, N = 3), res)
  o <- c(3, 8, 1, 5, 2, 7, 4, 6)
  expect_equal(associate_components(x[, o], md[o, ], N = 3), res)
  expect_equal(associate_components(as.data.frame(x), md, N = 3), res)
  # ... passed to prcomp
  R2 <- prcomp(x, scale. = TRUE)$rotation
  expect_equal(
    associate_components(x, md, N = 2, scale. = TRUE)["age", "PC1"],
    cor.test(md$age, R2[, 1])$p.value
  )
  expect_error(associate_components(x, md[8:1, ]), "do not match")
  md$age[2] <- NA
  expect_error(associate_components(x, md), "NA values")
})

test_that("associate_components gives exact Wilcoxon p for perfect separation", {
  set.seed(1)
  x <- matrix(rnorm(500 * 8), 500, 8, dimnames = list(NULL, paste0("s", 1:8)))
  x[1:100, 1:4] <- x[1:100, 1:4] + 5
  md <- data.frame(g = rep(c("A", "B"), each = 4), row.names = colnames(x))
  # two-sided exact p-value for complete separation of 4 vs 4
  expect_equal(associate_components(x, md, N = 2)["g", "PC1"], 2 / choose(8, 4))
})

test_that("associate_components works with N = 1", {
  # regression: rotation[, 1:1] dropped to a vector
  set.seed(1)
  x <- matrix(rnorm(400), 50, 8, dimnames = list(NULL, paste0("s", 1:8)))
  md <- data.frame(a = 1:8, row.names = colnames(x))
  r <- associate_components(x, md, N = 1)
  expect_equal(dim(r), c(1L, 1L))
  expect_equal(dimnames(r), list("a", "PC1"))
  expect_equal(r$PC1, cor.test(1:8, prcomp(x)$rotation[, 1])$p.value)
})

test_that("UMAP methods equal direct umap::umap", {
  skip_if_not_installed("umap")
  skip_if_not_installed("PCAtools")
  set.seed(1)
  m <- matrix(
    rnorm(30 * 5),
    30,
    5,
    dimnames = list(paste0("s", 1:30), paste0("f", 1:5))
  )
  md <- data.frame(Group = factor(rep(c("A", "B"), 15)), row.names = rownames(m))
  cfg <- umap::umap.defaults
  cfg$n_neighbors <- 5
  cfg$random_state <- 123
  cfg$n_epochs <- 50
  ref <- umap::umap(m, cfg)$layout
  run <- function(x, ...) {
    UMAP(x, ..., n_neighbors = 5, random_state = 123, n_epochs = 50)
  }

  um <- run(m, metadata = md)
  expect_s3_class(um, "data.frame")
  expect_equal(colnames(um), c("UMAP1", "UMAP2", "Group"))
  expect_equal(rownames(um), rownames(m))
  expect_s3_class(um$Group, "factor")
  expect_equal(unname(as.matrix(um[, 1:2])), unname(ref))
  expect_equal(
    unname(as.matrix(run(as.data.frame(m), metadata = md)[, 1:2])),
    unname(ref)
  )

  ud <- run(dist(m), metadata = md)
  cfgd <- cfg
  cfgd$input <- "dist"
  expect_equal(
    unname(as.matrix(ud[, 1:2])),
    unname(umap::umap(as.matrix(dist(m)), cfgd)$layout)
  )

  pr <- prcomp(m)
  expect_equal(
    unname(as.matrix(run(pr, metadata = md)[, 1:2])),
    unname(umap::umap(pr$x, cfg)$layout)
  )

  p <- PCAtools::pca(t(m), metadata = md)
  up <- run(p)
  expect_equal(colnames(up), c("UMAP1", "UMAP2", "Group"))
  expect_equal(
    unname(as.matrix(up[, 1:2])),
    unname(umap::umap(as.matrix(p$rotated), cfg)$layout)
  )

  expect_equal(colnames(run(m, n_components = 3)), paste0("UMAP", 1:3))
  cfgm <- cfg
  cfgm$metric <- "manhattan"
  expect_equal(
    unname(as.matrix(run(m, metric = "manhattan"))),
    unname(umap::umap(m, cfgm)$layout)
  )

  expect_error(run(m, metadata = md[30:1, , drop = FALSE]), "rownames")
  expect_error(run(dist(m), metadata = md[30:1, , drop = FALSE]), "Labels")
  expect_error(run(pr, metadata = md[30:1, , drop = FALSE]), "rownames")
})

test_that("UMAP returns a data.frame without metadata", {
  # regression: matrix was returned when metadata was NULL
  skip_if_not_installed("umap")
  skip_if_not_installed("PCAtools")
  set.seed(1)
  m <- matrix(
    rnorm(30 * 5),
    30,
    5,
    dimnames = list(paste0("s", 1:30), paste0("f", 1:5))
  )
  run <- function(x, ...) {
    UMAP(x, ..., n_neighbors = 5, random_state = 123, n_epochs = 50)
  }
  inputs <- list(
    matrix = m,
    data.frame = as.data.frame(m),
    dist = dist(m),
    prcomp = prcomp(m),
    pca = PCAtools::pca(t(m))
  )
  for (nm in names(inputs)) {
    u <- run(inputs[[nm]])
    expect_s3_class(u, "data.frame")
    expect_equal(colnames(u), c("UMAP1", "UMAP2"))
    expect_equal(rownames(u), rownames(m))
  }
})

test_that("UMAP keeps metadata column names as given", {
  # regression: data.frame() check.names mangled 'my var' -> 'my.var'
  skip_if_not_installed("umap")
  skip_if_not_installed("PCAtools")
  set.seed(1)
  m <- matrix(
    rnorm(30 * 5),
    30,
    5,
    dimnames = list(paste0("s", 1:30), paste0("f", 1:5))
  )
  md <- data.frame(
    `my var` = rep(c("A", "B"), 15),
    `1st` = 1:30,
    row.names = rownames(m),
    check.names = FALSE
  )
  run <- function(x, ...) {
    UMAP(x, ..., n_neighbors = 5, random_state = 123, n_epochs = 50)
  }
  expected <- c("UMAP1", "UMAP2", "my var", "1st")
  expect_equal(colnames(run(m, metadata = md)), expected)
  expect_equal(colnames(run(as.data.frame(m), metadata = md)), expected)
  expect_equal(colnames(run(dist(m), metadata = md)), expected)
  expect_equal(colnames(run(prcomp(m), metadata = md)), expected)
  expect_equal(colnames(run(PCAtools::pca(t(m), metadata = md))), expected)
})
