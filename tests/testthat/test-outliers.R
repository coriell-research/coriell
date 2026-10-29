test_that("outliers_by_iqr matches hand-computed Tukey fences", {
  # type-7 quartiles: q1 = 1.5, q3 = 5, IQR = 3.5 -> fences -3.75, 10.25
  x <- c(1, 1, 2, 2, 4, 6, 11)
  expect_identical(outliers_by_iqr(x), c(rep(FALSE, 6), TRUE))
  expect_identical(outliers_by_iqr(x, direction = "low"), rep(FALSE, 7))
  expect_identical(outliers_by_iqr(x, direction = "high"), c(rep(FALSE, 6), TRUE))
  # 11 is (11 - 5) / 3.5 = 1.714 IQRs above q3
  expect_false(any(outliers_by_iqr(x, threshold = 1.72)))
  expect_true(outliers_by_iqr(x, threshold = 1.71)[7])
})

test_that("outliers_by_iqr uses a strict inequality at the fence", {
  # q1 = 1.25, q3 = 3.75, IQR = 2.5 -> upper fence exactly 7.5
  expect_false(outliers_by_iqr(c(0, 1, 2, 3, 4, 7.5))[6])
  expect_true(outliers_by_iqr(c(0, 1, 2, 3, 4, 7.6))[6])
})

test_that("outliers_by_iqr agrees with boxplot.stats and keeps NA in place", {
  x <- c(1:20, -30, 50, NA)
  both <- outliers_by_iqr(x)
  expect_length(both, length(x))
  expect_true(is.na(both[23]))
  expect_identical(which(both), c(21L, 22L))
  expect_identical(which(outliers_by_iqr(x, direction = "low")), 21L)
  expect_identical(which(outliers_by_iqr(x, direction = "high")), 22L)
  expect_setequal(boxplot.stats(x)$out, x[which(both)])
})

test_that("outliers_by_mad matches hand-computed fences", {
  # median = 2, raw MAD = 1, scaled MAD = 1.4826
  # upper fences: 6.448 (k = 3), 4.965 (k = 2), 3.483 (k = 1)
  x <- c(1, 1, 2, 2, 4, 6, 11)
  expect_identical(outliers_by_mad(x), c(rep(FALSE, 6), TRUE))
  expect_identical(outliers_by_mad(x, threshold = 2), c(rep(FALSE, 5), TRUE, TRUE))
  expect_identical(outliers_by_mad(x, threshold = 1), c(rep(FALSE, 4), TRUE, TRUE, TRUE))
  expect_identical(outliers_by_mad(x, threshold = 1, direction = "low"), rep(FALSE, 7))
})

test_that("outliers_by_mad handles NA, direction and affine transforms", {
  set.seed(2)
  z <- c(rnorm(30), NA, 20, -20)
  r <- outliers_by_mad(z)
  expect_true(is.na(r[31]))
  expect_identical(which(r), c(32L, 33L))
  expect_identical(which(outliers_by_mad(z, direction = "low")), 33L)
  expect_identical(which(outliers_by_mad(z, direction = "high")), 32L)
  expect_identical(outliers_by_mad(5 * z + 100), r)
})

test_that("outliers_by_z hand example", {
  # c(-1, 0, 1) has mean 0 and sd 1, so z = x
  x <- c(-1, 0, 1)
  expect_identical(outliers_by_z(x, threshold = 0.5), c(TRUE, FALSE, TRUE))
  expect_identical(outliers_by_z(x, threshold = 0.5, direction = "high"), c(FALSE, FALSE, TRUE))
  expect_identical(outliers_by_z(x, threshold = 0.5, direction = "low"), c(TRUE, FALSE, FALSE))
  expect_identical(outliers_by_z(x, threshold = 1), c(FALSE, FALSE, FALSE))
})

test_that("outliers_by_z matches explicit z-scores and symmetries", {
  set.seed(3)
  y <- c(rnorm(100), 50)
  z <- (y - mean(y)) / sd(y)
  expect_identical(outliers_by_z(y), abs(z) > 3)
  expect_identical(which(outliers_by_z(y)), 101L)
  expect_identical(outliers_by_z(-y, direction = "low"), outliers_by_z(y, direction = "high"))
  expect_identical(outliers_by_z(3 * y - 7), outliers_by_z(y))

  r <- outliers_by_z(c(y, NA))
  expect_true(is.na(r[102]))
  expect_identical(which(r), 101L)
})
