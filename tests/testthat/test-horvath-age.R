test_that("Inverse transformation returns same values", {
  age <- 0:100
  transformed <- horvath_age(age)
  transformed2 <- horvath_age(age, adult_age = 30)

  expect_equal(age, horvath_age(transformed, inverse = TRUE))
  expect_equal(age, horvath_age(transformed2, inverse = TRUE, adult_age = 30))
})

test_that("horvath_age matches the published transformation", {
  # F(a) = log(a + 1) - log(A + 1) for a <= A, (a - A) / (A + 1) otherwise
  expect_equal(horvath_age(c(0, 20, 41, 62)), c(-log(21), 0, 1, 2))
  expect_equal(horvath_age(c(-0.5, 9)), c(log(0.5) - log(21), log(10) - log(21)))
  expect_equal(horvath_age(c(0, 15, 31), adult_age = 15), c(-log(16), 0, 1))
  expect_equal(horvath_age(c(-log(21), 0, 1, 2), inverse = TRUE), c(0, 20, 41, 62))
  expect_equal(horvath_age(-log(16), adult_age = 15, inverse = TRUE), 0)
})

test_that("horvath_age is monotone, continuous at adult_age and keeps shape", {
  a <- seq(-0.9, 120, by = 0.1)
  f <- horvath_age(a)
  expect_true(all(diff(f) > 0))
  expect_equal(horvath_age(f, inverse = TRUE), a)
  t <- seq(-3, 5, by = 0.01)
  expect_equal(horvath_age(horvath_age(t, inverse = TRUE)), t)
  # both one-sided derivatives at 20 equal 1 / 21
  eps <- 1e-6
  expect_equal((horvath_age(20) - horvath_age(20 - eps)) / eps, 1 / 21, tolerance = 1e-5)
  expect_equal((horvath_age(20 + eps) - horvath_age(20)) / eps, 1 / 21, tolerance = 1e-5)
  expect_equal(horvath_age(c(5, NA, 30)), c(log(6) - log(21), NA, 10 / 21))
  expect_equal(horvath_age(matrix(c(0, 20, 41, 62), 2)), matrix(c(-log(21), 0, 1, 2), 2))
})
