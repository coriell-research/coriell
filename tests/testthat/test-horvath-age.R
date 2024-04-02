test_that("Inverse transformation returns same values", {
  age <- 0:100
  transformed <- horvath_age(age)
  transformed2 <- horvath_age(age, adult_age = 30)
  
  expect_equal(age, horvath_age(transformed, inverse=TRUE))
  expect_equal(age, horvath_age(transformed2, inverse=TRUE, adult_age=30))
})
