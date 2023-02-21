test_that("Random RGB catches bad input", {
  expect_error(random_rgb_palette("A"))
  expect_error(random_rgb_palette(factor("A")))
  expect_error(random_rgb_palette(c(1, 2, 3)))
  expect_error(random_rgb_palette(c(1, 2, 3)))
  expect_error(random_rgb_palette(-1))
})

test_that("Distinct RGB catches bad input", {
  expect_error(distinct_rgb_palette("A"))
  expect_error(distinct_rgb_palette(factor("A")))
  expect_error(distinct_rgb_palette(c(1, 2, 3)))
  expect_error(distinct_rgb_palette(c(1, 2, 3)))
  expect_error(distinct_rgb_palette(-1))
})