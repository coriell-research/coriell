test_that("get_axis_limits returns data ranges and coord limits", {
  d <- data.frame(x = c(1, 5, 3), y = c(-2, 10, 4))
  p <- ggplot2::ggplot(d, ggplot2::aes(x, y)) + ggplot2::geom_point()
  expect_equal(get_axis_limits(p), list(x_min = 1, x_max = 5, y_min = -2, y_max = 10))
  expect_equal(
    get_axis_limits(p + ggplot2::coord_cartesian(xlim = c(0, 20))),
    list(x_min = 0, x_max = 20, y_min = -2, y_max = 10)
  )
  expect_equal(
    get_axis_limits(p + ggplot2::coord_cartesian(xlim = c(0, 20), ylim = c(-5, 5))),
    list(x_min = 0, x_max = 20, y_min = -5, y_max = 5)
  )
  expect_equal(get_axis_limits(p + ggplot2::geom_hline(yintercept = 50))$y_max, 50)

  set.seed(1)
  d2 <- data.frame(x = rnorm(50), y = runif(50))
  l <- get_axis_limits(ggplot2::ggplot(d2, ggplot2::aes(x, y)) + ggplot2::geom_point())
  expect_equal(c(l$x_min, l$x_max, l$y_min, l$y_max), c(range(d2$x), range(d2$y)))
})
