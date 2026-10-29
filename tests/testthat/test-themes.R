is_blank <- function(x) inherits(x, c("element_blank", "ggplot2::element_blank"))

test_that("theme_coriell sets the documented overrides", {
  th <- theme_coriell()
  expect_s3_class(th, "theme")
  expect_true(attr(th, "complete"))
  expect_identical(th$legend.position, "bottom")
  expect_identical(th$text$family, "Helvetica")
  expect_equal(th$plot.title$size, 22)
  expect_identical(th$plot.title$face, "bold")
  expect_equal(th$plot.title$hjust, 0)
  expect_equal(th$plot.subtitle$size, 18)
  expect_equal(th$legend.title$size, 16)
  expect_identical(th$legend.title$face, "bold")
  expect_equal(th$legend.text$size, 14)
  expect_equal(th$axis.title$size, 16)
  expect_equal(th$axis.text.x$size, 16)
  expect_equal(th$axis.text.y$size, 16)
  for (el in c("strip.text.x", "strip.text.y")) {
    expect_equal(th[[el]]$size, 16, info = el)
    expect_identical(th[[el]]$colour, "black", info = el)
    expect_identical(th[[el]]$face, "bold", info = el)
  }
  expect_true(is_blank(th$strip.background))
  expect_true(is_blank(th$panel.border))
})

test_that("Elements not overridden match theme_light", {
  th <- theme_coriell()
  tl <- ggplot2::theme_light()
  for (el in c("panel.grid.major", "panel.grid.minor", "axis.ticks", "panel.background")) {
    expect_equal(th[[el]], tl[[el]], info = el)
  }
})

test_that("A faceted plot with theme_coriell builds and resolves elements", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  df <- data.frame(x = 1:6, y = c(2, 4, 1, 3, 5, 6), g = rep(c("a", "b"), 3))
  p <- ggplot2::ggplot(df, ggplot2::aes(x, y, colour = g)) +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(~g) +
    ggplot2::labs(title = "T", subtitle = "s") +
    theme_coriell()
  expect_equal(nrow(ggplot2::ggplot_build(p)$data[[1]]), 6)
  expect_no_error(ggplot2::ggplotGrob(p))
  th <- theme_coriell()
  expect_equal(ggplot2::calc_element("plot.title", th)$size, 22)
  expect_equal(ggplot2::calc_element("axis.text.x.bottom", th)$size, 16)
  expect_identical(ggplot2::calc_element("text", th)$family, "Helvetica")
})
