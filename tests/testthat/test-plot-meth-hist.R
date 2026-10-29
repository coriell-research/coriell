# Record the base-graphics display list as (name, args) pairs
record_base_plot <- function(expr) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  grDevices::dev.control("enable")
  force(expr)
  dl <- grDevices::recordPlot()[[1]]
  lapply(dl, function(e) {
    f <- e[[2]][[1]]
    list(name = if (is.list(f)) f$name else NA_character_, args = e[[2]][-1])
  })
}
calls_named <- function(calls, nm) Filter(function(e) identical(e$name, nm), calls)

bin_labels <- paste(seq(0, 0.9, 0.1), seq(0.1, 1, 0.1), sep = " - ")

test_that("plot_meth_hist bar heights and labels match hand counts", {
  x <- c(rep(0.05, 30), rep(0.35, 20), rep(0.95, 50))
  calls <- record_base_plot(plot_meth_hist(x))
  counts <- c(30, 0, 0, 20, 0, 0, 0, 0, 0, 50)
  expect_equal(calls_named(calls, "C_rect")[[1]]$args[[4]], counts)
  txt <- calls_named(calls, "C_text")[[1]]$args
  expect_equal(txt[[2]], paste0(counts, "%"))
  expect_equal(txt[[1]]$y, counts)
  expect_equal(calls_named(calls, "C_axis")[[1]]$args[[3]], bin_labels)
})

test_that("plot_meth_hist drops NAs and matches table(cut()) on fixed 0.1 bins", {
  set.seed(10)
  x <- c(rbeta(997, 0.5, 0.5), NA, NA, NA)
  calls <- record_base_plot(plot_meth_hist(x))
  expected <- as.vector(table(cut(
    x[!is.na(x)],
    breaks = seq(0, 1, 0.1),
    include.lowest = TRUE
  )))
  expect_equal(calls_named(calls, "C_rect")[[1]]$args[[4]], expected)
  lab <- calls_named(calls, "C_text")[[1]]$args[[2]]
  expect_equal(as.numeric(sub("%", "", lab)), round(100 * expected / 997))
})

test_that("plot_meth_hist bins are right-closed and include 0", {
  r <- calls_named(
    record_base_plot(plot_meth_hist(c(0, 0.1, 0.1, 0.2, 1))),
    "C_rect"
  )[[1]]$args
  expect_equal(r[[4]], c(3, 1, 0, 0, 0, 0, 0, 0, 0, 1))
})

test_that("plot_meth_hist uses fixed [0, 1] bins for narrow data ranges", {
  x <- seq(0.4, 0.6, length.out = 1000)
  calls <- record_base_plot(plot_meth_hist(x))
  expect_equal(calls_named(calls, "C_axis")[[1]]$args[[3]], bin_labels)
  expected <- as.vector(table(cut(x, seq(0, 1, 0.1), include.lowest = TRUE)))
  expect_equal(calls_named(calls, "C_rect")[[1]]$args[[4]], expected)
})

test_that("plot_meth_hist applies colors, headroom and axis label", {
  x <- c(rep(0.05, 30), rep(0.95, 70))
  calls <- record_base_plot(
    plot_meth_hist(x, bar_color = "red", border_color = "blue", plot_title = "T")
  )
  expect_equal(calls_named(calls, "C_plot_window")[[1]]$args[[2]], c(0, 84))
  r <- calls_named(calls, "C_rect")[[1]]$args
  expect_equal(r$col, "red")
  expect_equal(r$border, "blue")
  expect_equal(calls_named(calls, "C_mtext")[[1]]$args[[1]], "Methylation Value")
})
