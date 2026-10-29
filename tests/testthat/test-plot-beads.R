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

test_that("plot_beads draws one bead per cell with fill from the matrix", {
  x <- matrix(c(1, 0, 0, 0, 1, 1), nrow = 2, byrow = TRUE)
  calls <- record_base_plot(plot_beads(x, m_col = "black", un_col = "white"))
  pts <- calls_named(calls, "C_plotXY")
  expect_length(pts, 1)
  xy <- pts[[1]]$args[[1]]
  bg <- pts[[1]]$args[[6]]
  got <- data.frame(read = xy$y, cpg = xy$x, bg = bg)
  got <- got[order(got$read, got$cpg), ]
  expect_equal(got$read, c(1, 1, 1, 2, 2, 2))
  expect_equal(got$cpg, c(1, 2, 3, 1, 2, 3))
  expect_equal(got$bg, c("black", "white", "white", "white", "black", "black"))
  expect_equal(sum(bg == "black"), sum(x == 1))
})

test_that("plot_beads draws one guide line per read and reverses the y-axis", {
  set.seed(1)
  x <- matrix(rbinom(20, 1, 0.5), nrow = 4)
  calls <- record_base_plot(plot_beads(x, dash_col = "grey50"))
  ab <- calls_named(calls, "C_abline")[[1]]$args
  expect_equal(ab[[3]], 1:4)
  expect_equal(ab[[6]], "grey50")
  win <- calls_named(calls, "C_plot_window")[[1]]$args
  expect_equal(win[[1]], c(1, 5))
  expect_equal(win[[2]], c(4, 1))
})

test_that("plot_beads handles single-read and single-CpG matrices", {
  xy <- calls_named(
    record_base_plot(plot_beads(matrix(c(1, 0, 1), nrow = 1))),
    "C_plotXY"
  )[[1]]$args[[1]]
  expect_equal(xy$x, 1:3)
  expect_equal(xy$y, rep(1, 3))

  xy <- calls_named(
    record_base_plot(plot_beads(matrix(c(1, 0, 1), ncol = 1))),
    "C_plotXY"
  )[[1]]$args[[1]]
  expect_equal(xy$x, rep(1, 3))
  expect_equal(xy$y, 1:3)
})

test_that("plot_beads applies custom colors and axis labels", {
  calls <- record_base_plot(
    plot_beads(diag(2), m_col = "red", un_col = "blue", x_lab = "Pos", y_lab = "Mol")
  )
  expect_equal(
    calls_named(calls, "C_plotXY")[[1]]$args[[6]],
    c("red", "blue", "blue", "red")
  )
  ttl <- calls_named(calls, "C_title")
  expect_equal(ttl[[1]]$args[[3]], "Pos")
  expect_equal(ttl[[2]]$args[[4]], "Mol")
})

test_that("plot_beads never shows NA CpGs as methylated", {
  x <- matrix(c(1, NA, 0, 1), 2)
  bg <- calls_named(record_base_plot(plot_beads(x)), "C_plotXY")[[1]]$args[[6]]
  expect_equal(bg, c("black", NA, "white", "black"))
})
