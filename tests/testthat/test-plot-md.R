md <- data.frame(
  gene = paste0("g", 1:8),
  logCPM = 1:8,
  logFC = c(2, -2, 1, -1, 3, -3, 0, 1.5),
  FDR = c(0.01, 0.01, 0.01, 0.01, 0.05, 0.2, 0.01, NA)
)

# Run base-graphics code on a null device and return its value and display list
record_base_plot <- function(code) {
  grDevices::pdf(NULL)
  grDevices::dev.control("enable")
  on.exit(grDevices::dev.off(), add = TRUE)
  res <- force(code)
  list(res = res, dl = grDevices::recordPlot()[[1]])
}

# Arguments of every call to the named graphics primitive in a display list
dl_args <- function(dl, name) {
  calls <- Filter(function(e) identical(e[[2]][[1]]$name, name), dl)
  lapply(calls, function(e) e[[2]][-1])
}

# plot_md ---------------------------------------------------------------------

test_that("plot_md classifies with strict cutoffs and maps x/y, lines and labels", {
  p <- plot_md(md, fdr = 0.05, lfc = 1)
  expect_s3_class(p, "ggplot")
  up <- ggplot2::layer_data(p, 3)
  dn <- ggplot2::layer_data(p, 2)
  un <- ggplot2::layer_data(p, 1)
  expect_equal(c(up$x, up$y), c(1, 2))
  expect_equal(c(dn$x, dn$y), c(2, -2))
  expect_equal(sort(un$x), 3:8)
  hl <- sapply(4:6, function(i) ggplot2::layer_data(p, i)$yintercept)
  expect_equal(hl, c(0, 1, -1))
  expect_equal(p$labels$x, "Average logCPM")
  expect_equal(p$labels$caption, "FDR = 0.05\nlfc cutoff = 1")
})

test_that("plot_md count labels are placed by offsets times axis limits", {
  p <- plot_md(md, fdr = 0.05, lfc = 1)
  up <- ggplot2::layer_data(p, 7)
  dn <- ggplot2::layer_data(p, 8)
  expect_equal(c(up$x, up$y), c(0.85 * 8, 0.8 * 3))
  expect_equal(c(dn$x, dn$y), c(0.85 * 8, 0.8 * -3))
  p2 <- plot_md(md, fdr = 0.05, lfc = 1, x_axis_limits = c(0, 10), y_axis_limits = c(-5, 5))
  expect_equal(ggplot2::layer_data(p2, 7)[, c("x", "y")], data.frame(x = 8.5, y = 4))
  expect_equal(ggplot2::layer_data(p2, 8)[, c("x", "y")], data.frame(x = 8.5, y = -4))
})

test_that("plot_md layer counts and labels match direct counts on real data", {
  d <- GSE161650_de
  p <- plot_md(d, fdr = 0.01, lfc = 1)
  n <- sapply(1:3, function(i) nrow(ggplot2::layer_data(p, i)))
  n_up <- sum(d$FDR < 0.01 & d$logFC > 1)
  n_down <- sum(d$FDR < 0.01 & d$logFC < -1)
  expect_equal(n, c(nrow(d) - n_up - n_down, n_down, n_up))
  expect_equal(ggplot2::layer_data(p, 7)$label, paste0(n_up, "\n", round(100 * n_up / nrow(d), 1), "%"))
  expect_equal(ggplot2::layer_data(p, 8)$label, paste0(n_down, "\n", round(100 * n_down / nrow(d), 1), "%"))
})

test_that("plot_md respects custom column names and label_sig labels only DE features", {
  skip_if_not_installed("ggrepel")
  d <- setNames(md, c("id", "AveExpr", "lfc_col", "padj"))
  p <- suppressWarnings(plot_md(
    d, x = "AveExpr", y = "lfc_col", sig_col = "padj", lab = "id",
    label_sig = TRUE, fdr = 0.05, lfc = 1, annotate_counts = FALSE
  ))
  expect_setequal(ggplot2::layer_data(p, 7)$label, c("g1", "g2"))
  expect_match(p$labels$caption, "^padj = 0.05")
})

test_that("plot_md does not modify a data.table input by reference", {
  dt <- data.table::as.data.table(md)
  nm <- names(dt)
  plot_md(dt)
  expect_identical(names(dt), nm)
})

test_that("plot_md percent labels are rounded once to lab_digits", {
  # 6 / 13 = 46.1538...; rounding to 2 then 1 digit would give 46.1
  d <- data.frame(logFC = c(rep(1, 6), rep(0, 7)), FDR = c(rep(0.01, 6), rep(0.5, 7)), logCPM = 1:13)
  expect_equal(ggplot2::layer_data(plot_md(d), 7)$label, "6\n46.2%")
  expect_equal(ggplot2::layer_data(plot_md(d, lab_digits = 3), 7)$label, "6\n46.154%")
})

test_that("plot_md percent denominator excludes NA rows, matching plot_md2", {
  # md has 8 rows, one with NA FDR: 1 up and 1 down out of 7 complete rows
  p <- plot_md(md, fdr = 0.05, lfc = 1)
  expect_equal(ggplot2::layer_data(p, 7)$label, "1\n14.3%")
  expect_equal(ggplot2::layer_data(p, 8)$label, "1\n14.3%")
  txt <- dl_args(record_base_plot(plot_md2(md, fdr = 0.05, lfc = 1))$dl, "C_text")
  expect_equal(c(txt[[1]][[2]], txt[[2]][[2]]), c("1\n14.3%", "1\n14.3%"))
})

# plot_md2 --------------------------------------------------------------------

test_that("plot_md2 classifies, filters NAs and sets default limits and lines", {
  rec <- record_base_plot(plot_md2(md, fdr = 0.05, lfc = 1))
  r <- rec$res
  expect_equal(nrow(r$data), 7)
  expect_equal(unname(r$counts), c(1, 1, 5))
  expect_equal(r$data$group[1:2], c("Up", "Down"))
  w <- dl_args(rec$dl, "C_plot_window")[[1]]
  expect_equal(w[[1]], c(1, 7))
  expect_equal(w[[2]], c(-3, 3))
  ab <- dl_args(rec$dl, "C_abline")
  expect_equal(ab[[1]][[3]], 0)
  expect_equal(ab[[2]][[3]], c(-1, 1))
  expect_equal(dl_args(rec$dl, "C_mtext")[[1]][[1]], "FDR < 0.05; lfc cutoff = 1.00")
})

test_that("plot_md2 drops non-finite values and FDR > 1", {
  d <- data.frame(
    logCPM = c(1, Inf, 3, 4),
    logFC = c(1, 1, -Inf, -1),
    FDR = c(0.01, 0.01, 0.01, 1.2)
  )
  r <- record_base_plot(plot_md2(d))$res
  expect_equal(nrow(r$data), 1)
  expect_equal(unname(r$counts), c(1, 0, 0))
})

test_that("plot_md2 y range includes 0 and degenerate ranges are widened", {
  rec <- record_base_plot(plot_md2(data.frame(logCPM = 1:3, logFC = c(1, 2, 3), FDR = 0.01)))
  expect_equal(dl_args(rec$dl, "C_plot_window")[[1]][[2]], c(0, 3))
  rec1 <- record_base_plot(plot_md2(data.frame(logCPM = 5, logFC = 0, FDR = 0.5)))
  w <- dl_args(rec1$dl, "C_plot_window")[[1]]
  expect_equal(w[[1]], c(4, 6))
  expect_equal(w[[2]], c(-1, 1))
})

test_that("plot_md2 count label positions and text", {
  rec <- record_base_plot(plot_md2(md, fdr = 0.05, lfc = 1))
  txt <- dl_args(rec$dl, "C_text")
  # par("usr") extends the limits (1, 7) and (-3, 3) by 4% on each side
  usr <- c(1 - 0.24, 7 + 0.24, -3 - 0.24, 3 + 0.24)
  expect_equal(c(txt[[1]][[1]]$x, txt[[1]][[1]]$y), c(0.85 * usr[2], 0.8 * usr[4]))
  expect_equal(c(txt[[2]][[1]]$x, txt[[2]][[1]]$y), c(0.85 * usr[2], 0.8 * usr[3]))
  expect_equal(c(txt[[1]][[2]], txt[[2]][[2]]), c("1\n14.3%", "1\n14.3%"))
})

test_that("plot_md2 label_sig selects the top n_labels by combined rank", {
  set.seed(2)
  k <- 50
  d <- data.frame(
    gene = paste0("g", 1:k),
    logCPM = runif(k, 0, 10),
    logFC = rnorm(k, 0, 2),
    FDR = runif(k, 0, 0.2)
  )
  rec <- record_base_plot(plot_md2(
    d, lab = "gene", label_sig = TRUE, n_labels = 4, fdr = 0.1, lfc = 1,
    annotate_counts = FALSE
  ))
  s <- d[d$FDR < 0.1 & abs(d$logFC) > 1, ]
  score <- rank(s$FDR) + rank(-abs(s$logFC))
  expect_equal(dl_args(rec$dl, "C_text")[[1]][[2]], s$gene[order(score)][1:4])
})

test_that("plot_md2 counts match direct counts on real data", {
  d <- GSE161650_de
  r <- record_base_plot(plot_md2(d, fdr = 0.01, lfc = 1))$res
  n_up <- sum(d$FDR < 0.01 & d$logFC > 1)
  n_down <- sum(d$FDR < 0.01 & d$logFC < -1)
  expect_equal(unname(r$counts), c(n_up, n_down, nrow(d) - n_up - n_down))
})
