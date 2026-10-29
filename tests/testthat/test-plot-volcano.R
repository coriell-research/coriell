toy <- data.frame(
  gene = paste0("g", 1:8),
  logCPM = 1:8,
  logFC = c(2, -2, 0.5, -0.5, 3, -3, 0, 1),
  PValue = c(0.001, 0.002, 0.003, 0.004, 0.1, 0.15, 0.005, 0.5),
  FDR = c(0.01, 0.01, 0.01, 0.01, 0.2, 0.2, 0.01, NA)
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

make_sim <- function() {
  set.seed(1)
  n <- 200
  pv <- c(runif(20, 0, 1e-4), runif(180))
  data.frame(
    gene = paste0("g", 1:n),
    logFC = c(rnorm(20, 0, 3), rnorm(180)),
    PValue = pv,
    FDR = p.adjust(pv, "BH")
  )
}

# plot_volcano ----------------------------------------------------------------

test_that("plot_volcano puts each point in the right layer using strict thresholds", {
  p <- plot_volcano(toy, fdr = 0.05, lfc = 1)
  expect_s3_class(p, "ggplot")
  ld <- lapply(1:3, function(i) ggplot2::layer_data(p, i))
  expect_equal(ld[[3]]$x, 2)
  expect_equal(unique(ld[[3]]$colour), "red2")
  expect_equal(ld[[2]]$x, -2)
  expect_equal(unique(ld[[2]]$colour), "royalblue2")
  expect_setequal(ld[[1]]$x, c(0.5, -0.5, 3, -3, 0, 1))
  expect_equal(unique(ld[[1]]$colour), "grey40")
  expect_equal(sum(sapply(ld, nrow)), nrow(toy))
})

test_that("plot_volcano y is -log10 of the chosen column", {
  p <- plot_volcano(toy, fdr = 0.05, lfc = 1, annotate_counts = FALSE)
  expect_equal(ggplot2::layer_data(p, 3)$y, -log10(0.01))
  p2 <- plot_volcano(toy, y = "PValue", fdr = 0.05, lfc = 1, annotate_counts = FALSE)
  expect_equal(ggplot2::layer_data(p2, 3)$y, -log10(0.001))
  expect_equal(p2$labels$y, "-log10(PValue)")
})

test_that("plot_volcano values exactly at the cutoffs are not DE", {
  d <- data.frame(logFC = c(1, -1, 2, -2), FDR = c(0.01, 0.01, 0.05, 0.05))
  p <- plot_volcano(d, fdr = 0.05, lfc = 1, annotate_counts = FALSE)
  expect_equal(nrow(ggplot2::layer_data(p, 1)), 4)
  expect_equal(nrow(ggplot2::layer_data(p, 2)), 0)
  expect_equal(nrow(ggplot2::layer_data(p, 3)), 0)
})

test_that("plot_volcano reference lines are at 0, +/-lfc and -log10(fdr)", {
  p <- plot_volcano(toy, fdr = 0.05, lfc = 1, annotate_counts = FALSE)
  expect_equal(ggplot2::layer_data(p, 4)$xintercept, 0)
  expect_equal(ggplot2::layer_data(p, 5)$xintercept, 1)
  expect_equal(ggplot2::layer_data(p, 6)$xintercept, -1)
  expect_equal(ggplot2::layer_data(p, 7)$yintercept, -log10(0.05))
})

test_that("plot_volcano layer counts and labels match counts on real data", {
  d <- GSE161650_de
  p <- plot_volcano(d, fdr = 0.01, lfc = 1)
  n <- sapply(1:3, function(i) nrow(ggplot2::layer_data(p, i)))
  n_up <- sum(d$FDR < 0.01 & d$logFC > 1)
  n_down <- sum(d$FDR < 0.01 & d$logFC < -1)
  expect_equal(n, c(nrow(d) - n_up - n_down, n_down, n_up))
  labs <- c(ggplot2::layer_data(p, 8)$label, ggplot2::layer_data(p, 9)$label)
  expect_equal(labs, c(
    paste0(n_down, "\n", round(100 * n_down / nrow(d), 1), "%"),
    paste0(n_up, "\n", round(100 * n_up / nrow(d), 1), "%")
  ))
})

test_that("plot_volcano count labels are placed by offsets times axis limits", {
  p <- plot_volcano(toy, fdr = 0.05, lfc = 1)
  down <- ggplot2::layer_data(p, 8)
  up <- ggplot2::layer_data(p, 9)
  expect_equal(c(down$x, down$y), c(0.8 * -3, 0.9 * 2))
  expect_equal(c(up$x, up$y), c(0.8 * 3, 0.9 * 2))
  pl <- plot_volcano(toy, fdr = 0.05, lfc = 1, x_axis_limits = c(-4, 4), y_axis_limits = c(0, 3))
  expect_equal(ggplot2::layer_data(pl, 8)$x, -3.2)
  expect_equal(ggplot2::layer_data(pl, 9)[, c("x", "y")], data.frame(x = 3.2, y = 2.7))
})

test_that("plot_volcano annotates an empty direction as 0", {
  d <- data.frame(logFC = c(-2, -1, 0.1), FDR = c(0.01, 0.01, 0.5))
  p <- plot_volcano(d)
  expect_equal(ggplot2::layer_data(p, 9)$label, "0\n0%")
  expect_equal(ggplot2::layer_data(p, 8)$label, "2\n66.7%")
})

test_that("plot_volcano percent labels are rounded once to lab_digits", {
  # 6 / 13 = 46.1538...; rounding to 2 then 1 digit would give 46.1
  d <- data.frame(logFC = c(rep(1, 6), rep(0, 7)), FDR = c(rep(0.01, 6), rep(0.5, 7)))
  expect_equal(ggplot2::layer_data(plot_volcano(d), 9)$label, "6\n46.2%")
  expect_equal(
    ggplot2::layer_data(plot_volcano(d, lab_digits = 3), 9)$label,
    "6\n46.154%"
  )
})

test_that("plot_volcano percent denominator excludes NA rows, matching plot_volcano2", {
  # toy has 8 rows, one with NA FDR: 1 up and 1 down out of 7 complete rows
  p <- plot_volcano(toy, fdr = 0.05, lfc = 1)
  expect_equal(ggplot2::layer_data(p, 8)$label, "1\n14.3%")
  expect_equal(ggplot2::layer_data(p, 9)$label, "1\n14.3%")
  txt <- dl_args(record_base_plot(plot_volcano2(toy, fdr = 0.05, lfc = 1))$dl, "C_text")
  expect_equal(c(txt[[1]][[2]], txt[[2]][[2]]), c("1\n14.3%", "1\n14.3%"))
})

test_that("plot_volcano caption has no extra spaces", {
  p <- plot_volcano(data.frame(logFC = c(1, -1), FDR = c(0.01, 0.5)), lfc = 1)
  expect_equal(p$labels$caption, "FDR = 0.05\nlfc cutoff = 1")
})

test_that("plot_volcano annotate_counts and label_sig control extra layers", {
  skip_if_not_installed("ggrepel")
  expect_length(plot_volcano(toy, annotate_counts = FALSE)$layers, 7)
  p <- plot_volcano(toy, fdr = 0.05, lfc = 1, lab = "gene", label_sig = TRUE, annotate_counts = FALSE)
  expect_s3_class(p$layers[[8]]$geom, "GeomTextRepel")
  expect_setequal(ggplot2::layer_data(p, 8)$label, c("g1", "g2"))
  expect_message(plot_volcano(toy, label_sig = TRUE), "lab = NULL")
})

test_that("plot_volcano does not modify a data.table input by reference", {
  dt <- data.table::as.data.table(toy)
  nm <- names(dt)
  plot_volcano(dt)
  expect_identical(names(dt), nm)
})

# plot_volcano2 ---------------------------------------------------------------

test_that("plot_volcano2 group counts and -log10 p match direct computation", {
  sim <- make_sim()
  r <- record_base_plot(plot_volcano2(sim, fdr = 0.1, lfc = 0.5))$res
  sig <- sim$FDR < 0.1
  expect_equal(unname(r$counts), c(
    sum(sig & sim$logFC > 0.5),
    sum(sig & sim$logFC < -0.5),
    sum(!(sig & abs(sim$logFC) > 0.5))
  ))
  expect_equal(sum(r$counts), nrow(sim))
  expect_equal(r$data$neg_log10_p, -log10(sim$PValue))
})

test_that("plot_volcano2 p-value cutoff is the largest raw p among significant features", {
  sim <- make_sim()
  rec <- record_base_plot(plot_volcano2(sim, fdr = 0.1, lfc = 0.5))
  r <- rec$res
  expect_equal(r$raw_p_cutoff, max(sim$PValue[sim$FDR < 0.1]))
  expect_equal(r$y_cutoff, -log10(r$raw_p_cutoff))
  ab <- dl_args(rec$dl, "C_abline")
  expect_equal(ab[[1]][[4]], 0)
  expect_equal(ab[[2]][[4]], c(-0.5, 0.5))
  expect_equal(ab[[3]][[3]], r$y_cutoff)
  r2 <- record_base_plot(plot_volcano2(sim, fdr = 0.1, lfc = 100))$res
  expect_equal(r2$raw_p_cutoff, r$raw_p_cutoff)
  expect_equal(sum(r2$counts[c("Up", "Down")]), 0)
})

test_that("plot_volcano2 with no significant features has NA cutoffs", {
  sim <- make_sim()
  rec <- record_base_plot(plot_volcano2(transform(sim, FDR = 1)))
  expect_true(is.na(rec$res$raw_p_cutoff))
  expect_true(is.na(rec$res$y_cutoff))
  expect_length(dl_args(rec$dl, "C_abline"), 1)
  expect_match(dl_args(rec$dl, "C_mtext")[[1]][[1]], "^No features with FDR < 0.1")
})

test_that("plot_volcano2 drops invalid rows and floors p = 0", {
  bad <- data.frame(
    logFC = c(1, NA, 2, -1, 1),
    PValue = c(0.01, 0.01, 1.5, 0.01, -0.1),
    FDR = c(0.02, 0.02, 0.02, NA, 0.02)
  )
  r <- record_base_plot(plot_volcano2(bad))$res
  expect_equal(nrow(r$data), 1)
  expect_equal(unname(r$counts), c(1, 0, 0))
  z <- record_base_plot(plot_volcano2(
    data.frame(logFC = c(1, -1), PValue = c(0, 0.5), FDR = c(0, 0.5))
  ))$res
  expect_true(all(is.finite(z$data$neg_log10_p)))
  expect_equal(z$data$neg_log10_p[1], -log10(.Machine$double.xmin))
  expect_error(
    record_base_plot(plot_volcano2(data.frame(logFC = NA_real_, PValue = 0.1, FDR = 0.1))),
    "No valid observations"
  )
})

test_that("plot_volcano2 count labels use the filtered denominator and par('usr')", {
  rec <- record_base_plot(plot_volcano2(toy, fdr = 0.05, lfc = 1))
  txt <- dl_args(rec$dl, "C_text")
  expect_equal(txt[[1]][[2]], "1\n14.3%")
  expect_equal(txt[[2]][[2]], "1\n14.3%")
  # default axes: x range (-3, 3) and y up to -log10(min p) = 3.15, extended by 4%
  usr_x <- c(-3, 3) + c(-1, 1) * 0.04 * 6
  usr_y4 <- -log10(0.001) * 1.05 * 1.04
  expect_equal(txt[[1]][[1]]$x, 0.8 * usr_x[1])
  expect_equal(txt[[2]][[1]]$x, 0.8 * usr_x[2])
  expect_equal(txt[[1]][[1]]$y, 0.9 * usr_y4)
})

test_that("plot_volcano2 labels the top n_labels features by combined rank", {
  sim <- make_sim()
  rec <- record_base_plot(plot_volcano2(
    sim, fdr = 0.1, lfc = 0.5, lab = "gene", label_sig = TRUE,
    n_labels = 3, annotate_counts = FALSE
  ))
  sig <- sim[sim$FDR < 0.1 & abs(sim$logFC) > 0.5, ]
  score <- rank(sig$PValue) + rank(-abs(sig$logFC))
  expect_equal(dl_args(rec$dl, "C_text")[[1]][[2]], sig$gene[order(score)][1:3])
  rec0 <- record_base_plot(plot_volcano2(
    sim, lab = "gene", label_sig = TRUE, n_labels = 0, annotate_counts = FALSE
  ))
  expect_length(dl_args(rec0$dl, "C_text"), 0)
})

test_that("plot_volcano2 forwards ... and uses user axis limits verbatim", {
  sim <- make_sim()
  rec <- record_base_plot(plot_volcano2(
    sim, xlab = "custom", main = "T", x_axis_limits = c(-5, 5), y_axis_limits = c(0, 10)
  ))
  ti <- dl_args(rec$dl, "C_title")[[1]]
  expect_equal(ti[[1]], "T")
  expect_equal(ti[[3]], "custom")
  w <- dl_args(rec$dl, "C_plot_window")[[1]]
  expect_equal(w[[1]], c(-5, 5))
  expect_equal(w[[2]], c(0, 10))
})

test_that("plot_volcano2 counts match direct counts on real data", {
  d <- GSE161650_de
  r <- record_base_plot(plot_volcano2(d, fdr = 0.01, lfc = 1))$res
  n_up <- sum(d$FDR < 0.01 & d$logFC > 1)
  n_down <- sum(d$FDR < 0.01 & d$logFC < -1)
  expect_equal(unname(r$counts), c(n_up, n_down, nrow(d) - n_up - n_down))
})
