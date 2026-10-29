# Replace the clipboard write so tests do not touch the system clipboard
mock_clip <- function(fun = function(content, ...) invisible(content), env = parent.frame()) {
  if (requireNamespace("clipr", quietly = TRUE)) {
    local_mocked_bindings(write_clip = fun, .package = "clipr", .env = env)
  }
}

base_libs <- c(
  "here", "coriell", "data.table", "ggplot2", "patchwork", "edgeR",
  "SummarizedExperiment"
)

test_that("Default chunk has the expected lines in order", {
  mock_clip()
  out <- suppressMessages(capture.output(setup_chunk()))
  expect_length(out, 14)
  expect_identical(out[1:7], paste0("suppressPackageStartupMessages(library(", base_libs, "))"))
  expect_identical(out[8:9], c("", ""))
  expect_identical(out[10:12], c(
    'datafiles <- here("results", "data-files", "01")',
    'figures <- here("results", "figures", "01")',
    'rdsfiles <- here("results", "rds-files", "01")'
  ))
  expect_identical(out[13], "")
  expect_identical(
    out[14],
    "invisible(sapply(c(datafiles, figures, rdsfiles), dir.create, showWarnings = F, recursive = T))"
  )
})

test_that("extra_libs are appended and de-duplicated; dir_name is used in all paths", {
  mock_clip()
  out <- suppressMessages(capture.output(
    setup_chunk(extra_libs = c("ggplot2", "limma", "limma"), dir_name = "07-dmr")
  ))
  libs <- sub(
    "^suppressPackageStartupMessages\\(library\\((.*)\\)\\)$", "\\1",
    grep("^suppress", out, value = TRUE)
  )
  expect_identical(libs, c(base_libs, "limma"))
  paths <- grep("<- here", out, value = TRUE)
  expect_length(paths, 3)
  expect_true(all(grepl('"07-dmr"\\)$', paths)))
})

test_that("Generated code parses and the path lines evaluate", {
  mock_clip()
  out <- suppressMessages(capture.output(setup_chunk()))
  exprs <- parse(text = out)
  expect_length(exprs, 11)
  env <- new.env()
  env$here <- function(...) file.path("ROOT", ...)
  for (e in exprs[8:10]) eval(e, env)
  expect_identical(env$datafiles, file.path("ROOT", "results", "data-files", "01"))
  expect_identical(env$figures, file.path("ROOT", "results", "figures", "01"))
  expect_identical(env$rdsfiles, file.path("ROOT", "results", "rds-files", "01"))
})

test_that("Chunk is returned invisibly, with or without quiet", {
  mock_clip()
  out <- suppressMessages(capture.output(v <- withVisible(setup_chunk())))
  expect_false(v$visible)
  expect_type(v$value, "character")
  expect_identical(v$value, out)

  out_q <- suppressMessages(capture.output(vq <- withVisible(setup_chunk(quiet = TRUE))))
  expect_identical(out_q, character(0))
  expect_false(vq$visible)
  expect_identical(vq$value, out)
})

test_that("Clipboard receives exactly the printed lines", {
  skip_if_not_installed("clipr")
  clip <- NULL
  mock_clip(function(content, ...) {
    clip <<- content
    invisible(content)
  })
  expect_message(out <- capture.output(setup_chunk(extra_libs = "limma")), "copied to clipboard")
  expect_identical(clip, out)
  clip <- NULL
  expect_message(out2 <- capture.output(setup_chunk(quiet = TRUE)), "copied")
  expect_identical(out2, character(0))
  expect_length(clip, 14)
})

test_that("Clipboard failure is a warning, and the chunk is still printed", {
  skip_if_not_installed("clipr")
  mock_clip(function(...) stop("nope"))
  expect_warning(out <- capture.output(setup_chunk()), "Could not copy to clipboard")
  expect_length(out, 14)
})
