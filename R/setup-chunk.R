#' Setup Chunk Generator
#'
#' This function prints R code to the console that can be copied and pasted
#' into an analysis script or RMarkdown file. It sets up standard libraries
#' and directory paths.
#'
#' @param extra_libs Character vector. Additional library names to load.
#' @param dir_name Character string. The specific subdirectory for results (e.g., "01", "02").
#'
#' @return Prints formatted code to console. Returns the character vector invisibly.
#' @export
#' @examples
#' \dontrun{
#'  setup_chunk()
#' }
#'
setup_chunk <- function(extra_libs = NULL, dir_name = "01") {
  base_libs <- c(
    "here",
    "coriell",
    "data.table",
    "ggplot2",
    "patchwork",
    "edgeR",
    "SummarizedExperiment"
  )

  all_libs <- unique(c(base_libs, extra_libs))
  lib_lines <- paste0("library(", all_libs, ")")
  path_lines <- c(
    paste0('datafiles <- here("results", "data-files", "', dir_name, '")'),
    paste0('figures <- here("results", "figures", "', dir_name, '")'),
    paste0('rdsfiles <- here("results", "rds-files", "', dir_name, '")')
  )
  dir_creation_line <- "invisible(sapply(c(datafiles, figures, rdsfiles), dir.create, showWarnings=F, recursive=T))"

  final_output <- c(
    lib_lines,
    "",
    "",
    path_lines,
    "",
    dir_creation_line
  )

  cat(final_output, sep = "\n")

  invisible(final_output)
}
