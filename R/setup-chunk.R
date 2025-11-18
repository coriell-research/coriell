#' Setup Chunk Generator
#'
#' This function prints R code to the console that can be copied and pasted into an analysis script
#' or quarto file. It sets up standard libraries and directory paths. If installed, \code{clipr} is
#' used to copy the code chuck to the clipboard.
#'
#' @param extra_libs Character vector. Additional library names to load.
#' @param dir_name Character string. The specific subdirectory for results (e.g., "01", "02").
#' @param quiet boolean. If TRUE do not print the code chunk to the console (only attempt to copy to clipboard).
#'
#' @return Prints formatted code to console. Returns the character vector invisibly.
#' @export
#' @examples
#' \dontrun{
#'  setup_chunk()
#' }
#'
setup_chunk <- function(extra_libs = NULL, dir_name = "01", quiet = FALSE) {
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
  lib_lines <- paste0("suppressPackageStartupMessages(library(", all_libs, "))")
  path_lines <- c(
    paste0('datafiles <- here("results", "data-files", "', dir_name, '")'),
    paste0('figures <- here("results", "figures", "', dir_name, '")'),
    paste0('rdsfiles <- here("results", "rds-files", "', dir_name, '")')
  )
  dir_creation_line <- "invisible(sapply(c(datafiles, figures, rdsfiles), dir.create, showWarnings = F, recursive = T))"

  final_output <- c(
    lib_lines,
    "",
    "",
    path_lines,
    "",
    dir_creation_line
  )

  if (requireNamespace("clipr", quietly = TRUE)) {
    tryCatch(
      {
        clipr::write_clip(final_output)
        message("\n[Code copied to clipboard]")
      },
      error = function(e) {
        warning(
          "Could not copy to clipboard: ",
          e$message,
          call. = FALSE
        )
      }
    )
  }

  if (isTRUE(quiet)) {
    return(invisible(NULL))
  }

  cat(final_output, sep = "\n")
}
