#' Read in and Filter Bismark Coverage Files
#'
#' This function will read in Bismark coverage files and optionally filter for
#' coverage and feature-wise variance.
#' @details The function can optionally return a list of CpG x Sample matrices of
#' the filtered data if return_mats = TRUE. The following matrices will be
#' returned for each of the retained CpG sites: Coverage, Percent Methylation
#' (Percent), count of methylated CpGs (Methylated) and count of unmethylated
#' CpGs (Unmethylated)
#' @import data.table
#' @param files Named vector of file paths to the bismark coverage files
#' @param coverage Minimum coverage for a CpG site. Default (1)
#' @param prop_samples Proportion of samples that site must be covered to be retained. Default (1)
#' @param remove_zero_var Should zero variance features be removed after filtering for coverage? Default (FALSE)
#' @param return_mats if TRUE, return matrices of the filtered data in a list. Else, return the filtered data.table in bismark format. Default (FALSE)
#' @return Either a data.table of filtered data in bismark format with an added
#' Coverage column (Me + Un) and, if remove_zero_var = TRUE, a Var column with the
#' per-site variance of Percent; or, if return_mats = TRUE, a list of matrices of
#' the filtered data with columns in the order of \code{files}. Note that
#' when prop_samples < 1, sites are filtered as a whole: rows from samples below
#' the coverage threshold are kept at sites that pass.
#' @export
#' @examples
#' # Write three small bismark coverage files
#' # (chrom, start, end, percent, methylated, unmethylated)
#' cov <- list(
#'   S1 = "chr1\t100\t100\t50\t5\t5\nchr1\t200\t200\t100\t20\t0\nchr2\t50\t50\t0\t0\t10\n",
#'   S2 = "chr1\t100\t100\t80\t8\t2\nchr1\t200\t200\t100\t30\t0\nchr2\t50\t50\t25\t1\t3\n",
#'   S3 = "chr1\t100\t100\t20\t2\t8\nchr1\t200\t200\t100\t20\t0\nchr2\t50\t50\t40\t4\t6\n"
#' )
#' files <- vapply(names(cov), function(n) {
#'   f <- file.path(tempdir(), paste0(n, ".bismark.cov"))
#'   writeLines(cov[[n]], f, sep = "")
#'   f
#' }, character(1))
#'
#' # Read in all files using default parameters -- returns a data.table
#' result <- read_bismark(files)
#'
#' # Filter for coverage >= 10 in all samples and remove zero variance features
#' dt <- read_bismark(files, coverage = 10, remove_zero_var = TRUE)
#'
#' # Same filtering but return a list of CpG x Sample matrices
#' l <- read_bismark(files, coverage = 10, return_mats = TRUE)
#' l$Coverage
#' l$Percent
#'
#' unlink(files)
read_bismark <- function(
  files,
  coverage = 1,
  prop_samples = 1,
  remove_zero_var = FALSE,
  return_mats = FALSE
) {
  stopifnot("Coverage must be >= 0" = coverage >= 0)
  stopifnot(
    "prop_samples must be between 0 and 1" = prop_samples > 0 &
      prop_samples <= 1
  )
  if (is.null(names(files))) {
    warning(
      "files are not named. Samples in the output data.table will be labelled with their index in files."
    )
  }

  message("\nReading in files...")
  col_names <- c("Chrom", "Start", "End", "Percent", "Me", "Un")
  dt <- rbindlist(lapply(files, fread, col.names = col_names), idcol = "Sample")
  dt[, Coverage := Me + Un]

  # Filter for coverage
  message(paste0(
    "Filtering for Coverage >= ",
    coverage,
    " in ",
    prop_samples * 100,
    "% of samples..."
  ))
  # Integer number of samples, tolerant of floating point error
  keep_n <- ceiling(length(files) * prop_samples - sqrt(.Machine$double.eps))
  keep <- dt[Coverage >= coverage, .N, by = .(Chrom, Start)][
    N >= keep_n,
    .(Chrom, Start)
  ]
  dt <- dt[keep, on = .(Chrom, Start), nomatch = 0L]

  # Filter for zero variance features
  if (remove_zero_var) {
    message("Removing any zero variance features...")
    hasVar <- dt[, .(Var = var(Percent)), by = .(Chrom, Start)][Var != 0]
    dt <- dt[hasVar, on = .(Chrom, Start), nomatch = 0L]
  }

  result <- dt
  if (return_mats) {
    message("Coercing data to matrices...")
    # Keep samples in input file order (dcast orders columns by factor levels)
    sample_ids <- if (is.null(names(files))) seq_along(files) else names(files)
    dt[, Sample := factor(Sample, levels = unique(sample_ids))]
    Coverage <- dcast(
      dt,
      Chrom + Start ~ Sample,
      value.var = "Coverage",
      fill = 0L
    )
    Percent <- dcast(
      dt,
      Chrom + Start ~ Sample,
      value.var = "Percent",
      fill = NA
    )
    Methylated <- dcast(dt, Chrom + Start ~ Sample, value.var = "Me", fill = NA)
    Unmethylated <- dcast(
      dt,
      Chrom + Start ~ Sample,
      value.var = "Un",
      fill = NA
    )

    # Add unique row names and remove extraneous columns
    Coverage[, Location := paste(Chrom, Start, sep = ".")][,
      c("Chrom", "Start") := NULL
    ]
    Percent[, Location := paste(Chrom, Start, sep = ".")][,
      c("Chrom", "Start") := NULL
    ]
    Methylated[, Location := paste(Chrom, Start, sep = ".")][,
      c("Chrom", "Start") := NULL
    ]
    Unmethylated[, Location := paste(Chrom, Start, sep = ".")][,
      c("Chrom", "Start") := NULL
    ]

    # Convert to matrices
    Coverage <- as.matrix(Coverage, rownames = "Location")
    Percent <- as.matrix(Percent, rownames = "Location")
    Methylated <- as.matrix(Methylated, rownames = "Location")
    Unmethylated <- as.matrix(Unmethylated, rownames = "Location")

    result <- list(
      Coverage = Coverage,
      Percent = Percent,
      Methylated = Methylated,
      Unmethylated = Unmethylated
    )
  }
  message("Done.")

  return(result)
}
