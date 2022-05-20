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
#' @return Either a data.table of filtered data with additional columns for Coverage and Variance or a list of matrices of filtered data
#' @export
#' @examples
#' \dontrun{
#' files <- c("P6_1.bismark.cov.gz", "P6_4.bismark.cov.gz", "P7_2.bismark.cov.gz", "P7_5.bismark.cov.gz", "P8_3.bismark.cov.gz", "P8_6.bismark.cov.gz")
#' names(files) <- gsub("\\.bismark\\.cov\\.gz", "", files)
#' 
#' # Read in all files using default parameters -- returns a data.table
#' result <- read_bismark(files)
#' 
#' # Filter for coverage >= 20 in at least 50% of samples and remove any zero variance features
#' dt <- read_bismark(files, coverage = 20, prop_samples = 0.5, remove_zero_var = TRUE)
#' 
#' # Same filtering as above but return a list of CpG x Sample matrices
#' l <- read_bismark(files, coverage = 20, prop_samples = 0.5, remove_zero_var = TRUE, return_mats = TRUE)
#' 
#' # To view the CpG x Sample Coverages:
#' l$Coverage
#' 
#' # Percent methylation
#' l$Percent
#' }
read_bismark <- function(files, coverage = 1, prop_samples = 1, remove_zero_var = FALSE, return_mats = FALSE) {
  stopifnot("Coverage must be >= 0" = coverage >= 0)
  stopifnot("prop_samples must be between 0 and 1" = prop_samples > 0 & prop_samples <= 1)
  if (is.null(names(files))) {
    warning("files are not named. Samples in the output data.table will be labelled with their index in files.")
  }
  
  message("\nReading in files...")
  col_names <- c("Chrom", "Start", "End", "Percent", "Me", "Un")
  dt <- rbindlist(lapply(files, fread, col.names = col_names), idcol = "Sample")
  dt[, Coverage := Me + Un]
  
  # Filter for coverage
  message(paste0("Filtering for Coverage >= ", coverage, " in ", prop_samples * 100, "% of samples..."))
  keep_n <- length(files) * prop_samples
  keep <- dt[Coverage >= coverage, .N, by = .(Chrom, Start)][N >= keep_n, .(Chrom, Start)]
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
    Coverage <- dcast(dt, Chrom + Start ~ Sample, value.var = "Coverage", fill = 0L)
    Percent <- dcast(dt, Chrom + Start ~ Sample, value.var = "Percent", fill = NA)
    Methylated <- dcast(dt, Chrom + Start ~ Sample, value.var = "Me", fill = NA)
    Unmethylated <- dcast(dt, Chrom + Start ~ Sample, value.var = "Un", fill = NA)
    
    # Add unique row names and remove extraneous columns
    Coverage[, Location := paste(Chrom, Start, sep = ".")][, c("Chrom", "Start") := NULL]
    Percent[, Location := paste(Chrom, Start, sep = ".")][, c("Chrom", "Start") := NULL]
    Methylated[, Location := paste(Chrom, Start, sep = ".")][, c("Chrom", "Start") := NULL]
    Unmethylated[, Location := paste(Chrom, Start, sep = ".")][, c("Chrom", "Start") := NULL]
    
    # Convert to matrices
    Coverage <- as.matrix(Coverage, rownames = "Location")
    Percent <- as.matrix(Percent, rownames = "Location")
    Methylated <- as.matrix(Methylated, rownames = "Location")
    Unmethylated <- as.matrix(Unmethylated, rownames = "Location")
    
    result <- list(Coverage = Coverage, Percent = Percent, Methylated = Methylated, Unmethylated = Unmethylated)
  }
  message("Done.")
  
  return(result)
}
