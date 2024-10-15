#' Convert pileups to data.tables
#' 
#' @param p pileup result
#' 
#' @import data.table
.pileup2dt <- function(p) {
  
  dt <- as.data.table(p)
  dt[, c("chrom", "start_end") := data.table::tstrsplit(which_label, ":", fixed = TRUE)][,
       c("start", "end") := data.table::tstrsplit(start_end, "-", fixed = TRUE)][,
       `:=`(start = as.integer(start), end = as.integer(end), start_end = NULL)]
  
  return(dt)
}

#' Return the pileup of reads over a set of GRanges
#' 
#' This function computes pileups for BAM files over a set of GRanges. It 
#' returns sum and average counts for all positions relative to the 5' end of 
#' the ranges. This function also allows for subsampling of ranges to speed up
#' calculations. 
#'
#' @param bamfile Path to an indexed BAM file
#' @param x GRanges to calculate pileups over
#' @param scan_bam_flag An instance of Rsamtools::scanBamFlag() used to to create 
#' a parameter object influencing what fields and which records are imported 
#' from a (binary) BAM file
#' @param sample_size Number of records ranges to be randomly subsampled from x. 
#' Default NULL, no subsampling.
#' @param max_depth integer; maximum number of overlapping alignments considered 
#' for each position in the pileup. See Rsamtools::PileupParam()
#' @param min_base_quality integer; minimum ‘QUAL’ value for each nucleotide 
#' in an alignment. Use phred2ASCIIOffset to help translate numeric or 
#' character values to these offsets.
#' @param min_mapq integer; minimum ‘MAPQ’ value for an alignment to be included 
#' in pileup.
#' @param min_nucleotide_depth integer(1); minimum count of each nucleotide 
#' (independent of other nucleotides) at a given position required for said 
#' nucleotide to appear in the result.
#' 
#' @details
#' This function takes in a BAM file and computes a pileup over the given ranges. 
#' The data.table which gets returned contains the sum and average count of 
#' reads at each position relative to the 5' end of the input ranges. The user 
#' can also randomly subsample the input ranges to speed up computation. 
#' For both subsampling and pileup computation, ranges are split by positive and 
#' negative strand, subsampled at the same depth (if set), and positions in the 
#' final data.table are returned relative to the 5' start in the input ranges.
#' 
#' @return data.table
#' @import data.table
#' @export
pileup_profile <- function(bamfile, x, scan_bam_flag, sample_size = NULL,
                           max_depth = 1e4, min_base_quality = 13, min_mapq = 1,
                           min_nucleotide_depth = 0) {
  
  if (!requireNamespace("Rsamtools", quietly = TRUE)) {
    stop("Rsamtools is required")
  }
  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    stop("GenomicRanges is required")
  }
  
  pparam <- Rsamtools::PileupParam(
    max_depth = max_depth,
    min_base_quality = min_base_quality,
    min_mapq = min_mapq,
    min_nucleotide_depth = min_nucleotide_depth,
    distinguish_nucleotides = FALSE,
    distinguish_strands = FALSE
  )
  
  if (all(GenomicRanges::strand(x) == "*")) {
    stop("Unstranded ranges are not yet supported!")
  }

  pos_gr <- x[GenomicRanges::strand(x) == "+"]
  neg_gr <- x[GenomicRanges::strand(x) == "-"]

  if (!is.null(sample_size)) {
    n <- floor(sample_size / 2)
    if (n > length(pos_gr) || n > length(neg_gr)) {
      n <- min(c(length(pos_gr), length(neg_gr)))
      message("subsample size > length of ranges on either strand. Setting n to ", n)
    }
    pos_gr <- pos_gr[sample.int(length(pos_gr), size = n)]
    neg_gr <- neg_gr[sample.int(length(neg_gr), size = n)]
  }

  pos_bparam <- Rsamtools::ScanBamParam(flag = scan_bam_flag, which = pos_gr)
  neg_bparam <- Rsamtools::ScanBamParam(flag = scan_bam_flag, which = neg_gr)

  pp <- .pileup2dt(Rsamtools::pileup(bamfile, scanBamParam = pos_bparam, pileupParam = pparam))
  np <- .pileup2dt(Rsamtools::pileup(bamfile, scanBamParam = neg_bparam, pileupParam = pparam))
  
  # Start counting from the 5' end
  pp[, relative_position := pos - start]
  np[, relative_position := end - pos]
  df <- rbind(pp, np)

  result <- df[, .(n = .N, total = sum(count), avg = mean(count)), by = relative_position]
  result[, normalized := coriell::minmax(avg)]
  data.table::setorder(result, relative_position)

  return(result)
}
