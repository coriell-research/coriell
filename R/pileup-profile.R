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
#' returns pileup stats relative to the 5' to 3' direction taking strandedness 
#' of the input range into account. 
#'
#' @param bamfile Path to an indexed BAM file
#' @param x GRanges to calculate pileups over
#' @param scan_bam_flag An instance of Rsamtools::scanBamFlag() used to to create 
#' a parameter object influencing what fields and which records are imported 
#' from a (binary) BAM file
#' @param max_depth maximum number of overlapping alignments considered 
#' for each position in the pileup. See Rsamtools::PileupParam()
#' @param min_base_quality minimum ‘QUAL’ value for each nucleotide 
#' in an alignment. Use phred2ASCIIOffset to help translate numeric or 
#' character values to these offsets.
#' @param min_mapq minimum ‘MAPQ’ value for an alignment to be included 
#' in pileup.
#' @param min_nucleotide_depth minimum count of each nucleotide 
#' (independent of other nucleotides) at a given position required for said 
#' nucleotide to appear in the result.
#' @param summarize Should the summary stats of the pileups at every relative
#' position over all ranges be returned? Default TRUE. If FALSE then a data.table
#' containing all measured ranges and relative positions is returned.
#' 
#' @details
#' For pileup computation, ranges are split by positive and negative strand 
#' and positions in the final data.table are returned relative to the 5' start
#' in the input ranges. If unstranded ranges are supplied then they are treated 
#' as positively stranded ranges. 
#' 
#' @return data.table
#' @import data.table
#' @export
#' 
#' @examples
#' \dontrun{
#' 
#' # Select GRanges to compute pileups over
#' gtf <- rtracklayer::import("/path/to/annotation.gtf")
#' genes <- gtf[gtf$type == "gene" & gtf$gene_biotype == "protein_coding", ]
#' 
#' # Define parameters for reading in BAM files
#' flags <- Rsamtools::scanBamFlag(
#'   isPaired = TRUE,
#'   isProperPair = TRUE,
#'   isUnmappedQuery = FALSE,
#'   hasUnmappedMate = FALSE
#'   )
#' 
#' # Compute pileups over all gene ranges  
#' result <- pileup_profile("/path/to/sorted.bam", x = genes, scan_bam_flag = flags)
#'
#' }
#' 
pileup_profile <- function(bamfile, x, scan_bam_flag, max_depth = 1e4, 
                           min_base_quality = 13, min_mapq = 1,
                           min_nucleotide_depth = 0, summarize = TRUE) {
  
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
  
  if (any(GenomicRanges::strand(x) == "*")) {
    message("Unstranded ranges detected. These will be treated as positively stranded.")
  }
  
  pos_gr <- x[GenomicRanges::strand(x) == "+" | GenomicRanges::strand(x) == "*"]
  neg_gr <- x[GenomicRanges::strand(x) == "-"]
  
  if (length(pos_gr) != 0) {
    pos_bparam <- Rsamtools::ScanBamParam(flag = scan_bam_flag, which = pos_gr)
    pp <- Rsamtools::pileup(bamfile, scanBamParam = pos_bparam, pileupParam = pparam)
    
    if (nrow(pp) > 0) {
      pp <- .pileup2dt(pp)
      pp[, relative_position := pos - start]
    } else {
      pp <- data.table::data.table()
    }
  } else {
    pp <- data.table::data.table()
  }
  
  if (length(neg_gr) != 0) {
    neg_bparam <- Rsamtools::ScanBamParam(flag = scan_bam_flag, which = neg_gr)
    np <- Rsamtools::pileup(bamfile, scanBamParam = neg_bparam, pileupParam = pparam)
    
    if (nrow(np > 0)) {
      np <- .pileup2dt(np)
      np[, relative_position := end - pos]
    } else {
      np <- data.table::data.table()
    }
  } else {
    np <- data.table::data.table()
  }
  
  dt <- rbind(pp, np)
  
  # Both pileups could be empty, if so return NULL
  if (nrow(dt) == 0) {
    return(NULL)
  }
  
  if (isFALSE(summarize)) {
    dt[, `:=`(chrom=NULL, start=NULL, end=NULL)]
    return(dt) 
  }
  
  result <- dt[, .(n = .N, total = sum(count), avg = mean(count)), by = relative_position]
  result[, normalized := coriell::minmax(avg)]
  data.table::setorder(result, relative_position)

  return(result)
}
