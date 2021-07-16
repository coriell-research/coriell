#' Process a quant file from REdiscoverTE pipeline
#'
#' Specialty function for processing quants files produced by the \href{https://www.nature.com/articles/s41467-019-13035-2}{REdiscoverTE}
#' pipeline. The function takes a quant file as input and calculates the counts
#' per gene and per repetitive element subfamily counts in a list. The user must 
#' provide the pre-processed gene and RepeatMasker annotations as data.tables.
#'
#' @param quant_file path to the quant file to be processed
#' @param gene_annot REdiscoverTE gene annotation data.frame
#' @param rmsk_annot REdiscoverTE RepeatMasker annotation data.frame
#' @return list containing five data.tables. One for the gene counts ("gene"), 
#' exonic RE counts ("exon"), intron RE counts ("intron"), intergenic RE counts ("intergenic"),
#' and intronic + intergenic RE counts ("re"). 
#' @import data.table
#' @export
#' @examples
#' \dontrun{
#' # read in annotation files
#' genes <- readRDS("path/to/gene-annotation.rds")
#' exonREs <- readRDS("path/to/exonRE-annotation.rds")
#' intronREs <- readRDS("path/to/intronRE-annotation.rds")
#' intergenicREs <- readRDS("path/to/intergenicRE-annotation.rds")
#'
#' # specify path to quant file
#' my_quant_file <- "/path/to/quant.sf"
#'
#' result <- process_quant_file(my_quant_file, genes, exonREs, intronREs, intergenicREs)
#'
#' # If you have a list of quant files, you can apply over all with multiple cores
#' my_quants <- list.files("path/to/quants_dir", pattern = "quant.sf", full.names = TRUE)
#'
#' processed <- parallel::mclapply(
#'   X = my_quants,
#'   FUN = process_quant_file,
#'   gene_annot = genes,
#'   exon_annot = exonREs,
#'   intron_annot = intronREs,
#'   intergenic_annot = intergenicREs,
#'   mc.cores = 12
#' )
#' }
process_quant_file <- function(quant_file, gene_annot, exon_annot, intron_annot, intergenic_annot) {
  DT <- data.table::fread(quant_file, sep = "\t")
  DT <- DT[NumReads > 0]

  # Inner join the gene annotation to the quants
  transcript_counts <-
    data.table::merge.data.table(
      x = DT,
      y = gene_annot,
      by.x = "Name",
      by.y = "md5",
      all.x = FALSE,
      all.y = FALSE
    )

  # summarize to the gene name level and round to integers
  gene_counts <- transcript_counts[, .(count = round(sum(NumReads, na.rm = TRUE), digits = 0)), by = .(symbol)]

  # Inner join counts with exonic RE annotation ---
  sample_exon_dt <-
    data.table::merge.data.table(
      x = DT,
      y = exon_annot,
      by.x = "Name",
      by.y = "md5",
      all.x = FALSE,
      all.y = FALSE
    )

  exon_counts <-
    sample_exon_dt[
      ,
      repElem := paste(repClass, repFamily, repName, sep = ".")
    ][,
      .(count = sum(NumReads, na.rm = TRUE)),
      by = "repElem"
    ][
      ,
      count := round(count, digits = 0)
    ]

  # Inner join counts with Intronic annotation ---
  sample_intron_dt <-
    data.table::merge.data.table(
      x = DT,
      y = intron_annot,
      by.x = "Name",
      by.y = "md5",
      all.x = FALSE,
      all.y = FALSE
    )

  intron_counts <- sample_intron_dt[
    ,
    repElem := paste(repClass, repFamily, repName, sep = ".")
  ][,
    .(count = sum(NumReads, na.rm = TRUE)),
    by = "repElem"
  ][
    ,
    count := round(count, digits = 0)
  ]

  # Inner join counts with Intronic annotation ---
  sample_intergenic_dt <-
    data.table::merge.data.table(
      x = DT,
      y = intergenic_annot,
      by.x = "Name",
      by.y = "md5",
      all.x = FALSE,
      all.y = FALSE
    )

  intergenic_counts <- sample_intergenic_dt[
    ,
    repElem := paste(repClass, repFamily, repName, sep = ".")
  ][,
    .(count = sum(NumReads, na.rm = TRUE)),
    by = "repElem"
  ][
    ,
    count := round(count, digits = 0)
  ]

  # combine the intronic and intergenic counts into a single data.table ---
  re_counts <- data.table::rbindlist(
    list(intron_counts, intergenic_counts)
  )[,
    .(count = sum(count)),
    by = "repElem"
  ]

  data_list <- list(
    "gene" = gene_counts,
    "exon" = exon_counts,
    "intron" = intron_counts,
    "intergenic" = intergenic_counts,
    "re" = re_counts
  )

  return(data_list)
}
