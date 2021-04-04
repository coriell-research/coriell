#' Process a quant file from REdiscoverTE pipeline
#'
#' Specialty function for processing quants files produced by the \href{https://www.nature.com/articles/s41467-019-13035-2}{REdiscoverTE}
#' pipeline. The function takes a quant file as input and calculates the counts
#' per gene and per repetitive element subfamily and returns both count data.frames
#' in a list. The user must provide the pre-processed gene and RepeatMasker
#' annotations.
#'
#' @param quant_file path to the quant file to be processed
#' @param gene_annot REdiscoverTE gene annotation data.frame
#' @param rmsk_annot REdiscoverTE RepeatMasker annotation data.frame
#' @return list containing two data.tables. One for the gene counts ("gene") and the other for repeat counts ("re")
#' @export
#' @examples
#' \dontrun{
#' # read in annotation files
#' gencode <- as.data.table(readRDS("/mnt/data/gdata/human/REdiscoverTE_hg38/REdiscoverTE_hg38_GFP/GENCODE.V26.Basic_Gene_Annotation_md5_GFP.rds"))
#' rmsk <- as.data.table(readRDS("/mnt/data/gdata/human/REdiscoverTE_hg38/rmsk_annotation.RDS"))
#'
#' # specify path to quant file
#' my_quant_file <- "/path/to/quant.sf"
#'
#' result <- process_quant_file(my_quant_file, gencode, rmsk)
#'
#' # If you have a list of quant files, you can apply over all with multiple cores
#' my_quants <- list.files("path/to/quants_dir", pattern = "quant.sf", full.names = TRUE)
#'
#' processed <- parallel::mclapply(
#'   X = my_quants,
#'   FUN = process_quant_file,
#'   gene_annot = gencode,
#'   rmsk_annot = rmsk,
#'   mc.cores = 12
#' )
#' }
process_quant_file <- function(quant_file, gene_annot, rmsk_annot) {
  DT <- data.table::fread(quant_file, sep = "\t")
  DT <- DT[NumReads > 0]

  transcript_counts <- data.table::merge.data.table(
    x = DT,
    y = gene_annot,
    by.x = "Name",
    by.y = "md5"
  )

  repElem_counts <- data.table::merge.data.table(
    x = DT,
    y = rmsk_annot,
    by.x = "Name",
    by.y = "md5"
  )

  gene_counts <- transcript_counts[, .(count = sum(NumReads, na.rm = TRUE)), by = .(symbol)]
  repElem_counts <- repElem_counts[, .(repElem = paste(repClass, repFamily, repName, sep = "."), NumReads)][, .(count = sum(NumReads, na.rm = TRUE)), by = .(repElem)]

  list("gene" = gene_counts, "re" = repElem_counts)
}
