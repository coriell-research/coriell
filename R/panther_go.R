#' Perform GO analysis with PANTHER
#'
#' This statistical test tool, compares a test gene list to a reference gene list,
#' and determines whether a particular class (e.g. molecular function, biological 
#' process, cellular component, PANTHER protein class, the PANTHER pathway or 
#' Reactome pathway) of genes is overrepresented or underrepresented.
#' 
#' @details
#' Sends a request to \href{http://pantherdb.org/}{PANTHER db} to perform over
#' representation analysis. This function excludes the option to import a reference
#' list and reference organism. By default, in this case, PANTHER will use all
#' of the genes of the given organism as the reference list.
#' 
#' @param gene_list character vector. Maximum of 100,000 identifiers. Can be any
#' of the following: Ensemble gene identifier, Ensemble protein identifier,
#' Ensemble transcript identifier, Entrez gene id, gene symbol, NCBI GI, HGNC Id,
#' International protein index id, NCBI UniGene id, UniProt accession and
#' UniProt id
#' @param organism character string. Taxon ID (e.g. "9606" for HUMAN, "10090" for
#' MOUSE, "10116" for RAT). To get list of available taxon IDs see:
#' \preformatted{curl -X GET "https://pantherdb.org/services/oai/pantherdb/supportedgenomes" -H  "accept: application/json"}
#' @param annot_dataset character string. One of c("biological_process",
#' "molecular_function", "cellular_component", "panther_go_slim_mf", "panther_go_slim_bp",
#' "panther_go_slim_cc", "panther_pc", "panther_pathway", "panther_reactome_pathway"). see:
#' \preformatted{curl -X POST "https://pantherdb.org/services/oai/pantherdb/supportedannotdatasets" -H "accept: application/json"}
#' for full descriptions.
#' @param enrichment_test_type character string. One of c("fisher", "binomial").
#' Default "fisher"
#' @param correction character string. One of c("fdr", "bonferroni", "none").
#' Default "fdr"
#' @return data.table of results from over representation analysis.
#' See \href{http://www.pantherdb.org/help/PANTHER_user_manual.pdf}{PANTHER user manual}
#' for column descriptions in "table".
#' @import data.table
#' @export
#' @examples
#' genes <- c(
#'   "CTNNB1", "ADAM17", "AXIN1", "AXIN2", "CCND2", "CSNK1E", "CTNNB1",
#'   "CUL1", "DKK1", "DKK4", "DLL1", "DVL2", "FRAT1", "FZD1", "FZD8",
#'   "GNAI1", "HDAC11", "HDAC2", "HDAC5", "HEY1", "HEY2", "JAG1",
#'   "JAG2", "KAT2A", "LEF1", "MAML1", "MYC", "NCOR2", "NCSTN",
#'   "NKD1", "NOTCH1", "NOTCH4", "NUMB", "PPARD", "PSEN2", "PTCH1",
#'   "RBPJ", "SKP2", "TCF7", "TP53", "WNT1", "WNT5B", "WNT6"
#' )
#'
#' go_results <- panther_go(genes, organism = "9606", annot_dataset = "biological_process")
#'
#' # to view results access the 'table' item
#' head(go_results$table)
#'
#' # to view the actual request (which allows you to parse results yourself) use
#' req <- go_results$request
panther_go <- function(gene_list,
                       organism,
                       annot_dataset,
                       enrichment_test_type = "fisher",
                       correction = "fdr") {
  if (!requireNamespace("httr2", quietly = TRUE)) {
    stop("httr2 package is required.")
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("jsonlite package is required.")
  }
  
  datasets <- c(
    "biological_process" = "GO:0008150",
    "molecular_function" = "GO:0003674",
    "cellular_component" = "GO:0005575",
    "panther_go_slim_mf" = "ANNOT_TYPE_ID_PANTHER_GO_SLIM_MF",
    "panther_go_slim_bp" = "ANNOT_TYPE_ID_PANTHER_GO_SLIM_BP",
    "panther_go_slim_cc" = "ANNOT_TYPE_ID_PANTHER_GO_SLIM_CC",
    "panther_pc" = "ANNOT_TYPE_ID_PANTHER_PC",
    "panther_pathway" = "ANNOT_TYPE_ID_PANTHER_PATHWAY",
    "panther_reactome_pathway" = "ANNOT_TYPE_ID_REACTOME_PATHWAY"
  )
  tests <- c("fisher" = "FISHER", "binomial" = "BINOMIAL")
  corrections <- c("fdr" = "FDR", "bonferroni" = "BONFERRONI", "none" = "NONE")
  
  stopifnot("Empty gene_list" = length(gene_list) > 0)
  stopifnot("Too many genes in gene list" = length(gene_list) <= 100000)
  stopifnot("Only one organism identifier should be provided" = length(organism) == 1)
  stopifnot("Incorrectly specified dataset" = annot_dataset %in% names(datasets))
  stopifnot("enrichment_test_type must be one of c('fisher', 'binomial')" = enrichment_test_type %in% names(tests))
  stopifnot("correction must be one of c('fdr', 'bonferroni', 'none')" = correction %in% names(corrections))
  
  base_url <- "https://pantherdb.org/services/oai/pantherdb/enrich/overrep?"
  gene_input <- paste(gene_list, collapse = ",")
  organism_input <- organism
  annot_input <- unname(datasets[annot_dataset])
  test_input <- unname(tests[enrichment_test_type])
  correction_input <- unname(corrections[correction])
  data <- list(
    geneInputList = gene_input,
    organism = organism_input,
    annotDataSet = annot_input,
    enrichmentTestType = test_input,
    correction = correction_input
  )
  
  req <- httr2::request(base_url)
  req <- httr2::req_body_json(req, data)
  return(req)
  resp <- httr2::req_perform(req)
  
  
  if (isTRUE(httr2::resp_is_error(resp)))
    stop("An error occured in the request")
  
  json <- httr2::resp_body_json(resp)
  return(json)
  parsed <- jsonlite::fromJSON(json, simplifyVector = FALSE)
  dt <- data.table::rbindlist(parsed$results$result)
  
  return(dt)
}
