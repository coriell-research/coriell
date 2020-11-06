#' Perform GO Analysis with PANTHER
#'
#' Sends a request to PANTHER GO db to perform over representation analysis. This function 
#' excludes the option to import a reference list and reference organism. By default,
#' in this case, PANTHER will use all of the genes of the given organism as the 
#' reference list.
#' @param gene_list character vector. Maximum of 100,000 Identifiers. Can be any 
#' of the following: Ensemble gene identifier, Ensemble protein identifier, 
#' Ensemble transcript identifier, Entrez gene id, gene symbol, NCBI GI, HGNC Id, 
#' International protein index id, NCBI UniGene id, UniProt accession and 
#' UniProt id 
#'@param organism character string. Taxon ID (e.g. "9606" for HUMAN, "10090" for 
#'MOUSE, "10116" for RAT). To get list of available taxon IDs see:
#'`curl -X GET "http://pantherdb.org/services/oai/pantherdb/supportedgenomes" -H  "accept: application/json"` 
#'@param annot_dataset character string. One of c("biological_process", 
#'"molecular_function", "cellular_component"). 
#'@param enrichment_test_type character string. One of c("fisher", "binomial").
#'Default "fisher"
#'@param correction character string. One of c("fdr", "bonferroni", "none").
#'Default "fdr"
#'@return list with results of GO analysis and input information. See PANTHER user manual 
#'\url{http://www.pantherdb.org/help/PANTHER_user_manual.pdf} for details
#'@export
#'@examples
#' library(coriell)
#'
#'
#' genes <- c("CTNNB1", "ADAM17", "AXIN1", "AXIN2", "CCND2", "CSNK1E", "CTNNB1", 
#'            "CUL1", "DKK1", "DKK4", "DLL1", "DVL2", "FRAT1", "FZD1", "FZD8", 
#'            "GNAI1", "HDAC11", "HDAC2", "HDAC5", "HEY1", "HEY2", "JAG1", 
#'            "JAG2", "KAT2A", "LEF1", "MAML1", "MYC", "NCOR2", "NCSTN", 
#'            "NKD1", "NOTCH1", "NOTCH4", "NUMB", "PPARD", "PSEN2", "PTCH1", 
#'            "RBPJ", "SKP2", "TCF7", "TP53", "WNT1", "WNT5B", "WNT6")
#'
#' go_results <- panther_go(genes, organism = "9606", annot_dataset = "biological_process")
#' 
#' # to view results access the 'results' item
#' head(go_results$results)
panther_go <- function(
                       gene_list,
                       organism,
                       annot_dataset = "biological_process",
                       enrichment_test_type = "fisher",
                       correction = "fdr") {
  datasets <- c("biological_process" = "GO%3A0008150", 
                "molecular_function" = "GO%3A0003674",
                "cellular_component" = "GO%3A0005575")
  tests <- c("fisher" = "FISHER", "binomial" = "BINOMIAL")
  corrections <- c("fdr" = "FDR", "bonferroni" = "BONFERRONI", "none" = "NONE")
  
  stopifnot("Empty gene_list" = length(gene_list) > 0)
  stopifnot("Too many genes in gene list" = length(gene_list) <= 100000)
  stopifnot("Only one organism identifier should be provided" = length(organism) == 1)
  stopifnot("annot_dataset must be one of c('biological_process', 'molecular_function', 'cellular_component')" = annot_dataset %in% names(datasets))
  stopifnot("enrichment_test_type must be one of c('fisher', 'binomial')" = enrichment_test_type %in% names(tests))
  stopifnot("correction must be one of c('fdr', 'bonferroni', 'none')" = correction %in% names(corrections))
  
  base_url <- "http://pantherdb.org/services/oai/pantherdb/enrich/overrep?"
  gene_input <- paste0("geneInputList=", paste(gene_list, collapse = ","))
  organism_input <- paste0("&organism=", organism)
  annot_input <- paste0("&annotDataSet=", datasets[annot_dataset])
  test_input <- paste0("&enrichmentTestType=", tests[enrichment_test_type])
  correction_input <- paste0("&correction=", corrections[correction])
  req_url <- paste0(base_url, gene_input, organism_input, annot_input, test_input, correction_input)
  
  r <- httr::GET(req_url)
  httr::stop_for_status(r, "ERROR: failed to execute request")
  httr::warn_for_status(r, "WARNING: request produced a warning response")
  
  parsed <- jsonlite::fromJSON(httr::content(r, "text"), simplifyVector = FALSE)
  
  results <- purrr::pluck(parsed, "results", "result") %>% 
    purrr::map(dplyr::as_tibble) %>% 
    dplyr::bind_rows(.id = "result_number") %>% 
    tidyr::unnest(term) %>%
    dplyr::group_by(result_number) %>% 
    dplyr::mutate(ind = rep(c("GO_term", "description"), length.out = dplyr::n())) %>% 
    tidyr::pivot_wider(names_from = ind, values_from = term)
  
  list("results" = results,
       "gene_query" = gene_list,
       "organism_query" = organism,
       "test_type" = tests[enrichment_test_type],
       "correction" = corrections[correction],
       "request_url" = req_url)
}
