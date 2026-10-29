# Perform GO analysis with PANTHER

This statistical test tool, compares a test gene list to a reference
gene list, and determines whether a particular class (e.g. molecular
function, biological process, cellular component, PANTHER protein class,
the PANTHER pathway or Reactome pathway) of genes is overrepresented or
underrepresented.

## Usage

``` r
panther_go(
  gene_list,
  organism,
  annot_dataset,
  ref_input_list = NULL,
  enrichment_test_type = "fisher",
  correction = "fdr",
  verbose = 0
)
```

## Arguments

- gene_list:

  character vector. Maximum of 100,000 identifiers. Can be any of the
  following: Ensemble gene identifier, Ensemble protein identifier,
  Ensemble transcript identifier, Entrez gene id, gene symbol, NCBI GI,
  HGNC Id, International protein index id, NCBI UniGene id, UniProt
  accession and UniProt id

- organism:

  character string. Taxon ID (e.g. "9606" for HUMAN, "10090" for MOUSE,
  "10116" for RAT). To get list of available taxon IDs see:

      curl -X GET "https://pantherdb.org/services/oai/pantherdb/supportedgenomes" -H  "accept: application/json"

- annot_dataset:

  character string. One of c("biological_process", "molecular_function",
  "cellular_component", "panther_go_slim_mf", "panther_go_slim_bp",
  "panther_go_slim_cc", "panther_pc", "panther_pathway",
  "panther_reactome_pathway"). see:

      curl -X POST "https://pantherdb.org/services/oai/pantherdb/supportedannotdatasets" -H "accept: application/json"

  for full descriptions.

- ref_input_list:

  character vector. Reference set of genes for the specified organism
  (e.g. all genes expressed in the experiment). Must contain all of the
  genes in `gene_list`. Maximum of 100,000 identifiers. If NULL
  (default) then PANTHER will use all genes for the specified organism.

- enrichment_test_type:

  character string. One of c("fisher", "binomial"). Default "fisher"

- correction:

  character string. One of c("fdr", "bonferroni", "none"). Default "fdr"

- verbose:

  integer. Verbosity level passed to
  [`httr2::req_perform()`](https://httr2.r-lib.org/reference/req_perform.html).
  Default 0

## Value

data.table of results from over representation analysis with one row per
term. The term returned by PANTHER is split into "term_id" and
"term_label" columns ("term_id" is NA for UNCLASSIFIED). See [PANTHER
user manual](http://www.pantherdb.org/help/PANTHER_user_manual.pdf) for
column descriptions.

## Details

Sends a request to [PANTHER db](http://pantherdb.org/) to perform over
representation analysis. By default (`ref_input_list = NULL`) PANTHER
will use all of the genes of the given organism as the reference list.
If a custom reference list is supplied then every identifier in
`gene_list` must also be present in `ref_input_list`, otherwise PANTHER
will refuse to process the request. Identifiers that PANTHER cannot map
are ignored.

## Examples

``` r
genes <- c(
  "CTNNB1", "ADAM17", "AXIN1", "AXIN2", "CCND2", "CSNK1E", "CTNNB1",
  "CUL1", "DKK1", "DKK4", "DLL1", "DVL2", "FRAT1", "FZD1", "FZD8",
  "GNAI1", "HDAC11", "HDAC2", "HDAC5", "HEY1", "HEY2", "JAG1",
  "JAG2", "KAT2A", "LEF1", "MAML1", "MYC", "NCOR2", "NCSTN",
  "NKD1", "NOTCH1", "NOTCH4", "NUMB", "PPARD", "PSEN2", "PTCH1",
  "RBPJ", "SKP2", "TCF7", "TP53", "WNT1", "WNT5B", "WNT6"
)

result <- panther_go(genes, "9606", "biological_process")
head(result)
#>    number_in_list fold_enrichment          fdr  expected number_in_reference
#>             <int>           <num>        <num>     <num>               <int>
#> 1:             30        7.302534 1.988037e-17 4.1081633                2013
#> 2:             11      110.000000 1.161814e-16 0.1000000                  49
#> 3:             12       78.400000 1.013277e-16 0.1530612                  75
#> 4:             31        5.111036 2.212363e-14 6.0653061                2972
#> 5:             30        5.390539 2.414309e-14 5.5653061                2727
#> 6:             12       47.040000 3.193091e-14 0.2551020                 125
#>          pValue plus_minus    term_id                         term_label
#>           <num>     <char>     <char>                             <char>
#> 1: 1.444689e-21          + GO:0048513           animal organ development
#> 2: 1.688561e-20          + GO:0007219            Notch signaling pathway
#> 3: 2.209020e-20          + GO:0060070    canonical Wnt signaling pathway
#> 4: 6.430819e-18          + GO:0007275 multicellular organism development
#> 5: 8.772287e-18          + GO:0048731                 system development
#> 6: 1.392235e-17          + GO:0016055              Wnt signaling pathway
```
