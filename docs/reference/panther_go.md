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

  Reference set of genes for the specified organism. If NULL (default)
  then PANTHER will use all genes for the specified organism.

- enrichment_test_type:

  character string. One of c("fisher", "binomial"). Default "fisher"

- correction:

  character string. One of c("fdr", "bonferroni", "none"). Default "fdr"

## Value

data.table of results from over representation analysis. See [PANTHER
user manual](http://www.pantherdb.org/help/PANTHER_user_manual.pdf) for
column descriptions in "table".

## Details

Sends a request to [PANTHER db](http://pantherdb.org/) to perform over
representation analysis. This function excludes the option to import a
reference list and reference organism. By default, in this case, PANTHER
will use all of the genes of the given organism as the reference list.

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
#> 2:             30        7.302534 1.988037e-17 4.1081633                2013
#> 3:             11      110.000000 1.161814e-16 0.1000000                  49
#> 4:             11      110.000000 1.161814e-16 0.1000000                  49
#> 5:             12       78.400000 1.013277e-16 0.1530612                  75
#> 6:             12       78.400000 1.013277e-16 0.1530612                  75
#>          pValue                            term plus_minus
#>           <num>                          <list>     <char>
#> 1: 1.444689e-21                      GO:0048513          +
#> 2: 1.444689e-21        animal organ development          +
#> 3: 1.688561e-20                      GO:0007219          +
#> 4: 1.688561e-20         Notch signaling pathway          +
#> 5: 2.209020e-20                      GO:0060070          +
#> 6: 2.209020e-20 canonical Wnt signaling pathway          +
```
