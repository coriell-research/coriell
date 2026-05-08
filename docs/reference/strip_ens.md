# Strip version IDs from ENSEMBL identifiers

This function is simply a wrapper around gsub for removing trailing IDs
from ENSEMBL identifiers. All entries must start with 'ENS'. NAs are
tolerated.

## Usage

``` r
strip_ens(x)
```

## Arguments

- x:

  character vector containing ENSEMBL IDs

## Value

character vector with trailing IDs removed

## Examples

``` r
ids <- c("ENSG0000001.12", "ENSMUSG00021.3", "ENST00000556.2")
strip_ens(ids)
#> [1] "ENSG0000001"  "ENSMUSG00021" "ENST00000556"
```
