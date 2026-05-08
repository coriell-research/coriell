# Summarize RNA-seq expression results

Summarize a results data.frame. Return data.frame of counts of
up/down/non-DE genes based on log-fold-change and significance values.

## Usage

``` r
summarize_dge(df, fdr_col = "FDR", lfc_col = "logFC", fdr = 0.05, lfc = 0)
```

## Arguments

- df:

  dataframe of results. Must have columns containing significance values
  and log-fold changes.

- fdr_col:

  Column name of the data.frame containing the significance level
  values.

- lfc_col:

  Column name of data.frame containing the log-fold change values.

- fdr:

  numeric. FDR or significance value below which genes are considered
  significant.

- lfc:

  numeric. abs(log-fold change) value above which genes are considered
  significant.

## Value

data.frame with columns describing the count and percentages of
up/down/unperturbed genes

## Examples

``` r
summarize_dge(GSE161650_de)
#>     Direction    N Percent
#> 1          Up 2376   21.60
#> 2        Down 3020   27.45
#> 3 Unperturbed 5606   50.95
```
