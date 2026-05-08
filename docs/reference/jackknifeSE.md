# Perform jackknife resampling on all columns of a SummarizeExperiment object

This function provides a simple wrapper to perform jackknife resampling
on all columns of a SummarizedExperiment object and returns the results
of each resample in a list. This function is designed to be used to
assess the robustness of p-value combination techniques of the included
[`meta_de()`](https://coriell-research.github.io/coriell/reference/meta_de.md)
function but in theory any arbitrary function which operates on the
columns of a SummarizedExperiment object could be used.

## Usage

``` r
jackknifeSE(x, FUN, ...)
```

## Arguments

- x:

  SummarizedExperiment object to perform jackknife resampling of columns
  on

- FUN:

  Function to perform on each resample.

- ...:

  Additional arguments passed to FUN

## Value

List of jackknife resampled results

## Examples

``` r
# Define three differential expression dataset data.frames
exp1 <- data.frame(
  feature_id = c("geneA", "geneB", "geneC"),
  PValue = c(0.01, 0.5, 0.05),
  FDR = c(0.02, 0.5, 0.07),
  logFC = c(1.2, -2.5, 3.7),
  logCPM = c(12, 9, 0)
)

exp2 <- data.frame(
  feature_id = c("geneA", "geneB", "geneD"),
  PValue = c(0.07, 0.3, 0.8),
  FDR = c(0.08, 0.4, 1.0),
  logFC = c(1.5, -2.0, 3.0),
  logCPM = c(14, 10, 2)
)

exp3 <- data.frame(
  feature_id = c("geneA", "geneB", "geneC", "geneD"),
  PValue = c(0.03, 0.3, 0.01, 0.8),
  FDR = c(0.08, 0.4, 0.04, 0.9),
  logFC = c(1.5, -2.0, 3.0, 4.1),
  logCPM = c(14, 10, 1, 2.1)
)

# Combine into a single list
l <- list(experiment1 = exp1, experiment2 = exp2, experiment3 = exp3)

# Convert the data to a SummarizedExperiment
se <- dfs2se(l)

# Perform the jackknife using meta_de on each subset of the data
metafun <- function(x) { meta_de(x, metapod::parallelFisher) }
result <- jackknifeSE(se, FUN = metafun)

# Combine the results from calling meta_de on each resample and show
result <- data.table::rbindlist(result, idcol = "Jackknife")
head(result[order(Feature)])
#>    Jackknife Feature Combined.Pval Direction Rep.logFC Rep.Pval Median.logFC
#>        <int>  <char>         <num>    <char>     <num>    <num>        <num>
#> 1:         1   geneA   0.015048218        up       1.5     0.03         1.50
#> 2:         2   geneA   0.002733518        up       1.2     0.01         1.35
#> 3:         3   geneA   0.005785101        up       1.2     0.01         1.35
#> 4:         1   geneB   0.306715105      down      -2.0     0.30        -2.00
#> 5:         2   geneB   0.434567998      down      -2.0     0.30        -2.25
#> 6:         3   geneB   0.434567998      down      -2.0     0.30        -2.25
#>    Mean.logFC Min.logFC Max.logFC Meta.logFC Meta.Pval Meta.z
#>         <num>     <num>     <num>     <lgcl>    <lgcl> <lgcl>
#> 1:       1.50       1.5       1.5         NA        NA     NA
#> 2:       1.35       1.2       1.5         NA        NA     NA
#> 3:       1.35       1.2       1.5         NA        NA     NA
#> 4:      -2.00      -2.0      -2.0         NA        NA     NA
#> 5:      -2.25      -2.5      -2.0         NA        NA     NA
#> 6:      -2.25      -2.5      -2.0         NA        NA     NA
```
