# Perform p-value combination for sets of differential expression tests

This function performs p-value combination for all genes and estimates
summary statistics for average effect sizes for all experiments in the
input. The p-value combination methods are performed using the `metapod`
package. The direction of effect is summarized using the representative
values determined by the given p-value combination function. See the
`metapod` package for details. If standard error values are available
for the logFC then these values can be used to compute fixed effect
meta-analysis statistics.

## Usage

``` r
meta_de(
  x,
  FUN,
  pval = "PValue",
  lfc = "logFC",
  se = NULL,
  impute_missing = TRUE,
  ...
)
```

## Arguments

- x:

  SummarizedExperiment object containing combined differential
  expression results from different studies. The SE object must contain
  at least two assays, one for the P-values to combine and the other for
  the effect sizes to compute (e.g. logFC).

- FUN:

  One of the 'parallel' functions provided by `metapod`. One of
  "parallelBerger", "parallelFisher", "parallelHolmMin",
  "parallelPearson", "parallelSimes", "parallelStouffer", or
  "parallelWilkinson".

- pval:

  assay name in SE object containing the P-values to combine.

- lfc:

  assay name in the SE object containing the logFC values to combine.

- se:

  assay name in SE object containing standard error of logFC values.
  Default NULL.

- impute_missing:

  TRUE/FALSE should missing values in the logFC and P-Value assays be
  imputed prior p-value combination? Default TRUE, missing p-values are
  imputed with 1 and missing logFCs are imputed with 0.

- ...:

  Additional arguments passed to FUN. See the `metapod` package for
  details.

## Value

data.table with summary stats of the p-value combination of all
experiments. Please see the documentation in the `metapod` package for
more details. The returned columns, "Rep.LogFC" and "Rep.Pval" contain
the results of extracting the representative effect and p-value from all
influential tests. These are individual tests in the data that are
particularly important for calculating the combined effects. The
"meta-analysis" results returned when standard errors are provided
assume fixed effects. This may not be appropriate.

## Details

`SummarizedExperiment` object.

## Examples

``` r
exp1 <- data.frame(
  feature_id = c("geneA", "geneB", "geneC"),
  PValue = c(0.01, 0.04, 0.5),
  FDR = c(0.02, 0.06, 0.6),
  logFC = c(2.0, -1.5, 0.1),  # geneA UP, geneB DOWN
  SE = c(0.5, 0.4, 1.0),
  logCPM = c(12, 11, 5)
)

# Experiment 2 (Similar to Exp 1)
exp2 <- data.frame(
  feature_id = c("geneA", "geneB", "geneD"),
  PValue = c(0.02, 0.05, 0.8),
  FDR = c(0.03, 0.07, 0.9),
  logFC = c(2.2, -1.4, 0.2),  # Very close to Exp 1
  SE = c(0.6, 0.45, 1.2),
  logCPM = c(12.5, 10.5, 4)
)

# Experiment 3 (Similar to Exp 1 & 2)
exp3 <- data.frame(
  feature_id = c("geneA", "geneB", "geneC"),
  PValue = c(0.005, 0.01, 0.4),
  FDR = c(0.01, 0.02, 0.5),
  logFC = c(1.9, -1.6, -0.1), # Very close to Exp 1 & 2
  SE = c(0.4, 0.35, 0.9),
  logCPM = c(11.8, 11.2, 5.2)
)

# Combine into a single list
l <- list(experiment1 = exp1, experiment2 = exp2, experiment3 = exp3)

# Convert the data to a SummarizedExperiment
se <- dfs2se(l, import = c("PValue", "logFC", "SE"))

# Perform p-value combination across experiments for each gene
# Perform meta-analysis on SEs by passing assay name
result <- meta_de(se, metapod::parallelFisher, se = "SE")
head(result)
#>    Feature Combined.Pval Direction Rep.logFC Rep.Pval Median.logFC  Mean.logFC
#>     <char>         <num>    <char>     <num>    <num>        <num>       <num>
#> 1:   geneA  0.0001102497        up       1.9    0.005          2.0  2.03333333
#> 2:   geneB  0.0014070716      down      -1.6    0.010         -1.5 -1.50000000
#> 3:   geneC  0.7809166219     mixed      -0.1    0.400          0.0  0.00000000
#> 4:   geneD  0.9984320588        up       0.2    0.800          0.0  0.06666667
#>    Min.logFC Max.logFC  Meta.logFC    Meta.Pval      Meta.z
#>        <num>     <num>       <num>        <num>       <num>
#> 1:       1.9       2.2  1.99466951 6.040596e-13  7.19956273
#> 2:      -1.6      -1.4 -1.51666531 2.525667e-11 -6.67186861
#> 3:      -0.1       0.1 -0.01049724 9.874803e-01 -0.01569177
#> 4:       0.0       0.2  0.20000000 8.676323e-01  0.16666667
```
