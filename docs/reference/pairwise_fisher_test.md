# Perform pairwise fisher.tests on an input matrix

This function will perform a
[`fisher.test()`](https://rdrr.io/r/stats/fisher.test.html) for all
columns relative to the given reference column of the input 2 x m matrix
and return a data.frame of the results.

## Usage

``` r
pairwise_fisher_test(x, ref)
```

## Arguments

- x:

  numeric 2 x m data.frame, matrix, or table

- ref:

  column name of the reference column of the input matrix

## Value

data.frame of results for each test.

## Details

The function performs a fisher test for each test column relative to the
reference column of the input matrix. The input matrix must be a 2 x m
matrix. It is the user's responsibility to ensure the levels of the rows
align with the desired odds ratio calculation. e.g. R will calculate the
odds ratio as
`(x[1, "test_col"] / x[2, "test_col"]) / (x[1, "ref_col"] / x[2, "ref_col"])`.
It is the user's responsibility to make sure `x[1, ]` and `x[2, ]` are
in the desired order.

## Examples

``` r
# From chisq.test docs:
M <- as.table(rbind(c(762, 327, 468), c(484, 239, 477)))
dimnames(M) <- list(
  gender = c("F", "M"),
  party = c("Democrat", "Independent", "Republican")
)
pairwise_fisher_test(M, ref = "Independent")
#>                  comparison odds_ratio    ci_low   ci_high     p_value
#> 1   Democrat vs Independent  1.1505995 0.9349959 1.4151128 0.178627092
#> 2 Republican vs Independent  0.7172647 0.5779167 0.8894628 0.002018961
#>   p_value_adj
#> 1 0.357254185
#> 2 0.004037923
```
