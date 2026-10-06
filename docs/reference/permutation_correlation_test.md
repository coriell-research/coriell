# Perform a permutation correlation test on every row of a matrix

Compute a correlation value for every row of X against the vector y and
n random permutations of y. If the number of possible permutations is
less than the the argument n_perm then an exact test is performed
instead. In both cases the function returns a data.frame of the original
data with additional columns for the test statistic, empirical p-value,
and FDR corrected empirical p-value.

## Usage

``` r
permutation_correlation_test(X, y, n_perm = 10000, n_core = 1, ...)
```

## Arguments

- X:

  numeric matrix or data.frame that can be converted to a numeric matrix

- y:

  numeric vector of values to correlate with rows of X

- n_perm:

  integer. The desired number of permutations to sample from. Default
  (10,000)

- n_core:

  integer. The number of cores to use for processing. Default (1)

- ...:

  Additional arguments to pass to \`cor\` function

## Examples

``` r
# generate example data
X <- matrix(runif(1e3 * 10), nrow = 1e3, ncol = 10)
y <- 1:10
dimnames(X) <- list(paste("feature", 1:1e3, sep = "."), paste("sample", 1:10, sep = "."))

# correlate each row of X with 1,000 random permutations of vector y
res <- permutation_correlation_test(X, y, n_perm = 1e3, n_core = 8, method = "spearman")

head(res)
#>            sample.1  sample.2  sample.3   sample.4  sample.5  sample.6
#> feature.1 0.8083622 0.8383412 0.1309728 0.08066313 0.4570966 0.8303105
#> feature.2 0.1862091 0.6485302 0.5776029 0.73453357 0.6566667 0.3617257
#> feature.3 0.6550341 0.7636185 0.6608876 0.49824185 0.3085832 0.5314102
#> feature.4 0.6083527 0.1596689 0.1662174 0.57505885 0.9627974 0.5769995
#> feature.5 0.8929753 0.7293678 0.4999551 0.20110452 0.7460042 0.0604303
#> feature.6 0.1571574 0.1973273 0.2653899 0.12052447 0.6094003 0.7215492
#>            sample.7   sample.8   sample.9 sample.10         cor empirical.p
#> feature.1 0.5623403 0.39973106 0.76813744 0.8876729  0.16363636       0.295
#> feature.2 0.7827487 0.79909107 0.05683166 0.7807953  0.35757576       0.147
#> feature.3 0.5906869 0.47597199 0.18499534 0.5399923 -0.60000000       0.037
#> feature.4 0.1440349 0.56141709 0.13912811 0.3610715 -0.34545455       0.160
#> feature.5 0.8143201 0.19901102 0.91171531 0.5939837 -0.05454545       0.433
#> feature.6 0.6647761 0.07689694 0.27017436 0.1215346  0.01818182       0.440
#>                 FDR
#> feature.1 0.4851974
#> feature.2 0.4608434
#> feature.3 0.4062500
#> feature.4 0.4655172
#> feature.5 0.4963274
#> feature.6 0.4963274
```
