# Perform a permutation correlation test on every row of a matrix

Compute a correlation value for every row of X against the vector y and
n_perm random permutations of y. If the number of distinct permutations
of y is less than or equal to n_perm then an exact test using all
permutations is performed instead. In both cases the function returns a
data.frame of the original data with additional columns for the test
statistic, permutation p-value, and FDR corrected p-value.

## Usage

``` r
permutation_correlation_test(
  X,
  y,
  n_perm = 10000,
  method = c("pearson", "spearman"),
  alternative = c("two.sided", "greater", "less"),
  chunk_size = 1000,
  n_core = 1
)
```

## Arguments

- X:

  numeric matrix or data.frame that can be converted to a numeric matrix

- y:

  numeric vector of values to correlate with rows of X

- n_perm:

  integer. The desired number of permutations to sample from. Default
  (10,000)

- method:

  character. One of "pearson" (default) or "spearman"

- alternative:

  character. One of "two.sided" (default), "greater" or "less"

- chunk_size:

  integer. Number of permutations to process at once. Controls memory
  usage, nrow(X) x chunk_size correlations are held in memory at a time.
  Default (1000)

- n_core:

  integer. The number of cores to use for processing chunks of
  permutations. Default (1)

## Value

data.frame containing the original values in X along with columns 'cor',
'empirical.p', and 'FDR'

## Details

Exact p-values are computed as the proportion of all distinct
permutations (including the observed ordering) with a correlation at
least as extreme as the observed value. When permutations are randomly
sampled, p-values are computed with
[`statmod::permp()`](https://rdrr.io/pkg/statmod/man/permp.html) (if
installed), which gives exact p-values for permutations sampled with
replacement (Phipson & Smyth 2010). Otherwise the conservative estimate
(b + 1) / (n_perm + 1) is used.

Correlations are computed by standardizing the rows of X once and
performing matrix multiplication against chunks of permuted y. Rows of X
containing missing values or with zero variance return NA.

## References

Phipson B, Smyth GK (2010). Permutation P-values should never be zero:
calculating exact P-values when permutations are randomly drawn.
Statistical Applications in Genetics and Molecular Biology 9, Article
39.

## Examples

``` r
# generate example data
X <- matrix(runif(1e3 * 10), nrow = 1e3, ncol = 10)
y <- 1:10
dimnames(X) <- list(paste("feature", 1:1e3, sep = "."), paste("sample", 1:10, sep = "."))

# correlate each row of X with 1,000 random permutations of vector y
res <- permutation_correlation_test(X, y, n_perm = 1e3, method = "spearman")

head(res)
#>              sample.1   sample.2  sample.3  sample.4  sample.5  sample.6
#> feature.1 0.141592341 0.42516774 0.4146605 0.2634410 0.3691033 0.9263442
#> feature.2 0.635351351 0.70836595 0.2823288 0.6278566 0.3701241 0.9765515
#> feature.3 0.987373579 0.65636057 0.1578583 0.4449980 0.3436975 0.8037980
#> feature.4 0.459281017 0.68779913 0.3004809 0.2736790 0.3555447 0.8287041
#> feature.5 0.005317752 0.13298787 0.6195136 0.5304429 0.5265921 0.4056231
#> feature.6 0.202410232 0.02137243 0.3041225 0.2431722 0.2127211 0.1878105
#>             sample.7   sample.8  sample.9 sample.10         cor empirical.p
#> feature.1 0.68855969 0.28051254 0.7078689 0.6548090  0.53939394   0.1298700
#> feature.2 0.12599496 0.95200286 0.3859197 0.5272708 -0.11515152   0.7432566
#> feature.3 0.01177438 0.56333502 0.9124741 0.8295821  0.07878788   0.8311687
#> feature.4 0.26823535 0.52765485 0.9441315 0.5284524  0.29696970   0.3886113
#> feature.5 0.31443307 0.08967741 0.7765870 0.6204135  0.49090909   0.1728270
#> feature.6 0.53255829 0.12447965 0.3868276 0.3544599  0.46666667   0.1698300
#>                 FDR
#> feature.1 0.9204271
#> feature.2 0.9999999
#> feature.3 0.9999999
#> feature.4 0.9911808
#> feature.5 0.9257751
#> feature.6 0.9257751
```
