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
#>             sample.1  sample.2   sample.3  sample.4  sample.5   sample.6
#> feature.1 0.88776958 0.5287509 0.52534493 0.4714122 0.5600067 0.05079063
#> feature.2 0.60973271 0.7499077 0.13934146 0.7577911 0.9422215 0.15216337
#> feature.3 0.72935218 0.6258708 0.09341653 0.6242353 0.9244045 0.21640077
#> feature.4 0.30077489 0.7821802 0.72793197 0.2810700 0.5259212 0.38674096
#> feature.5 0.34060712 0.2926655 0.32752380 0.3647919 0.9685592 0.52973834
#> feature.6 0.02557859 0.1189060 0.62369842 0.7502761 0.5291743 0.09080587
#>            sample.7  sample.8  sample.9  sample.10         cor empirical.p
#> feature.1 0.3677616 0.1206362 0.8202867 0.61420133 -0.12727273       0.357
#> feature.2 0.7759316 0.8384177 0.2240106 0.15639422  0.01818182       0.470
#> feature.3 0.9740257 0.4975912 0.9488369 0.02835764 -0.07878788       0.394
#> feature.4 0.3350775 0.1998726 0.1103619 0.07449522 -0.70909091       0.014
#> feature.5 0.8874792 0.9321687 0.8972810 0.45942028  0.64848485       0.022
#> feature.6 0.4261961 0.4126758 0.5614466 0.86681060  0.44242424       0.084
#>                 FDR
#> feature.1 0.4862385
#> feature.2 0.4876289
#> feature.3 0.4862385
#> feature.4 0.4102564
#> feature.5 0.4735294
#> feature.6 0.4735294
```
