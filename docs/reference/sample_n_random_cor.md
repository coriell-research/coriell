# Generate a null distribution of correlation values

Selects n random rows from a numeric matrix with replacement. For each
random row, permute vector y and perform correlation.

## Usage

``` r
sample_n_random_cor(X, y, n = 10000, ...)
```

## Arguments

- X:

  numeric matrix or data.frame that can be converted to numeric matrix.

- y:

  numeric vector. Numeric vector of values used to correlate with each
  row of df

- n:

  integer. Number of random correlations to return. Default (1e4)

- ...:

  Additional arguments passed to \`cor\` function

## Value

numeric vector of correlation values. Vector can contain NAs if row
variance == 0.

## Examples

``` r
# generate example data
X <- matrix(runif(1e3 * 6), nrow = 1e3, ncol = 6)
y <- 1:6
dimnames(X) <- list(paste("feature", 1:1e3, sep = "."), paste("sample", 1:6, sep = "."))

# sample random correlations from permuted data
null_dist <- sample_n_random_cor(X, y, n = 1e3, method = "spearman")

head(null_dist)
#> feature.333 feature.718 feature.797 feature.335 feature.122  feature.71 
#> -0.42857143 -0.20000000 -0.42857143 -0.08571429 -0.48571429  0.77142857 
```
