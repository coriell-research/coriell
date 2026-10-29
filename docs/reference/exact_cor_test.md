# Perform an exact correlation test on every row of a matrix

Perform an exact correlation permutation test. Computes correlation
value for every row of matrix X against all distinct permutations of
vector y.

## Usage

``` r
exact_cor_test(
  X,
  y,
  method = c("pearson", "spearman"),
  alternative = c("two.sided", "greater", "less")
)
```

## Arguments

- X:

  numeric matrix or data.frame that can be converted to a numeric matrix

- y:

  numeric vector of data to correlate

- method:

  character. One of "pearson" (default) or "spearman"

- alternative:

  character. One of "two.sided" (default), "greater" or "less"

## Value

data.frame containing the original values in X along with columns
containing the correlation of Xi and Y, the exact permutation p-value,
and the FDR corrected p-value
