# Perform and exact correlation test on every row of a matrix

Perform an exact correlation permutation test. Computes correlation
value for every row of matrix X against all permutations of vector y.

## Usage

``` r
exact_cor_test(X, y, ...)
```

## Arguments

- X:

  numeric matrix or data.frame that can be converted to a numeric matrix

- y:

  numeric vector of data to correlate

- ...:

  arguments passed to the \`cor\` function

## Value

data.frame containing the original values in X along with columns
containing the correlation of Xi and Y, the empirical p-value of the
permutation test, and the FDR corrected empirical p-value
