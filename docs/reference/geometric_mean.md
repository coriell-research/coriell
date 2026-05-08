# Geometric mean of a vector

Calculate the geometric mean of a vector. This function is valid for
non-negative, non-NA containing vectors.

## Usage

``` r
geometric_mean(x, zero_propagate = FALSE, ignore_zero = TRUE)
```

## Arguments

- x:

  numeric vector of non-negative values

- zero_propagate:

  logical. Should zeros be included in the calculation. Default FALSE.

- ignore_zero:

  logical. Should zero values be ignored in the calculation of the mean?
  Default TRUE.

## Value

geometric mean of the vector

## Examples

``` r
# Normal case
geometric_mean(c(2, 5, 95, 5))
#> [1] 8.301822

# Default ignores 0s entirely
geometric_mean(c(2, 5, 95, 5, 0, 0, 0, 0))
#> [1] 8.301822

# Ignore zero = FALSE -- zero is used in mean calculation
geometric_mean(c(2, 5, 95, 5, 0, 0, 0, 0), ignore_zero = FALSE)
#> [1] 2.881288

# Case with NA -- Returns NA
geometric_mean(c(NA, 1, 2, 3))
#> [1] NA

# Case with 0 propagation -- Returns 0
geometric_mean(c(0, 1, 2, 3), zero_propagate = TRUE)
#> [1] 0

# Case with negative -- Returns NaN
geometric_mean(c(-1, 2, 3))
#> [1] NaN
```
