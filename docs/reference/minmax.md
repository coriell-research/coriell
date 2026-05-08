# Min-Max normalize a value within a given range

Min-Max normalize a value within a given range

## Usage

``` r
minmax(x, min = NA, max = NA)
```

## Arguments

- x:

  numeric. Value to be transformed (normalized).

- min:

  numeric. Minimum value of range. Default NA, use the min of the
  supplied vector

- max:

  numeric. Maximum value of range. Default NA, use the max of the
  supplied vector

## Value

Normalized value (0-1) of x in the range given by min, max.

## Examples

``` r
minmax(c(0, 1, 10, 100))
#> [1] 0.00 0.01 0.10 1.00

# Use a scale outside of the range of x
minmax(c(0, 1, 10, 100), min = 0, max = 1000)
#> [1] 0.000 0.001 0.010 0.100
```
