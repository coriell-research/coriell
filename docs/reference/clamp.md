# Limit values to a given range

This function is useful for input validation. If the value of x is in
the given range then the function returns x. If the value is outside of
the range then the function returns either the max or the min value of
the range.

## Usage

``` r
clamp(x, min, max)
```

## Arguments

- x:

  numeric. Value to validate

- min:

  numeric. Minimum value of range.

- max:

  numeric. Maximum value of range.

## Value

The clamped value.

## Examples

``` r
clamp(c(-1, 121, 10, 15), min = 0, max = 100)
#> [1]   0 100  10  15
```
