# Linear interpolation of a value

Linear interpolation of a value

## Usage

``` r
lerp(x, min, max)
```

## Arguments

- min:

  numeric. Minimum of the desired range.

- max:

  numeric. Maximum of the desired range.

- x.:

  numeric. Normalized input value to be transformed. Value must be
  between 0-1.

## Value

Value of x within the desired range given by min, max.

## Examples

``` r
lerp(c(0.1, 0.25, 0.5, 0.75), min = 0, max = 100)
#> [1] 10 25 50 75
```
