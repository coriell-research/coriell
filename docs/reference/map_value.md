# Map a value in one range to a value in another

Map a value in one range to a value in another

## Usage

``` r
map_value(x, old_min, old_max, new_min, new_max)
```

## Arguments

- x:

  numeric. Value to be transformed

- old_min:

  numeric. Minimum value of source range.

- old_max:

  numeric. Maximum value of source range.

- new_min:

  numeric. Minimum value of desired range.

- new_max:

  numeric. Maximum value of desired range.

## Value

Value of x mapped to the desired range.

## Examples

``` r
# Map values from the range 0-10 to the range 0-1
map_value(c(0, 1, 2, 3), 0, 10, 0, 1)
#> [1] 0.0 0.1 0.2 0.3
```
