# Calculate threshold value on ranked input

Calculate a threshold value of a vector based on the unimodal threshold
method described by [Rosin
2001](https://users.cs.cf.ac.uk/Paul.Rosin/resources/papers/unimodal2.pdf%22)
for image thresholding. Conceptually, the method attempts to draw a line
from the peak of the distribution to the tail and calculates the maximum
distance from a point to that line. The point where the distance is
maximized is the threshold.

## Usage

``` r
rank_threshold(x, show = FALSE)
```

## Arguments

- x:

  numeric vector of values.

- show_plot:

  logical. default (FALSE). Show the scatter plot for the given
  distribution.

## Value

named vector containing "x" and "y" position of calculated threshold

## Details

The function takes a numeric vector and ranks the values using
[`data.table::frank`](https://rdrr.io/pkg/data.table/man/frank.html).
with `ties.method = "first"`. The returned values for the threshold "x"
and "y", indicate the rank of the data point and the value at that rank,
respectively.

## Examples

``` r
# simulate exponential counts
x <- (1:100)^exp(3)
rank_threshold(x, show = TRUE)

#>            x            y 
#> 8.500000e+01 5.667807e+38 
```
