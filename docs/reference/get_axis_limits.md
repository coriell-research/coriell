# Get the axis limits of a ggplot2 object

Internal helper function used to get the axis limits of a plot

## Usage

``` r
get_axis_limits(p)
```

## Arguments

- p:

  `ggplot` plot object

## Value

list with elements for x and y limits. These are the coordinate limits
if set with `coord_cartesian()`, otherwise the trained range of the
data. Limits set on scales (e.g. `xlim()`) are ignored.
