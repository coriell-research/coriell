# Generate a distinct RGB color palette

This function uses [`kmeans()`](https://rdrr.io/r/stats/kmeans.html)
over the RGB colorspace to generate N distinct RGB colors.

## Usage

``` r
distinct_rgb_palette(n, alpha = 1, ...)
```

## Arguments

- n:

  numeric. Number of colors to generate

- alpha:

  numeric. Transparency level of the color palette (0-1). Default 1.0

- ...:

  arguments passed to [`kmeans()`](https://rdrr.io/r/stats/kmeans.html)

## Value

vector of distinct RGB colors

## Details

This function uses *very* inaccurate defaults for the
[`kmeans()`](https://rdrr.io/r/stats/kmeans.html) function in the
interest of speed. It's usually not a problem if
[`kmeans()`](https://rdrr.io/r/stats/kmeans.html) does not converge
(colors are distinct enough for most purposes). If you get warnings, or
you don't like the colors produced, you can modify the default arguments
to the [`kmeans()`](https://rdrr.io/r/stats/kmeans.html) function by
passing in additional arguments to `...`. For example, increasing
iterations can be done by passing `iter.max = 100`. Changing the kmeans
algorithm can be done by specifying `algorithm = "MacQueen"`. Changing
these arguments may eliminate warnings or produce more distinct colors.
