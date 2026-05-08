# Generate random RGB color palette

This function will randomly generate N colors from the rgb colorspace.
It makes no attempt at ensuring the colors are distinct or well
separated.

## Usage

``` r
random_rgb_palette(n, alpha = 1)
```

## Arguments

- n:

  numeric. Number of random colors to generate

- alpha:

  numeric. Transparency level of color palette (0-1). Default 1.0

## Value

vector of hex values for RGB colors of length n
