# Flag outliers by MAD

Flag outliers based on the median absolute deviation

## Usage

``` r
outliers_by_mad(x, threshold = 3, direction = c("both", "low", "high"))
```

## Arguments

- x:

  numeric vector

- threshold:

  threshold for number of MADs used to determine outlier. default 3

- direction:

  return TRUE if the value is above or below the outlier cutoff. default
  "both", samples above and below the threshold are called outliers.

## Value

boolean vector indicating which values of the input vector are flagged
as outliers

## Examples

``` r
x <- c(1, 1, 2, 2, 4, 6, 11)

outliers_by_mad(x)
#> [1] FALSE FALSE FALSE FALSE FALSE FALSE  TRUE

# Using direction="low" disregards outliers above threshold, for example
outliers_by_mad(x, direction="low")
#> [1] FALSE FALSE FALSE FALSE FALSE FALSE FALSE
```
