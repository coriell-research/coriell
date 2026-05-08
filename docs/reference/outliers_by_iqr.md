# Flag outliers by IQR

Flag outliers based on the interquartile range

## Usage

``` r
outliers_by_iqr(x, threshold = 1.5, direction = c("both", "low", "high"))
```

## Arguments

- x:

  numeric vector

- threshold:

  threshold for the upper and lower limits, A.K.A. Tukey's fences

- direction:

  return TRUE if the value is above or below the outlier cutoff. default
  "both", samples above and below the threshold are called outliers.

## Value

boolean vector indicating which values of the input vector are flagged
as outliers

## Examples

``` r
x <- c(1, 1, 2, 2, 4, 6, 11)

outliers_by_iqr(x)
#> [1] FALSE FALSE FALSE FALSE FALSE FALSE  TRUE

# Using direction="low" disregards outliers above threshold, for example
outliers_by_iqr(x, direction="low")
#> [1] FALSE FALSE FALSE FALSE FALSE FALSE FALSE
```
