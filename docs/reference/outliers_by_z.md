# Flag outliers by z-score

Flag outliers that exceed a given z-score value.

## Usage

``` r
outliers_by_z(x, threshold = 3, direction = c("both", "low", "high"))
```

## Arguments

- x:

  numeric vector

- threshold:

  threshold for number of SDs used to determine outlier. default 3

- direction:

  return TRUE if the value is above or below the outlier cutoff. default
  "both", samples above and below the threshold are called outliers.

## Value

boolean vector indicating which values of the input vector are flagged
as outliers

## Examples

``` r
x <- c(rnorm(10), 100)

outliers_by_z(x)
#>  [1] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE  TRUE

# Using direction="low" disregards outliers above threshold, for example
outliers_by_z(x, direction="low")
#>  [1] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
```
