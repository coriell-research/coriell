# Remove principal components from data

Reconstruct data after removing principal components. This function will
reconstruct a truncated version of the original data matrix after
removing the specified principal components. The source code for this
function was stolen from [stack
exchange](https://stats.stackexchange.com/questions/57467/how-to-perform-dimensionality-reduction-with-pca-in-r/57478#57478)

## Usage

``` r
remove_components(x, components = 1, ...)
```

## Arguments

- x:

  data.frame or matrix of original data.

- components:

  numeric vector of components to remove from original data. Default (1)

- ...:

  Additional arguments passed to
  [`prcomp()`](https://rdrr.io/r/stats/prcomp.html)

## Value

data.frame of values in the original units after removing components

## Examples

``` r
# Remove first two components from dataset
trunc <- remove_components(USArrests, components = 1:2, scale. = TRUE, center = TRUE)

# original data
head(USArrests)
#>            Murder Assault UrbanPop Rape
#> Alabama      13.2     236       58 21.2
#> Alaska       10.0     263       48 44.5
#> Arizona       8.1     294       80 31.0
#> Arkansas      8.8     190       50 19.5
#> California    9.0     276       91 40.6
#> Colorado      7.9     204       78 38.7

# reconstructed data
head(trunc)
#>              Murder  Assault UrbanPop     Rape
#> Alabama    8.879093 171.0042 68.24625 17.99226
#> Alaska     3.558807 152.5293 53.64856 36.33858
#> Arizona    5.370958 220.7384 63.64209 20.95841
#> Arkansas   7.107685 179.4374 64.56869 21.94987
#> California 5.949990 178.4936 61.64172 25.48834
#> Colorado   6.181011 146.4461 59.61149 29.53625
```
