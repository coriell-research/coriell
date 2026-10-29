# Perform simple imputation on rows of a matrix

Impute NA values for each row of a matrix.

## Usage

``` r
impute(x, fun = median)
```

## Arguments

- x:

  numeric matrix or data.frame that can be converted to a numeric matrix

- fun:

  Imputation function to apply to rows of the matrix. Default median

## Value

Matrix with imputed values

## Examples

``` r
# Create a matrix of values with NAs
X <- matrix(runif(25), 5, dimnames = list(paste0("CpG", 1:5), paste0("Sample", 1:5)))
X[sample.int(25, 5)] <- NA
X
#>        Sample1   Sample2    Sample3    Sample4   Sample5
#> CpG1 0.1342512 0.4513373         NA 0.04433255 0.3041900
#> CpG2 0.3221015 0.4765508         NA 0.44896594 0.5151142
#> CpG3 0.3789878 0.7516720 0.45423241 0.37770845 0.2071946
#> CpG4 0.7828940 0.5002779 0.72605308         NA 0.8537786
#> CpG5 0.9746990        NA 0.07100958         NA 0.5298018

# Impute missing values with row medians
impute(X)
#>        Sample1   Sample2    Sample3    Sample4   Sample5
#> CpG1 0.1342512 0.4513373 0.21922061 0.04433255 0.3041900
#> CpG2 0.3221015 0.4765508 0.46275838 0.44896594 0.5151142
#> CpG3 0.3789878 0.7516720 0.45423241 0.37770845 0.2071946
#> CpG4 0.7828940 0.5002779 0.72605308 0.75447356 0.8537786
#> CpG5 0.9746990 0.5298018 0.07100958 0.52980179 0.5298018

# Impute missing values with arbitrary function
impute(X, fun = function(x) { -1} )
#>        Sample1    Sample2     Sample3     Sample4   Sample5
#> CpG1 0.1342512  0.4513373 -1.00000000  0.04433255 0.3041900
#> CpG2 0.3221015  0.4765508 -1.00000000  0.44896594 0.5151142
#> CpG3 0.3789878  0.7516720  0.45423241  0.37770845 0.2071946
#> CpG4 0.7828940  0.5002779  0.72605308 -1.00000000 0.8537786
#> CpG5 0.9746990 -1.0000000  0.07100958 -1.00000000 0.5298018
```
