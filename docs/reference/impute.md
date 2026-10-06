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
#>         Sample1   Sample2   Sample3   Sample4    Sample5
#> CpG1 0.08410487 0.4434200        NA        NA 0.07100958
#> CpG2         NA 0.1342512 0.4513373 0.8447113 0.04433255
#> CpG3 0.20129167 0.3221015 0.4765508        NA 0.44896594
#> CpG4 0.15412326 0.3789878 0.7516720        NA 0.37770845
#> CpG5 0.29067720 0.7828940 0.5002779 0.7260531 0.84237378

# Impute missing values with row medians
impute(X)
#>         Sample1   Sample2    Sample3    Sample4    Sample5
#> CpG1 0.08410487 0.4434200 0.08410487 0.08410487 0.07100958
#> CpG2 0.29279427 0.1342512 0.45133730 0.84471128 0.04433255
#> CpG3 0.20129167 0.3221015 0.47655081 0.38553374 0.44896594
#> CpG4 0.15412326 0.3789878 0.75167202 0.37834812 0.37770845
#> CpG5 0.29067720 0.7828940 0.50027788 0.72605308 0.84237378

# Impute missing values with arbitrary function
impute(X, fun = function(x) { -1} )
#>          Sample1   Sample2    Sample3    Sample4    Sample5
#> CpG1  0.08410487 0.4434200 -1.0000000 -1.0000000 0.07100958
#> CpG2 -1.00000000 0.1342512  0.4513373  0.8447113 0.04433255
#> CpG3  0.20129167 0.3221015  0.4765508 -1.0000000 0.44896594
#> CpG4  0.15412326 0.3789878  0.7516720 -1.0000000 0.37770845
#> CpG5  0.29067720 0.7828940  0.5002779  0.7260531 0.84237378
```
