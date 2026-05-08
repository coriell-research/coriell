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
#>          Sample1   Sample2   Sample3    Sample4    Sample5
#> CpG1          NA 0.7403324 0.1575701         NA 0.09950194
#> CpG2 0.872853523 0.5677240 0.1141094 0.09897127 0.78451347
#> CpG3 0.672521725 0.1186733 0.1891373 0.06796681 0.21496728
#> CpG4 0.004080585        NA 0.8993360         NA         NA
#> CpG5 0.265511174 0.6003532 0.3032485 0.08691927 0.36929772

# Impute missing values with row medians
impute(X)
#>          Sample1   Sample2   Sample3    Sample4    Sample5
#> CpG1 0.157570087 0.7403324 0.1575701 0.15757009 0.09950194
#> CpG2 0.872853523 0.5677240 0.1141094 0.09897127 0.78451347
#> CpG3 0.672521725 0.1186733 0.1891373 0.06796681 0.21496728
#> CpG4 0.004080585 0.4517083 0.8993360 0.45170827 0.45170827
#> CpG5 0.265511174 0.6003532 0.3032485 0.08691927 0.36929772

# Impute missing values with row mins
impute(X, min)
#>          Sample1     Sample2   Sample3     Sample4     Sample5
#> CpG1 0.099501943 0.740332381 0.1575701 0.099501943 0.099501943
#> CpG2 0.872853523 0.567723954 0.1141094 0.098971266 0.784513465
#> CpG3 0.672521725 0.118673346 0.1891373 0.067966805 0.214967275
#> CpG4 0.004080585 0.004080585 0.8993360 0.004080585 0.004080585
#> CpG5 0.265511174 0.600353186 0.3032485 0.086919273 0.369297724
```
