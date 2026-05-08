# Centered Log-ratio transformation

Calculate the centered log-ratio transformation of a vector

## Usage

``` r
clr(x, base = 2)
```

## Arguments

- x:

  numeric vector

- base:

  integer base of the log function. Default 2.

## Examples

``` r
x <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 1500)
clr(x)
#>  [1] -2.90198798 -1.90198798 -1.31702548 -0.90198798 -0.58005989 -0.31702548
#>  [7] -0.09463306  0.09801202  0.26793702  7.64875880
```
