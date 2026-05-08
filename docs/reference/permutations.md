# Generate all permutations of a vector

Given the input vector, generate a matrix of permutations where each row
represents a permutation of the data. Stolen from:
<https://stackoverflow.com/a/34287541>

## Usage

``` r
permutations(x)
```

## Arguments

- x:

  vector

## Value

factorial(x) x length(x) matrix

## Examples

``` r
v <- letters[1:4]
p <- permutations(v)
p
#>       [,1] [,2] [,3] [,4]
#>  [1,] "a"  "b"  "c"  "d" 
#>  [2,] "a"  "b"  "d"  "c" 
#>  [3,] "a"  "c"  "b"  "d" 
#>  [4,] "a"  "c"  "d"  "b" 
#>  [5,] "a"  "d"  "b"  "c" 
#>  [6,] "a"  "d"  "c"  "b" 
#>  [7,] "b"  "a"  "c"  "d" 
#>  [8,] "b"  "a"  "d"  "c" 
#>  [9,] "b"  "c"  "a"  "d" 
#> [10,] "b"  "c"  "d"  "a" 
#> [11,] "b"  "d"  "a"  "c" 
#> [12,] "b"  "d"  "c"  "a" 
#> [13,] "c"  "a"  "b"  "d" 
#> [14,] "c"  "a"  "d"  "b" 
#> [15,] "c"  "b"  "a"  "d" 
#> [16,] "c"  "b"  "d"  "a" 
#> [17,] "c"  "d"  "a"  "b" 
#> [18,] "c"  "d"  "b"  "a" 
#> [19,] "d"  "a"  "b"  "c" 
#> [20,] "d"  "a"  "c"  "b" 
#> [21,] "d"  "b"  "a"  "c" 
#> [22,] "d"  "b"  "c"  "a" 
#> [23,] "d"  "c"  "a"  "b" 
#> [24,] "d"  "c"  "b"  "a" 
```
