# Convert a list of vectors to a binary matrix

Convert a list of sets into a binary matrix.

## Usage

``` r
list_to_matrix(sets)
```

## Arguments

- sets:

  named list of sets (vectors) to be converted to binary matrix.

## Value

binary matrix with a column for each set. rownames of the matrix
represent the union of all sets. A '1' indicates the inclusion of the
element in the set.

## Examples

``` r
sets <- list(
  "set1" = letters[1:5],
  "set2" = letters[2:6],
  "set3" = letters[1:7]
)

list_to_matrix(sets)
#>   set1 set2 set3
#> a    1    0    1
#> b    1    1    1
#> c    1    1    1
#> d    1    1    1
#> e    1    1    1
#> f    0    1    1
#> g    0    0    1
```
