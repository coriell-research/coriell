# Get unique pairwise intersections of a list of vectors

This function takes in a list of vectors and performs pairwise set
intersections for all unique pairs of vectors in the list.

## Usage

``` r
pairwise_intersections(x, universe_size = NULL)
```

## Arguments

- x:

  List of vectors to perform intersections on

- universe_size:

  Size of the universe of features each set was drawn from. default
  NULL. If supplied then a fisher's exact test is performed to assess
  the significance of the overlap. Note, the same universe size is used
  for all pairwise comparisons.

## Value

data.table

a data.table containing columns for the sets being compared, a list
column which contains the actual values in the intersection, a column
with the intersection size, and a column with the Jaccard index, and
optionally a P-Value column containing the p-value from a fisher's exact
test if universe_size is given.

## Examples

``` r
l <- list(
  Set1 = c("A", "B", "C"),
  Set2 = c("B", "C", "D"),
  Set3 = c("X", "Y", "Z"),
  Set4 = LETTERS
)

pairwise_intersections(l)
#>      Set1   Set2 Set1Size Set2Size IntersectionSize UnionSize   Jaccard
#>    <char> <char>    <num>    <num>            <num>     <num>     <num>
#> 1:   Set1   Set2        3        3                2         4 0.5000000
#> 2:   Set1   Set3        3        3                0         6 0.0000000
#> 3:   Set1   Set4        3       26                3        26 0.1153846
#> 4:   Set2   Set3        3        3                0         6 0.0000000
#> 5:   Set2   Set4        3       26                3        26 0.1153846
#> 6:   Set3   Set4        3       26                3        26 0.1153846
#>    Intersection
#>          <list>
#> 1:          B,C
#> 2:             
#> 3:        A,B,C
#> 4:             
#> 5:        B,C,D
#> 6:        X,Y,Z
```
