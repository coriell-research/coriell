# Give means of rows of matrix based on column grouping variable

Inspired by
[`edgeR::sumTechReps`](https://rdrr.io/pkg/edgeR/man/sumTechReps.html)
and [`base::rowsum()`](https://rdrr.io/r/base/rowsum.html), this
function takes the average of the values in each group given by the
group argument for each row of the data matrix.

## Usage

``` r
colmean(x, group, na.rm = FALSE)
```

## Arguments

- x:

  variable (gene) by sample numeric matrix

- group:

  Factor specifying the grouping level to by averaged

- na.rm:

  Logical (TRUE or FALSE). Should NA (including NaN) values be
  discarded?

## Value

variable x nLevels(group) matrix

## Examples

``` r
# Specify the Group levels
Group <- gl(n = 2, k = 3, labels = c("DMSO", "THZ1"))

# Take the average of every gene by treatment group
by_group <- colmean(GSE161650_lc, group = Group)
by_group[1:5, ]
#>             DMSO       THZ1
#> A1BG   5.3591868  6.7678862
#> AAAS   3.7607083  3.8849627
#> AACS   2.2630657  3.5369248
#> AADAT -0.0950194 -0.6603692
#> AAED1  1.6294052  1.8470007
```
