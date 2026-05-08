# Calculate associations between variables and principal components

Calculate associations between metadata variables and PCA rotations.
This function is inspired by `methylKit::assocComp` but is designed to
work with arbitrary input data.

## Usage

``` r
associate_components(x, metadata, N = 10, ...)
```

## Arguments

- x:

  data.frame or matrix of values passed to
  [`prcomp()`](https://rdrr.io/r/stats/prcomp.html)

- metadata:

  data.frame of metadata variables to test associations.
  rownames(metadata) must be identical to colnames(x).

- N:

  First N PCs to test associations between. Default 10.

- ...:

  Additional arguments passed to the
  [`prcomp()`](https://rdrr.io/r/stats/prcomp.html) function.

## Value

data.frame with rows for each metadata variable, columns for each PC,
and p-values from the given test in cells.

## Details

This function returns the p-values from testing the association between
a given metadata column and all PC rotations. For numeric values the
p-value returned is computed using the
[`cor.test()`](https://rdrr.io/r/stats/cor.test.html) function. For
factor variables and variables that can be converted to factor variables
the function will return p-values from the
[`wilcox.test()`](https://rdrr.io/r/stats/wilcox.test.html) function or
the [`kruskal.test()`](https://rdrr.io/r/stats/kruskal.test.html)
function (when the metadata variable has \> 2 levels).

## Examples

``` r
# Specify metadata
metadata <- data.frame(
  age = c(30, 80, 34, 30, 80, 40),
  treatment = factor(c(rep("Treatment", 3), rep("Control", 3))),
  class = factor(c(rep("A", 2), rep("B", 2), rep("C", 2))),
  row.names = c(paste0("trt", 1:3), paste0("ctrl", 1:3))
)

# Create values to perform PCA on
df <- data.frame(replicate(6, runif(1000, 0, 100)))
colnames(df) <- c(paste0("trt", 1:3), paste0("ctrl", 1:3))

# Test for associations
res <- associate_components(df, metadata)

# Show results
head(res)
#>                 PC1       PC2       PC3       PC4       PC5       PC6
#> age       0.3669306 0.6047674 0.2282258 0.6979298 0.4756648 0.1893065
#> treatment 1.0000000 0.7000000 1.0000000 0.4000000 0.4000000 0.1000000
#> class     0.5647181 0.3678794 0.5647181 1.0000000 0.1017014 0.1800923
```
