# Remove low variance features from a matrix

This function removes the p lowest variance features from a matrix. The
function expects features in rows and samples in columns (e.g. an
expression matrix).

## Usage

``` r
remove_var(x, ...)

# Default S3 method
remove_var(x)

# S3 method for class 'matrix'
remove_var(x, p)

# S3 method for class 'data.frame'
remove_var(x, p)

# S3 method for class 'SummarizedExperiment'
remove_var(x, p, assay = "counts")

# S3 method for class 'DelayedArray'
remove_var(x, p)
```

## Arguments

- x:

  matrix, numeric data.frame, or SummarizedExperiment object

- p:

  proportion of low variance features to remove

- assay:

  if SummarizedExperiment, what assay to use. Default = "counts"

## Value

matrix with lowest variance features removed

## Examples

``` r
# Remove 80% lowest variance features
removed <- remove_var(GSE161650_lc, p = 0.8)
#> Removing 80% lowest variance features...
```
