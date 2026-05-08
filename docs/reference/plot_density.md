# Show density distributions for columns of data in a matrix

Create a density plot of values for columns in the given expression
matrix.

## Usage

``` r
plot_density(x, ...)

# Default S3 method
plot_density(x)

# S3 method for class 'matrix'
plot_density(x, metadata = NULL, colBy = NULL, ...)

# S3 method for class 'data.frame'
plot_density(x, metadata = NULL, colBy = NULL, ...)

# S3 method for class 'SummarizedExperiment'
plot_density(x, assay = "counts", colBy = NULL, ...)
```

## Arguments

- x:

  gene by sample matrix or `SummarizedExperiment` object

- ...:

  Additional parameters passed to
  [`ggplot2::geom_density()`](https://ggplot2.tidyverse.org/reference/geom_density.html)

- metadata:

  data.frame containing metadata per sample. rownames of metadata must
  match the colnames of the input matrix. Default NULL, each sample in
  the matrix will be plotted.

- colBy:

  metadata column used to color density lines. Default NULL, each sample
  in the matrix will be plotted.

- assay:

  assay of the SummarizedExperiment object. Default "counts"

## Value

ggplot object

## Examples

``` r
# Create metadata for plotting
metadata <- data.frame(row.names = colnames(GSE161650_lc))
metadata$Group <- rep(c("DMSO", "THZ1"), each = 3)

# Plot the density by sample
plot_density(GSE161650_lc)


# Color each sample by their Group in metadata
plot_density(GSE161650_lc, metadata, colBy = "Group")

```
