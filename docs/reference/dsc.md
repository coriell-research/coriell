# Calculate the Dispersion Separability Criterion (DSC)

This function computes the [Dispersion Separability Criterion
(DSC)](https://ibl.mdanderson.org/public-software/tcga-batch-effects/)
metric which was defined by the MD Anderson Bioinformatics team for
computing the strength of batch effects in high-throughput sequencing
data from TCGA.

## Usage

``` r
dsc(x, g)
```

## Arguments

- x:

  A numeric matrix (features in rows, samples in columns)

- g:

  A vector containing the batch assignment for each sample

## Value

DSC metric

## Examples

``` r
set.seed(123)

# Small DSC when samples drawn from same distribution
x1 <- matrix(rnorm(500 * 100), nrow = 500, ncol = 100)
batches <- rep(c("Batch_1", "Batch_2"), each = 50)
dsc(x = x1, g = batches)
#> [1] 0.09999285

# Large DSC when samples have batch effect present
x2 <- x1
x2[, 51:100] <- x2[, 51:100] + 50
dsc(x2, batches)
#> [1] 25.24697
```
