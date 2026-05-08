# Compute TPM normalized counts

Convert counts to TPM values. This function simply performs the
computations for creating TPM normalized counts given a count matrix and
any vector of lengths. This is useful for quickly computing scaling
factors for arbitrary features and lengths. For RNA-seq, it's better to
use a transcript aligner + tximport to get TPMs. Refer to
<https://bioconductor.org/packages/release/bioc/html/tximport.html> for
more information.

## Usage

``` r
tpm(x, l, ...)

# Default S3 method
tpm(x, l, ...)

# S3 method for class 'matrix'
tpm(x, l)

# S3 method for class 'data.frame'
tpm(x, l)

# S3 method for class 'SummarizedExperiment'
tpm(x, l, assay = "counts")

# S3 method for class 'DelayedArray'
tpm(x, l)
```

## Arguments

- x:

  matrix, numeric data.frame, or SummarizedExperiment object

- l:

  vector of feature lengths

- assay:

  if SummarizedExperiment, what assay to use. Default = "counts"

## Value

matrix containing TPM normalized counts

## Examples

``` r
# Generate some fake data with fake positive feature lengths
counts <- coriell::simulate_counts()$table
gene_lengths <- rbinom(nrow(counts), 1000, 0.9)

tpms <- tpm(counts, gene_lengths)
```
