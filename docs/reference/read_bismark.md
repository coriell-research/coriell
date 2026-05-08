# Read in and Filter Bismark Coverage Files

This function will read in Bismark coverage files and optionally filter
for coverage and feature-wise variance.

## Usage

``` r
read_bismark(
  files,
  coverage = 1,
  prop_samples = 1,
  remove_zero_var = FALSE,
  return_mats = FALSE
)
```

## Arguments

- files:

  Named vector of file paths to the bismark coverage files

- coverage:

  Minimum coverage for a CpG site. Default (1)

- prop_samples:

  Proportion of samples that site must be covered to be retained.
  Default (1)

- remove_zero_var:

  Should zero variance features be removed after filtering for coverage?
  Default (FALSE)

- return_mats:

  if TRUE, return matrices of the filtered data in a list. Else, return
  the filtered data.table in bismark format. Default (FALSE)

## Value

Either a data.table of filtered data with additional columns for
Coverage and Variance or a list of matrices of filtered data

## Details

The function can optionally return a list of CpG x Sample matrices of
the filtered data if return_mats = TRUE. The following matrices will be
returned for each of the retained CpG sites: Coverage, Percent
Methylation (Percent), count of methylated CpGs (Methylated) and count
of unmethylated CpGs (Unmethylated)

## Examples

``` r
if (FALSE) { # \dontrun{
files <- c("P6_1.bismark.cov.gz", "P6_4.bismark.cov.gz", "P7_2.bismark.cov.gz", "P7_5.bismark.cov.gz", "P8_3.bismark.cov.gz", "P8_6.bismark.cov.gz")
names(files) <- gsub("\\.bismark\\.cov\\.gz", "", files)

# Read in all files using default parameters -- returns a data.table
result <- read_bismark(files)

# Filter for coverage >= 20 in at least 50% of samples and remove any zero variance features
dt <- read_bismark(files, coverage = 20, prop_samples = 0.5, remove_zero_var = TRUE)

# Same filtering as above but return a list of CpG x Sample matrices
l <- read_bismark(files, coverage = 20, prop_samples = 0.5, remove_zero_var = TRUE, return_mats = TRUE)

# To view the CpG x Sample Coverages:
l$Coverage

# Percent methylation
l$Percent
} # }
```
