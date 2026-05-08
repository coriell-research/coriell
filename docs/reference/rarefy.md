# Rarefy (subsample) a matrix

This function will randomly subsample counts from rows of the input
matrix such that the colSums all have even depth.

## Usage

``` r
rarefy(x, depth)
```

## Arguments

- x:

  numeric matrix or data.frame that can be converted to a numeric
  matrix. Samples in columns features in rows.

- depth:

  desired sampling depth applied to each library.

## Examples

``` r
set.seed(123)
m <- matrix(
  sample.int(100),
  nrow = 10,
  dimnames = list(
    Gene = paste0("gene.", 1:10),
    Sample = paste0("sample.", 1:10)
  )
)

colSums(m)
#>  sample.1  sample.2  sample.3  sample.4  sample.5  sample.6  sample.7  sample.8 
#>       499       548       551       545       377       548       521       545 
#>  sample.9 sample.10 
#>       565       351 

rarefied <- rarefy(m, depth = 100)
colSums(rarefied)
#>  sample.1  sample.2  sample.3  sample.4  sample.5  sample.6  sample.7  sample.8 
#>       100       100       100       100       100       100       100       100 
#>  sample.9 sample.10 
#>       100       100 
```
