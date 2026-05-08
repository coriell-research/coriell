# Convert EdgeR results object to a data.frame

Create a data.frame from an `edgeR` results object. This function calls
[`edgeR::topTags()`](https://rdrr.io/pkg/edgeR/man/topTags.html) on the
object and extracts the `table` data.frame with all features. This
function returns all rows unsorted by default i.e.
`topTags(..., n=Inf, sort.by="none")`.

## Usage

``` r
edger_to_df(x, ...)
```

## Arguments

- x:

  `edgeR` results object to be converted

- ...:

  Additional arguments passed to
  [`edgeR::topTags()`](https://rdrr.io/pkg/edgeR/man/topTags.html)

## Value

data.frame

## Examples

``` r
library(edgeR)
#> Loading required package: limma
#> 
#> Attaching package: ‘edgeR’
#> The following object is masked from ‘package:coriell’:
#> 
#>     tpm
library(coriell)

# create some fake data
x <- data.frame(
  ctl1 = rnbinom(1000, size = 0.4, prob = 1e-5),
  ctl2 = rnbinom(1000, size = 0.4, prob = 1e-5),
  trt1 = rnbinom(1000, size = 0.4, prob = 1e-5),
  trt2 = rnbinom(1000, size = 0.4, prob = 1e-5),
  row.names = paste0("gene", 1:1000)
)

# run edger pipeline
group <- factor(c(1, 1, 2, 2))
y <- DGEList(counts = x, group = group)
y <- calcNormFactors(y)
#> calcNormFactors has been renamed to normLibSizes
design <- model.matrix(~group)
y <- estimateDisp(y, design)

# To perform quasi-likelihood F-tests:
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = 2)

# convert the results object to a dataframe -- do not filter the results
res_df <- edger_to_df(qlf)

head(res_df)
#>   feature_id        logFC    logCPM            F      PValue       FDR
#> 1      gene1  2.156023221  8.191727 7.354276e-01 0.391232269 0.8748293
#> 2      gene2 -9.220337249 10.963950 7.212490e+00 0.007299639 0.5163010
#> 3      gene3 -0.975198607 11.104462 1.600519e-01 0.689151498 0.9314002
#> 4      gene4  4.466613174  9.785699 2.758125e+00 0.096918341 0.6592418
#> 5      gene5  0.888850697  6.506156 1.377554e-01 0.710562065 0.9361819
#> 6      gene6 -0.003484469 10.391769 2.169606e-06 0.998824897 0.9988249
```
