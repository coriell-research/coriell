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
#>   feature_id      logFC    logCPM          F       PValue        FDR
#> 1      gene1  -1.402956 10.105751  0.3630773 5.468704e-01 0.94424152
#> 2      gene2 -13.101124  6.149114 16.0052202 6.546745e-05 0.03273372
#> 3      gene3  -4.409160  9.226123  2.7029654 1.003190e-01 0.69418285
#> 4      gene4  -3.552844  9.709009  1.9407926 1.637373e-01 0.73755547
#> 5      gene5  -1.195803 10.592819  0.2623450 6.085702e-01 0.94424152
#> 6      gene6   1.236282 10.423917  0.2821887 5.953292e-01 0.94424152
```
