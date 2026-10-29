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
design <- model.matrix(~group)
y <- estimateDisp(y, design)

# To perform quasi-likelihood F-tests:
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = 2)

# convert the results object to a dataframe -- do not filter the results
res_df <- edger_to_df(qlf)

head(res_df)
#>   feature_id      logFC    logCPM          F       PValue        FDR
#> 1      gene1  -1.882510 10.465423  0.6251378 0.4292381688 0.90176086
#> 2      gene2 -15.139394  8.183815 14.9523685 0.0001137826 0.05689128
#> 3      gene3  -3.444070  8.318262  1.9320795 0.1646857976 0.74669133
#> 4      gene4  -2.251411  8.561397  0.8863766 0.3465748027 0.86952748
#> 5      gene5  -1.210155 10.599488  0.2665029 0.6057444133 0.94302877
#> 6      gene6   1.224415 10.424128  0.2748761 0.6001375716 0.94302877
```
