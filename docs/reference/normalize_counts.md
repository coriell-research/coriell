# Normalize counts in a DGEList object

This function wraps additional normalization methods (qsmooth and cyclic
loess) for a given DGEList object in addition to the standard
normalization methods implemented by `edgeR::normlibSizes()`.

## Usage

``` r
normalize_counts(
  x,
  method = c("TMM", "TMMwsp", "RLE", "upperquartile", "none", "qsmooth", "loess"),
  refColumn = NULL,
  logratioTrim = 0.3,
  sumTrim = 0.05,
  doWeighting = TRUE,
  Acutoff = -1e+10,
  p = 0.75,
  group_factor = NULL,
  batch = NULL,
  norm_factors = NULL,
  window = 0.05,
  prior_count = 0.1,
  weights = NULL,
  span = 0.3,
  iterations = 4L,
  min.weight = 1e-05,
  max.weight = 1e+05,
  equal.weights.as.null = TRUE,
  loess_method = c("weightedLowess", "loess", "locfit")
)
```

## Arguments

- x:

  DGEList object

- method:

  Normalization method. One of "TMM", "TMMwsp", "RLE", "upperquartile",
  "none", "qsmooth", or "loess".

- refColumn:

  column to use as reference for method="TMM". Can be a column number or
  a numeric vector of length nrow(object).

- logratioTrim:

  the fraction (0 to 0.5) of observations to be trimmed from each tail
  of the distribution of log-ratios (M-values) before computing the
  mean. Used by method="TMM" for each pair of samples.

- doWeighting:

  logical, whether to use (asymptotic binomial precision) weights when
  computing the mean M-values. Used by method="TMM" for each pair of
  samples.

- Acutoff:

  minimum cutoff applied to A-values. Count pairs with lower A-values
  are ignored. Used by method="TMM" for each pair of samples.

- p:

  numeric value between 0 and 1 specifying which quantile of the counts
  should be used by method="upperquartile".

- group_factor:

  a group level continuous or categorial covariate associated with each
  sample or column in the object. The order of the group_factor must
  match the order of the columns in object. Used in qsmooth
  normalization

- batch:

  (Optional) batch covariate (multiple batches are not allowed). If
  batch covariate is provided, `Combat()` from sva is used prior to
  qsmooth normalization to remove batch effects. See `Combat()` for more
  details. Used in qsmooth.

- norm_factors:

  optional normalization scaling factors. Used in qsmooth.

- window:

  window size for running median which is a fraction of the number of
  rows in object. Default is 0.05. Used in qsmooth.

- prior_count:

  prior count to add to 0 counts when computing offsets. Default = 0.1.
  Used in qsmooth.

- weights:

  numeric vector of non-negative prior weights. Missing values are
  treated as zero. Default = NULL, Used in loess.

- span:

  positive numeric value between 0 and 1 specifying proportion of data
  to be used in the local regression moving window. Larger numbers give
  smoother fits. Default = 0.3. Used in loess.

- iterations:

  number of local regression fits. Values greater than 1 produce robust
  fits. Default = 4L. Used in loess.

- min.weight:

  minimum weight. Any lower weights will be reset. Default = 1e-5. Used
  in loess.

- max.weight:

  maximum weight. Any higher weights will be reset. Default = 1e5. Used
  in loess.

- equal.weights.as.null:

  should equal weights be treated as if weights were NULL, so that
  lowess is called? Applies even if all weights are all zero. Default =
  TRUE. Used in loess.

- loess_method:

  method used for weighted lowess. Possibilities are "weightedLowess",
  "loess" or "locfit". Used in loess.

- sumTrimthe:

  fraction (0 to 0.5) of observations to be trimmed from each tail of
  the distribution of A-values before computing the mean. Used by
  method="TMM" for each pair of samples.

## Value

DGEList with norm.factors or offset matrix

## Details

If "TMM", "TMMwsp", "RLE", "upperquartile", or "none" is selected, then
`edgeR::normLibSizes(x, ...)` is used. If method = "qsmooth" then
`qsmooth::qsmooth(x, group_factor, ...)` is used. If method = "loess"
then `csaw::normOffsets(x, ...)` is used.

## Examples

``` r
counts <- matrix(rnbinom(1e3, mu = 10, size = 20), ncol = 20)

y <- edgeR::DGEList(
  counts = counts,
  group = gl(n = 2, k = 10, labels = c("control", "treatment"))
)

# TMM normalization
y.tmm <- normalize_counts(y, method = "TMM")

# Qsmooth normalization
y.qs <- normalize_counts(y, method = "qsmooth", group_factor = y$samples$group)

# Cyclic loess
y.loess <- normalize_counts(y, method = "loess")
#> Error in normalize_counts(y, method = "loess"): csaw package is required for method='loess'
```
