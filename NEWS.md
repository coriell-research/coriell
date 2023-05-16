## coriell 0.12.0

- Added `UMAP()` and `plot_umap()` functions. The `UMAP()` functions accepts 
PCA objects from `PCAtools`, `prcomp`, or a distance matrix or raw data matrix
and exposes the `umap.defaults` as function arguments. 
- The `plot_umap()` function provides a simple plotting method for the 
data.frame produced by the `UMAP()` function.

## coriell 0.11.0

- Internal changes to the `quickmap()` function. Avoid `pheatmap` scaling in 
favor of vectorized scaling. Speed up removeVar calculations with 
`Rfast::rowVars()` or `matrixStats::rowVars()` if available. Speed up 
clustering and distance matrix calculations by performing distance matrix
calculations with `rdist::rdist()` and clustering with `fastcluster::hclust()` 
if available. Round values in fix_extreme to better maintain original scale
limits.

## coriell 0.10.0

- Moved most packages to "Suggests" instead of "Imports" to reflect that this
package is a collection of helpers. This reduces dependencies upon install.
- Inclusion of and update for `plot_boxplot()`, `plot_density()`, and 
`plot_parallel()`. These functions are now generics that work with matrix, 
data.frame, and `SummarizedExperiment` classes

## coriell 0.9.0

- Potential breaking changes to `plot_volcano()` and `plot_md()`. For 
`plot_volcano()` set the default value for the labels to NULL and removed the
removed the labels altogether for `plot_md()`
- Added new plotting functions for expression matrices: `plot_boxplot()`, 
`plot_density()`, and `plot_parallel()`

## coriell 0.8.0

- Added a new ggplot2 theme, `theme_coriell()`
- Set defaults on `plot_volcano()` and `plot_md()` to have consistent colors
- Added new argument to `quickmap()` that removes low variance features before
plotting.

## coriell 0.7.0

- Added function for performing pairwise fisher tests relative to a reference. 
`pairwise_fisher_test()`
- Added arguments to `quickmap()` to enable fixing the colors at the extreme ends
of the data.
- Added lab_size arguments to `plot_volcano()` and `plot_md()`
- Added `rarefy()` function. Replaces `subsample_counts()`
- Removed some old functions.

## coriell 0.6.0

- Added new function `read_bismark()` that reads in a list of Bismark coverage files
and optionally filters by coverage and variance.

## coriell 0.5.0

- Updated `plot_md()`, `plot_volcano()`, and `summarize_dge()` to remove `dplyr()` 
dependency. 
**The changes to these functions are breaking**. Arguments for column names must now be quoted.
- `plot_volcano()` and `plot_md()` now support additional arguments for modifying the point
size, shape, and color. See function documentation.

## coriell 0.4.0

- Added a function for calling outlier values in columns of a numeric matrix by 
the IQR method
- Removed `process_quant_file()` function. Switched to using `tximport()` in all 
pipelines.

## coriell 0.3.0

- Added meta-analysis functions that rip off `MetaVolcanoR` but are much faster. 
  - The `meta_vote()` function implements a vote-counting strategy for determining 
  common differentially expressed genes
  - The `meta_pcombine()` function combines p-values and logFCs across studies.
  - The `plot_metavolcano()` function provides a plotting function specific to the 
  `meta_vote()` results.
- Eliminated export of `magittr` pipe. Now `coriell` doesn't export the pipe.

## coriell 0.2.0

- `edger_to_df()` returns a data.frame instead of a tibble
- New functions for generating random color palettes: `distinct_rgb_palette()`
and `random_rgb_palette()`
- Added function for defining threshold based on ranked data: `rank_threshold()`. 
Inspired by unimodal thresholding algorithm from image analysis.
- `panther_go()` now returns a `data.table` of the raw, unlisted data returned from 
the request. The original version pivoted the table wider using columns for the
GO term and the description of the GO term. This result gave inaccurate results
when using a different pathway in the function call. 
- General move towards reducing the number of dependencies in functions by either
removing outside packages and switching to base R or moving to `data.table()`.
- Added utility functions for transforming numeric values.

## coriell 0.1.0

- Updated internals for `permutation_correlation_test()` that applies permuted 
vector in a vectorized fashion over the entire matrix. 
- Removed filtering capability from `edger_to_df()` to allow for any `EdgeR` 
results objects to be used as input.
