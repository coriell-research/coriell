## coriell 0.17.0

- *Breaking change* to `quickmap()`. Removed all additional processing, 
threshold setting, and value clipping except for optional variance removal. 
The function now simply relays arguments to `pheatmap::pheatmap()`. 
- Make `edger_to_df()` return unsorted results. this is useful for returning 
results back in the same order that they are present in the original object
- Make the default theme for plotting function `theme_coriell()`
- Add function for stripping ENSEMBL version IDs from character vector

## coriell 0.16.0

- Changed `theme_coriell()` to drop panel borders and removed angle from x-axis 
text. 
- Changed `plot_volcano()` and `plot_md()` defaults to use unaliased points, 
smaller text sizes for labels, and moved labels further to the sides. 
- Added some significance testing to the `pairwise_intersections()`. 
- Added a new function for plotting pairwise correlations between columns in a 
matrix, `plot_cor_pairs()`

## coriell 0.15.0

- Added an argument for `panther_go()` to include the reference gene list when performing over representation testing. There may still be 
bugs here.
- Updated site. Removed some vignettes with older workflows. TODO: add
new vignettes with current best practices for RNAseq.
- Added new arguments for `plot_volcano()` and `plot_md()` to allow setting axis
limits prior to determining annotation placement. These functions now also use
`theme_coriell()` by default.

## coriell 0.14.0

- Updated `meta_de()` function to operate strictly on SummarizedExperiment objects. This function is significantly faster than the previous version.
- Added a helper function, `dfs2se()` to convert a list of data.frames to a SummarizedExperiment object for use by `meta_de()`
- Added a helper function to perform jackknife resampling on the columns of a SummarizedExperiment, `jackknifeSE()`.

## coriell 0.13.0

- Removed meta-analysis functions, `meta_vote()`, `meta_pcombine()` and 
`plot_metavolcano()`, in favor of the newer `meta_de()` function which
provides an interface to `metapod` for combining p-values in a more robust way.
- Removed a redundant scaling step in `quickmap()` when calculating breaks and
made some changes to the way the `fix_extreme` argument behaves.
- Updated RNA-seq article with more analysis steps and helper functions
- Updated `panther_go()` to use `httr2`

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
