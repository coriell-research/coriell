# Changelog

## coriell 0.20.0

- *Breaking change*:
  [`UMAP()`](https://coriell-research.github.io/coriell/reference/UMAP.md)
  now always returns a data.frame, as documented. Previously a matrix
  was returned when no metadata was supplied. Metadata column names are
  no longer altered by `check.names`.
- *Breaking change*:
  [`plot_meth_hist()`](https://coriell-research.github.io/coriell/reference/plot_meth_hist.md)
  now uses fixed bins on \[0, 1\] (`breaks = seq(0, 1, 0.1)`) and errors
  on values outside that range.
- [`remove_var()`](https://coriell-research.github.io/coriell/reference/remove_var.md)
  now keeps `round(nrow(x) * (1 - p))` features. Floating-point
  truncation previously dropped an extra feature (e.g. 10 rows with p =
  0.8 kept 1).
- Fixed
  [`colmean()`](https://coriell-research.github.io/coriell/reference/colmean.md)
  dividing by the wrong group sizes when groups were not in factor-level
  order, when factors had unused levels, or when `na.rm = TRUE`.
- Fixed
  [`meta_de()`](https://coriell-research.github.io/coriell/reference/meta_de.md)
  fixed-effect estimates being biased towards 0 when a logFC was missing
  but its SE was present.
- Fixed exact
  [`exact_cor_test()`](https://coriell-research.github.io/coriell/reference/exact_cor_test.md)
  /
  [`permutation_correlation_test()`](https://coriell-research.github.io/coriell/reference/permutation_correlation_test.md)
  p-values when `y` contains values that are nearly, but not exactly,
  tied.
- Fixed
  [`list_to_matrix()`](https://coriell-research.github.io/coriell/reference/list_to_matrix.md)
  for numeric and factor sets,
  [`pairwise_intersections()`](https://coriell-research.github.io/coriell/reference/pairwise_intersections.md)
  with duplicated elements or a single set,
  [`rank_threshold()`](https://coriell-research.github.io/coriell/reference/rank_threshold.md)
  with tied maxima,
  [`dsc()`](https://coriell-research.github.io/coriell/reference/dsc.md)
  with NA batch labels, and
  [`rarefy()`](https://coriell-research.github.io/coriell/reference/rarefy.md)
  with single-read columns or a single-row input.
- Fixed
  [`read_bismark()`](https://coriell-research.github.io/coriell/reference/read_bismark.md)
  dropping sites due to floating-point error in `prop_samples`, and
  `return_mats = TRUE` now keeps columns in the order of `files`.
- [`plot_volcano()`](https://coriell-research.github.io/coriell/reference/plot_volcano.md)
  and
  [`plot_md()`](https://coriell-research.github.io/coriell/reference/plot_md.md)
  percent annotations are now rounded once and computed out of features
  with non-missing values, matching
  [`plot_volcano2()`](https://coriell-research.github.io/coriell/reference/plot_volcano2.md)
  and
  [`plot_md2()`](https://coriell-research.github.io/coriell/reference/plot_md2.md).
- [`plot_dist()`](https://coriell-research.github.io/coriell/reference/plot_dist.md)
  gives the same result whether or not `rdist` is installed.
- Fixed errors in
  [`plot_parallel()`](https://coriell-research.github.io/coriell/reference/plot_parallel.md)
  for SummarizedExperiment objects;
  [`plot_boxplot()`](https://coriell-research.github.io/coriell/reference/plot_boxplot.md),
  [`plot_density()`](https://coriell-research.github.io/coriell/reference/plot_density.md),
  and
  [`plot_parallel()`](https://coriell-research.github.io/coriell/reference/plot_parallel.md)
  for matrices without rownames;
  [`remove_components()`](https://coriell-research.github.io/coriell/reference/remove_components.md)
  when keeping a single PC or removing the last PC;
  [`associate_components()`](https://coriell-research.github.io/coriell/reference/associate_components.md)
  with `N = 1`;
  [`sample_n_random_cor()`](https://coriell-research.github.io/coriell/reference/sample_n_random_cor.md)
  with `n = 1`; and
  [`dfs2se()`](https://coriell-research.github.io/coriell/reference/dfs2se.md)
  with columns named “rows” or “cols”.
- [`dfs2se()`](https://coriell-research.github.io/coriell/reference/dfs2se.md)
  now errors on duplicated feature ids within an experiment instead of
  silently filling the assay with counts.
- [`edger_to_df()`](https://coriell-research.github.io/coriell/reference/edger_to_df.md)
  returns a 0-row data.frame instead of NULL when filtering leaves no
  rows.
  [`setup_chunk()`](https://coriell-research.github.io/coriell/reference/setup_chunk.md)
  invisibly returns the chunk text.
  [`clamp()`](https://coriell-research.github.io/coriell/reference/clamp.md)
  returns `numeric(0)` for empty input.
- [`plot_density2()`](https://coriell-research.github.io/coriell/reference/plot_density2.md)
  respects `show_legend`.
  [`plot_cor_pairs()`](https://coriell-research.github.io/coriell/reference/plot_cor_pairs.md)
  uses `scatter_colramp`, and `cor_digits` now sets the number of
  decimal places.
- S3 methods now accept `...` consistently with their generics.
- Documentation fixes and missing Suggests (clipr, csaw, DelayedArray,
  fgsea, KernSmooth, limma, qsmooth) added to DESCRIPTION.
- Greatly expanded the test suite to cover all exported functions.

## coriell 0.19.0

- *Breaking change*: removed unused functions `random_rgb_palette()`,
  `distinct_rgb_palette()`, `pairwise_fisher_test()`,
  `simulate_counts()`, and `env2global()`.

## coriell 0.18.0

- [`plot_volcano()`](https://coriell-research.github.io/coriell/reference/plot_volcano.md)
  and
  [`plot_md()`](https://coriell-research.github.io/coriell/reference/plot_md.md)
  get a “raster” argument. Setting raster=TRUE will allow for the point
  layers of the plots to be rasterized. this is useful when saving plots
  as vector graphics.
- *Breaking change* to
  [`panther_go()`](https://coriell-research.github.io/coriell/reference/panther_go.md).
  Results no longer contain duplicated rows for each term. The nested
  “term” column is replaced by “term_id” and “term_label” columns.
- [`panther_go()`](https://coriell-research.github.io/coriell/reference/panther_go.md)
  now raises an informative error when PANTHER cannot process the
  request (e.g. when genes in `gene_list` are missing from
  `ref_input_list`) instead of silently returning an empty table.

## coriell 0.17.0

- *Breaking change* to
  [`quickmap()`](https://coriell-research.github.io/coriell/reference/quickmap.md).
  Removed all additional processing, threshold setting, and value
  clipping except for optional variance removal. The function now simply
  relays arguments to
  [`pheatmap::pheatmap()`](https://rdrr.io/pkg/pheatmap/man/pheatmap.html).
- Make
  [`edger_to_df()`](https://coriell-research.github.io/coriell/reference/edger_to_df.md)
  return unsorted results. this is useful for returning results back in
  the same order that they are present in the original object
- Make the default theme for plotting function
  [`theme_coriell()`](https://coriell-research.github.io/coriell/reference/theme_coriell.md)
- Add function for stripping ENSEMBL version IDs from character vector

## coriell 0.16.0

- Changed
  [`theme_coriell()`](https://coriell-research.github.io/coriell/reference/theme_coriell.md)
  to drop panel borders and removed angle from x-axis text.
- Changed
  [`plot_volcano()`](https://coriell-research.github.io/coriell/reference/plot_volcano.md)
  and
  [`plot_md()`](https://coriell-research.github.io/coriell/reference/plot_md.md)
  defaults to use unaliased points, smaller text sizes for labels, and
  moved labels further to the sides.
- Added some significance testing to the
  [`pairwise_intersections()`](https://coriell-research.github.io/coriell/reference/pairwise_intersections.md).
- Added a new function for plotting pairwise correlations between
  columns in a matrix,
  [`plot_cor_pairs()`](https://coriell-research.github.io/coriell/reference/plot_cor_pairs.md)

## coriell 0.15.0

- Added an argument for
  [`panther_go()`](https://coriell-research.github.io/coriell/reference/panther_go.md)
  to include the reference gene list when performing over representation
  testing. There may still be bugs here.
- Updated site. Removed some vignettes with older workflows. TODO: add
  new vignettes with current best practices for RNAseq.
- Added new arguments for
  [`plot_volcano()`](https://coriell-research.github.io/coriell/reference/plot_volcano.md)
  and
  [`plot_md()`](https://coriell-research.github.io/coriell/reference/plot_md.md)
  to allow setting axis limits prior to determining annotation
  placement. These functions now also use
  [`theme_coriell()`](https://coriell-research.github.io/coriell/reference/theme_coriell.md)
  by default.

## coriell 0.14.0

- Updated
  [`meta_de()`](https://coriell-research.github.io/coriell/reference/meta_de.md)
  function to operate strictly on SummarizedExperiment objects. This
  function is significantly faster than the previous version.
- Added a helper function,
  [`dfs2se()`](https://coriell-research.github.io/coriell/reference/dfs2se.md)
  to convert a list of data.frames to a SummarizedExperiment object for
  use by
  [`meta_de()`](https://coriell-research.github.io/coriell/reference/meta_de.md)
- Added a helper function to perform jackknife resampling on the columns
  of a SummarizedExperiment,
  [`jackknifeSE()`](https://coriell-research.github.io/coriell/reference/jackknifeSE.md).

## coriell 0.13.0

- Removed meta-analysis functions, `meta_vote()`, `meta_pcombine()` and
  `plot_metavolcano()`, in favor of the newer
  [`meta_de()`](https://coriell-research.github.io/coriell/reference/meta_de.md)
  function which provides an interface to `metapod` for combining
  p-values in a more robust way.
- Removed a redundant scaling step in
  [`quickmap()`](https://coriell-research.github.io/coriell/reference/quickmap.md)
  when calculating breaks and made some changes to the way the
  `fix_extreme` argument behaves.
- Updated RNA-seq article with more analysis steps and helper functions
- Updated
  [`panther_go()`](https://coriell-research.github.io/coriell/reference/panther_go.md)
  to use `httr2`

## coriell 0.12.0

- Added
  [`UMAP()`](https://coriell-research.github.io/coriell/reference/UMAP.md)
  and
  [`plot_umap()`](https://coriell-research.github.io/coriell/reference/plot_umap.md)
  functions. The
  [`UMAP()`](https://coriell-research.github.io/coriell/reference/UMAP.md)
  functions accepts PCA objects from `PCAtools`, `prcomp`, or a distance
  matrix or raw data matrix and exposes the `umap.defaults` as function
  arguments.
- The
  [`plot_umap()`](https://coriell-research.github.io/coriell/reference/plot_umap.md)
  function provides a simple plotting method for the data.frame produced
  by the
  [`UMAP()`](https://coriell-research.github.io/coriell/reference/UMAP.md)
  function.

## coriell 0.11.0

- Internal changes to the
  [`quickmap()`](https://coriell-research.github.io/coriell/reference/quickmap.md)
  function. Avoid `pheatmap` scaling in favor of vectorized scaling.
  Speed up removeVar calculations with
  [`Rfast::rowVars()`](https://rdrr.io/pkg/Rfast/man/colVars.html) or
  [`matrixStats::rowVars()`](https://rdrr.io/pkg/matrixStats/man/rowVars.html)
  if available. Speed up clustering and distance matrix calculations by
  performing distance matrix calculations with
  [`rdist::rdist()`](https://rdrr.io/pkg/rdist/man/rdist.html) and
  clustering with
  [`fastcluster::hclust()`](https://rdrr.io/pkg/fastcluster/man/hclust.html)
  if available. Round values in fix_extreme to better maintain original
  scale limits.

## coriell 0.10.0

- Moved most packages to “Suggests” instead of “Imports” to reflect that
  this package is a collection of helpers. This reduces dependencies
  upon install.
- Inclusion of and update for
  [`plot_boxplot()`](https://coriell-research.github.io/coriell/reference/plot_boxplot.md),
  [`plot_density()`](https://coriell-research.github.io/coriell/reference/plot_density.md),
  and
  [`plot_parallel()`](https://coriell-research.github.io/coriell/reference/plot_parallel.md).
  These functions are now generics that work with matrix, data.frame,
  and `SummarizedExperiment` classes

## coriell 0.9.0

- Potential breaking changes to
  [`plot_volcano()`](https://coriell-research.github.io/coriell/reference/plot_volcano.md)
  and
  [`plot_md()`](https://coriell-research.github.io/coriell/reference/plot_md.md).
  For
  [`plot_volcano()`](https://coriell-research.github.io/coriell/reference/plot_volcano.md)
  set the default value for the labels to NULL and removed the removed
  the labels altogether for
  [`plot_md()`](https://coriell-research.github.io/coriell/reference/plot_md.md)
- Added new plotting functions for expression matrices:
  [`plot_boxplot()`](https://coriell-research.github.io/coriell/reference/plot_boxplot.md),
  [`plot_density()`](https://coriell-research.github.io/coriell/reference/plot_density.md),
  and
  [`plot_parallel()`](https://coriell-research.github.io/coriell/reference/plot_parallel.md)

## coriell 0.8.0

- Added a new ggplot2 theme,
  [`theme_coriell()`](https://coriell-research.github.io/coriell/reference/theme_coriell.md)
- Set defaults on
  [`plot_volcano()`](https://coriell-research.github.io/coriell/reference/plot_volcano.md)
  and
  [`plot_md()`](https://coriell-research.github.io/coriell/reference/plot_md.md)
  to have consistent colors
- Added new argument to
  [`quickmap()`](https://coriell-research.github.io/coriell/reference/quickmap.md)
  that removes low variance features before plotting.

## coriell 0.7.0

- Added function for performing pairwise fisher tests relative to a
  reference. `pairwise_fisher_test()`
- Added arguments to
  [`quickmap()`](https://coriell-research.github.io/coriell/reference/quickmap.md)
  to enable fixing the colors at the extreme ends of the data.
- Added lab_size arguments to
  [`plot_volcano()`](https://coriell-research.github.io/coriell/reference/plot_volcano.md)
  and
  [`plot_md()`](https://coriell-research.github.io/coriell/reference/plot_md.md)
- Added
  [`rarefy()`](https://coriell-research.github.io/coriell/reference/rarefy.md)
  function. Replaces `subsample_counts()`
- Removed some old functions.

## coriell 0.6.0

- Added new function
  [`read_bismark()`](https://coriell-research.github.io/coriell/reference/read_bismark.md)
  that reads in a list of Bismark coverage files and optionally filters
  by coverage and variance.

## coriell 0.5.0

- Updated
  [`plot_md()`](https://coriell-research.github.io/coriell/reference/plot_md.md),
  [`plot_volcano()`](https://coriell-research.github.io/coriell/reference/plot_volcano.md),
  and
  [`summarize_dge()`](https://coriell-research.github.io/coriell/reference/summarize_dge.md)
  to remove `dplyr()` dependency. **The changes to these functions are
  breaking**. Arguments for column names must now be quoted.
- [`plot_volcano()`](https://coriell-research.github.io/coriell/reference/plot_volcano.md)
  and
  [`plot_md()`](https://coriell-research.github.io/coriell/reference/plot_md.md)
  now support additional arguments for modifying the point size, shape,
  and color. See function documentation.

## coriell 0.4.0

- Added a function for calling outlier values in columns of a numeric
  matrix by the IQR method
- Removed `process_quant_file()` function. Switched to using
  `tximport()` in all pipelines.

## coriell 0.3.0

- Added meta-analysis functions that rip off `MetaVolcanoR` but are much
  faster.
  - The `meta_vote()` function implements a vote-counting strategy for
    determining common differentially expressed genes
  - The `meta_pcombine()` function combines p-values and logFCs across
    studies.
  - The `plot_metavolcano()` function provides a plotting function
    specific to the `meta_vote()` results.
- Eliminated export of `magittr` pipe. Now `coriell` doesn’t export the
  pipe.

## coriell 0.2.0

- [`edger_to_df()`](https://coriell-research.github.io/coriell/reference/edger_to_df.md)
  returns a data.frame instead of a tibble
- New functions for generating random color palettes:
  `distinct_rgb_palette()` and `random_rgb_palette()`
- Added function for defining threshold based on ranked data:
  [`rank_threshold()`](https://coriell-research.github.io/coriell/reference/rank_threshold.md).
  Inspired by unimodal thresholding algorithm from image analysis.
- [`panther_go()`](https://coriell-research.github.io/coriell/reference/panther_go.md)
  now returns a `data.table` of the raw, unlisted data returned from the
  request. The original version pivoted the table wider using columns
  for the GO term and the description of the GO term. This result gave
  inaccurate results when using a different pathway in the function
  call.
- General move towards reducing the number of dependencies in functions
  by either removing outside packages and switching to base R or moving
  to `data.table()`.
- Added utility functions for transforming numeric values.

## coriell 0.1.0

- Updated internals for
  [`permutation_correlation_test()`](https://coriell-research.github.io/coriell/reference/permutation_correlation_test.md)
  that applies permuted vector in a vectorized fashion over the entire
  matrix.
- Removed filtering capability from
  [`edger_to_df()`](https://coriell-research.github.io/coriell/reference/edger_to_df.md)
  to allow for any `EdgeR` results objects to be used as input.
