## coriell 0.7.0

- Added function for performing pairwise fisher tests relative to a reference. 
`pairwise_fisher_test()`
- Added arguments to `quickmap()` to enable fixing the colors at the extreme ends
of the data.
- Added lab_size arguments to `plot_volcano()` and `plot_md()`
- Removed some old functions.

## coriell 0.6.0

- Added new function `read_bismark` that reads in a list of Bismark coverage files
and optionally filters by coverage and variance.

## coriell 0.5.0

- Updated `plot_md`, `plot_volcano`, and `summarize_dge` to remove `dplyr` dependency. 
**The changes to these functions are breaking**. Arguments for column names must now be quoted.
- `plot_volcano` and `plot_md` now support additional arguments for modifying the point
size, shape, and color. See function documentation.

## coriell 0.4.0

- Added a function for calling outlier values in columns of a numeric matric by 
the IQR method
- Removed `process_quant_file` function. Switched to using `tximport` in all 
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

