## coriell 0.2.0

- `edger_to_df()` returns a data.frame instead of a tibble
- New functions for generating random color palettes: `distinct_rgb_palette()`
and `random_rgb_palette()`
- Function for defining threshold based on ranked data: `unimodal_threshold()`. 
Inspired by unimodal thresholding algorithm from image analysis. Function name
might change in later version to something more descriptive (i.e. `rank_threshold()`)
- General move towards reducing the number of dependencies in functions by either
removing outside packages and switching to base R or moving to `data.table()`. 

## coriell 0.1.0

- Updated internals for `permutation_correlation_test()` that applies permuted 
vector in a vectorized fashion over the entire matrix. 
- Removed filtering capability from `edger_to_df()` to allow for any `EdgeR` 
results objects to be used as input.

