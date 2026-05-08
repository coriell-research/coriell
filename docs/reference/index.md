# Package index

## Plotting

- [`plot_beads()`](https://coriell-research.github.io/coriell/reference/plot_beads.md)
  : Create a methylation bead plot
- [`plot_boxplot()`](https://coriell-research.github.io/coriell/reference/plot_boxplot.md)
  : Show boxplots for columns of data in a matrix
- [`plot_boxplot2()`](https://coriell-research.github.io/coriell/reference/plot_boxplot2.md)
  : Show boxplots for columns of data in a matrix
- [`plot_cor_pairs()`](https://coriell-research.github.io/coriell/reference/plot_cor_pairs.md)
  : Plot pairwise correlations between columns of a matrix
- [`plot_density()`](https://coriell-research.github.io/coriell/reference/plot_density.md)
  : Show density distributions for columns of data in a matrix
- [`plot_density2()`](https://coriell-research.github.io/coriell/reference/plot_density2.md)
  : Show density distributions for columns of data in a matrix
- [`plot_dist()`](https://coriell-research.github.io/coriell/reference/plot_dist.md)
  : Plot the distance between all columns of a matrix
- [`plot_enrichment()`](https://coriell-research.github.io/coriell/reference/plot_enrichment.md)
  : GSEA enrichment plot
- [`plot_md()`](https://coriell-research.github.io/coriell/reference/plot_md.md)
  : Create an MD plot from expression data
- [`plot_meth_hist()`](https://coriell-research.github.io/coriell/reference/plot_meth_hist.md)
  : Create a binned histogram of methylation percentage statistics
- [`plot_parallel()`](https://coriell-research.github.io/coriell/reference/plot_parallel.md)
  : Parallel coordinates plot of row data for each column in a matrix
- [`plot_parallel2()`](https://coriell-research.github.io/coriell/reference/plot_parallel2.md)
  : Parallel coordinates plot of row data for each column in a matrix
- [`plot_umap()`](https://coriell-research.github.io/coriell/reference/plot_umap.md)
  : Plot results of UMAP
- [`plot_volcano()`](https://coriell-research.github.io/coriell/reference/plot_volcano.md)
  : Create a volcano plot from expression data
- [`quickmap()`](https://coriell-research.github.io/coriell/reference/quickmap.md)
  : Heatmap with sensible defaults for RNA-seq expression data
- [`theme_coriell()`](https://coriell-research.github.io/coriell/reference/theme_coriell.md)
  : ggplot2 theme for coriell package
- [`random_rgb_palette()`](https://coriell-research.github.io/coriell/reference/random_rgb_palette.md)
  : Generate random RGB color palette
- [`distinct_rgb_palette()`](https://coriell-research.github.io/coriell/reference/distinct_rgb_palette.md)
  : Generate a distinct RGB color palette

## Differential expression helpers

Convenience functions useful when performing differential expression
analysis

- [`edger_to_df()`](https://coriell-research.github.io/coriell/reference/edger_to_df.md)
  : Convert EdgeR results object to a data.frame
- [`normalize_counts()`](https://coriell-research.github.io/coriell/reference/normalize_counts.md)
  : Normalize counts in a DGEList object
- [`panther_go()`](https://coriell-research.github.io/coriell/reference/panther_go.md)
  : Perform GO analysis with PANTHER
- [`simulate_counts()`](https://coriell-research.github.io/coriell/reference/simulate_counts.md)
  : Generate simulated RNA-seq data for testing purposes
- [`strip_ens()`](https://coriell-research.github.io/coriell/reference/strip_ens.md)
  : Strip version IDs from ENSEMBL identifiers
- [`summarize_dge()`](https://coriell-research.github.io/coriell/reference/summarize_dge.md)
  : Summarize RNA-seq expression results

## Meta-analysis

Functions for performing and evaluating meta-analysis on differential
expression results

- [`dfs2se()`](https://coriell-research.github.io/coriell/reference/dfs2se.md)
  : Convert a list of differential expression data.frames to a
  SummarizedExperiment
- [`meta_de()`](https://coriell-research.github.io/coriell/reference/meta_de.md)
  : Perform p-value combination for sets of differential expression
  tests
- [`jackknifeSE()`](https://coriell-research.github.io/coriell/reference/jackknifeSE.md)
  : Perform jackknife resampling on all columns of a SummarizeExperiment
  object

## Outlier detection

Functions for simple outlier detection

- [`outliers_by_iqr()`](https://coriell-research.github.io/coriell/reference/outliers_by_iqr.md)
  : Flag outliers by IQR
- [`outliers_by_mad()`](https://coriell-research.github.io/coriell/reference/outliers_by_mad.md)
  : Flag outliers by MAD
- [`outliers_by_z()`](https://coriell-research.github.io/coriell/reference/outliers_by_z.md)
  : Flag outliers by z-score
- [`rank_threshold()`](https://coriell-research.github.io/coriell/reference/rank_threshold.md)
  : Calculate threshold value on ranked input

## Dimensionality reduction

Functions to perform and examine dimensionality reduction

- [`associate_components()`](https://coriell-research.github.io/coriell/reference/associate_components.md)
  : Calculate associations between variables and principal components
- [`remove_components()`](https://coriell-research.github.io/coriell/reference/remove_components.md)
  : Remove principal components from data
- [`remove_var()`](https://coriell-research.github.io/coriell/reference/remove_var.md)
  : Remove low variance features from a matrix
- [`UMAP()`](https://coriell-research.github.io/coriell/reference/UMAP.md)
  : Perform UMAP

## Utility functions

Data transformations and general utility functions

- [`clamp()`](https://coriell-research.github.io/coriell/reference/clamp.md)
  : Limit values to a given range
- [`clr()`](https://coriell-research.github.io/coriell/reference/clr.md)
  : Centered Log-ratio transformation
- [`colmean()`](https://coriell-research.github.io/coriell/reference/colmean.md)
  : Give means of rows of matrix based on column grouping variable
- [`env2global()`](https://coriell-research.github.io/coriell/reference/env2global.md)
  : Extract variable from an environment and remove that environment
- [`geometric_mean()`](https://coriell-research.github.io/coriell/reference/geometric_mean.md)
  : Geometric mean of a vector
- [`horvath_age()`](https://coriell-research.github.io/coriell/reference/horvath_age.md)
  : Transform ages using Horvath's method
- [`impute()`](https://coriell-research.github.io/coriell/reference/impute.md)
  : Perform simple imputation on rows of a matrix
- [`lerp()`](https://coriell-research.github.io/coriell/reference/lerp.md)
  : Linear interpolation of a value
- [`list_to_matrix()`](https://coriell-research.github.io/coriell/reference/list_to_matrix.md)
  : Convert a list of vectors to a binary matrix
- [`map_value()`](https://coriell-research.github.io/coriell/reference/map_value.md)
  : Map a value in one range to a value in another
- [`minmax()`](https://coriell-research.github.io/coriell/reference/minmax.md)
  : Min-Max normalize a value within a given range
- [`pairwise_intersections()`](https://coriell-research.github.io/coriell/reference/pairwise_intersections.md)
  : Get unique pairwise intersections of a list of vectors
- [`permutations()`](https://coriell-research.github.io/coriell/reference/permutations.md)
  : Generate all permutations of a vector
- [`rarefy()`](https://coriell-research.github.io/coriell/reference/rarefy.md)
  : Rarefy (subsample) a matrix
- [`read_bismark()`](https://coriell-research.github.io/coriell/reference/read_bismark.md)
  : Read in and Filter Bismark Coverage Files
- [`sample_n_random_cor()`](https://coriell-research.github.io/coriell/reference/sample_n_random_cor.md)
  : Generate a null distribution of correlation values
- [`setup_chunk()`](https://coriell-research.github.io/coriell/reference/setup_chunk.md)
  : Setup Chunk Generator
- [`tpm()`](https://coriell-research.github.io/coriell/reference/tpm.md)
  : Compute TPM normalized counts

## Tests

Statistical tests (use with caution)

- [`pairwise_fisher_test()`](https://coriell-research.github.io/coriell/reference/pairwise_fisher_test.md)
  : Perform pairwise fisher.tests on an input matrix
- [`permutation_correlation_test()`](https://coriell-research.github.io/coriell/reference/permutation_correlation_test.md)
  : Perform a permutation correlation test on every row of a matrix

## Package data

Built-in example datasets

- [`GSE161650_de`](https://coriell-research.github.io/coriell/reference/GSE161650_de.md)
  : Differential expression results from GSE161650
- [`GSE161650_lc`](https://coriell-research.github.io/coriell/reference/GSE161650_lc.md)
  : Normalized log2 counts from GSE161650
