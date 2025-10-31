# coriell

This package contains helper functions for common bioinformatics tasks (and some 
not-so-common tasks). For additional notes and use cases, check out the 
[Coriell bioinformatics notes](https://coriell-research.github.io/coriell-bioinformatics-notes/)
Bookdown book. 

## Installation

The latest version can be installed from Github using:

```R
if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("coriell-research/coriell")
```

## Examples

### Plotting functions

- [Volcano plot from differential expression results](https://coriell-research.github.io/coriell/reference/plot_volcano.html)
- [MA plot from differential expression results](https://coriell-research.github.io/coriell/reference/plot_md.html)
- [Heatmap with sensible defaults for RNA-seq](https://coriell-research.github.io/coriell/reference/quickmap.html)
- [RLE boxplots of an expression matrix](https://coriell-research.github.io/coriell/reference/plot_boxplot.html)
- [Density plot of an expression matrix](https://coriell-research.github.io/coriell/reference/plot_density.html)
- [Parallel coordinates plot of expression matrix](https://coriell-research.github.io/coriell/reference/plot_parallel.html)
- [Plot pairwise correlations between all samples in a matrix](https://coriell-research.github.io/coriell/reference/plot_cor_pairs.html)
- [Plot heatmap of pairwise distances between all samples in a matrix](https://coriell-research.github.io/coriell/reference/plot_dist.html)
- [Plot GSEA enrichment plots](https://coriell-research.github.io/coriell/reference/plot_enrichment.html)
- [Plot histogram of methylation beta-values](https://coriell-research.github.io/coriell/reference/plot_meth_hist.html)

### Differential Expression and Meta-Analysis

- [Perform Meta-Analysis across sets of differential expression results](https://coriell-research.github.io/coriell/reference/meta_de.html)
- [Perform jackknife resampling on meta-analysis results to assess stability](https://coriell-research.github.io/coriell/reference/jackknifeSE.html)

### Dimensionality reduction

- [Calculate associations between variables and principal components](https://coriell-research.github.io/coriell/reference/associate_components.html)
- [Remove the effect of principal components from data](https://coriell-research.github.io/coriell/reference/remove_components.html)
- [Perform UMAP on various objects](https://coriell-research.github.io/coriell/reference/UMAP.html)

### Miscellaneous 

The package also contains many other convenience functions so be sure to check 
out the [reference](https://coriell-research.github.io/coriell/reference/index.html) 
page and articles as well.

## Built-in datasets

These datasets are built into the package for testing purposes and are used to 
illustrate functionality in reference pages.

* `GSE161650_de` : Differential expression results of THZ1 vs DMSO from GSE161650 
* `GSE161650_lc` : Normalized log2 counts from THZ1 vs DMSO replicates from GSE161650 

The head of `GSE161650_de` looks like:

```
>   feature_id     logFC unshrunk.logFC    logCPM       PValue          FDR
> 1        JUN  5.759233       5.759908  9.079350 2.919097e-14 1.990666e-10
> 2       IER5  3.931325       3.931420 10.158336 3.365738e-14 1.990666e-10
> 3    GADD45B  5.813030       5.814071  8.432432 6.435666e-14 2.537583e-10
> 4       IER2  4.457981       4.458016 12.223835 1.528890e-13 4.521309e-10
> 5     PIK3R3 -4.018325      -4.019603  7.124408 2.122484e-13 4.752097e-10
> 6     HEXIM1  4.497345       4.497844  8.275561 2.696985e-13 4.752097e-10
```

And `GSE161650_lc`:

```
>          DMSO.1     DMSO.2     DMSO.3    THZ1.1    THZ1.2     THZ1.3
> A1BG  5.3323512  5.4576081  5.2876011  6.703752 6.8090471  6.7908595
> AAAS  3.8738768  3.8839857  3.5242625  3.768811 4.1406003  3.7454773
> AACS  2.1381539  2.3748462  2.2761971  3.769163 3.4352003  3.4064114
> AADAT 0.8240013 -0.2391457 -0.8699138 -1.342273 0.1967574 -0.8355924
> AAED1 1.2814586  1.5332476  2.0735095  2.443188 2.1270932  0.9707206
> AAGAB 6.8747238  6.7396922  6.6670762  6.756664 6.5840839  6.7662317
```

See the package documentation `?GSE161650_de` and `?GSE161650_lc` for 
citation information.
