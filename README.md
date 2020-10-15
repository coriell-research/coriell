# coriell

Helper functions for common bioinformatics tasks. If you find yourself reusing old code over and over, let me know and we'll 
create a function and put it in this package. 

I will be adding tests, however, **please test functions and let me know when things are broken** by either raising an issue on Github or contacting me.

## Installation

Since this package will be constantly changing be sure to install the latest version from Github using:

```R
# Make sure you have devtools installed
install.packages("devtools")

# Then install using devtools::install_github
devtools::install_github("coriell-research/coriell")
```

## Examples

- [Correlate methylation data with age](https://github.com/coriell-research/coriell#perform-correlation-permutation-test-using-multiple-cores)
- [Create a null distribution of correlations](https://github.com/coriell-research/coriell#extract-random-correlations-from-a-dataframe)
- [Convert edgeR results into tidy dataframes](https://github.com/coriell-research/coriell#convert-edger-results-into-a-tidy-dataframe)
- [Summarize results from differential expression analysis](https://github.com/coriell-research/coriell#summarize-results-from-differential-expression-test)
- [Create volcano plot from differential expression results](https://github.com/coriell-research/coriell#create-volcano-plot-from-differential-expression-results)
- [Create md plot from differential expression results](https://github.com/coriell-research/coriell#create-md-plot-from-differential-expression-results)
- [Heatmap with sensible defaults](https://github.com/coriell-research/coriell#heatmap-with-sensible-defaults)
- [Z-score a dataframe](https://github.com/coriell-research/coriell#z-score-a-dataframe)
- [Convert a list of sets into a binary matrix](https://github.com/coriell-research/coriell#convert-a-list-of-sets-into-a-binary-matrix)
- [Get statistics for all pairwise combinations of a list of sets](https://github.com/coriell-research/coriell#get-statistics-for-all-pairwise-combinations-of-a-list-of-sets)

### Perform correlation permutation test using multiple cores

Note: if the number of possible permutations of your data is less than the desired number
of permutations given by the parameter `n_perm` then an exact test on all of the real permutations
will be performed instead of random sampling. For example, if you only have 6 samples (720 possible permutations)
but set `n_perm = 1000` only 720 permutations (i.e. the exact test) will be tested.
The `permutation_correlation_test` function will display a message if this occurs.

```R
library(coriell)
library(methylKit)

# define age dataframe -- ages matching column order
ages = data.frame(age = c(30, 80, 34, 30, 80, 40, 35, 80))

# simulate a methylation dataset
sim_meth <- dataSim(
  replicates = 8,
  sites = 1000,
  treatment = c(rep(1, 4), rep(0, 4)),
  covariates = ages,
  sample.ids = c(paste0("test", 1:4), paste0("ctrl", 1:4))
)

# extract the methylation as percentages and coerce to data.frame
perc_meth <- as.data.frame(percMethylation(sim_meth))

head(perc_meth)
>       test1     test2     test3     test4    ctrl1     ctrl2     ctrl3      ctrl4
> 1  48.14815  54.05405  45.45455  57.89474 32.46753 24.418605  4.687500  30.769231
> 2   0.00000   0.00000   0.00000   0.00000  0.00000  0.000000  0.000000   0.000000
> 3  38.20598  16.34615  12.50000  60.00000 25.42373  3.636364 18.421053  32.335329
> 4 100.00000 100.00000 100.00000 100.00000 86.15385 72.340426 86.666667 100.000000
> 5   0.00000   0.00000   0.00000   0.00000  0.00000  0.000000  0.000000  10.526316
> 6  26.98413  12.63158  23.07692  10.00000 22.22222  0.000000  3.174603   4.166667

# permutation testing -----------------------------------------------------------
# perform permutation testing using 4 cores and 10000 permutations
res <- permutation_correlation_test(perc_meth, 
                                    y = ages$age, 
                                    n_cores = 4, 
                                    n_perm = 10000,
                                    cor_method = "spearman", 
                                    p_adjust_method = "fdr")

head(res)
>       test1     test2     test3     test4    ctrl1     ctrl2     ctrl3      ctrl4   spearman empirical_p       fdr
> 1  48.14815  54.05405  45.45455  57.89474 32.46753 24.418605  4.687500  30.769231 -0.3437200      0.1925 0.4812073
> 2   0.00000   0.00000   0.00000   0.00000  0.00000  0.000000  0.000000   0.000000         NA          NA        NA
> 3  38.20598  16.34615  12.50000  60.00000 25.42373  3.636364 18.421053  32.335329 -0.3559957      0.1903 0.4812073
> 4 100.00000 100.00000 100.00000 100.00000 86.15385 72.340426 86.666667 100.000000 -0.3093992      0.2195 0.4812073
> 5   0.00000   0.00000   0.00000   0.00000  0.00000  0.000000  0.000000  10.526316  0.4252433      0.3769 0.4812073
> 6  26.98413  12.63158  23.07692  10.00000 22.22222  0.000000  3.174603   4.166667 -0.2946172      0.2334 0.4812073
```

### Extract random correlations from a data.frame

```R
# using the same perc_meth data.frame and ages as defined above
# get 1,000,000 random correlations from the dataset
cors <- sample_n_random_cor(df = perc_meth, 
                            y = ages$age,
                            n = 1000000,
                            cor_method = "spearman")

# simple histogram of correlation values -- drop NAs if present
hist(cors[!is.na(cors)])
```

### Convert edgeR results into a tidy dataframe

```R
library(edgeR)
library(coriell)

# simulate expression data using coriell::simulate_counts()
x <- simulate_counts()

# run edger pipeline
group <-  factor(rep(c("ctl", "trt"), each = 3))
y <- DGEList(counts = x$table, group = group)
y <- calcNormFactors(y)
design <- model.matrix(~group)
y <- estimateDisp(y, design)

# To perform quasi-likelihood F-tests
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = 2)

# -----------------------------------------------------------------------------
# convert to tidy dataframe
res_df <- edger_to_df(qlf, fdr = 1, lfc = 0)

head(res_df)
> # A tibble: 6 x 6
>   feature_id logFC unshrunk.logFC logCPM   PValue        FDR
>   <chr>      <dbl>          <dbl>  <dbl>    <dbl>      <dbl>
> 1 gene.51    -4.43          -4.43   6.37 2.93e-10 0.00000561
> 2 gene.26    -6.63          -6.68   3.45 8.37e-9   0.0000801 
> 3 gene.19    -7.01          -7.05   3.41 1.51e-8   0.0000962 
> 4 gene.100   -4.62          -4.63   4.40 3.24e-8    0.000132  
> 5 gene.74    -6.05          -6.07   7.04 4.06e-8    0.000132  
> 6 gene.53    -4.67          -4.68   3.91 4.14e-8    0.000132 
```

### Summarize results from differential expression test

Return a table of up/down/non-de counts and their percentages.

```R
# for edgeR results default values can be used
res_df <- edger_to_df(qlf)

summarize_dge(res_df, fdr = 0.05)
> # A tibble: 3 x 3
>   dge         n   perc
>   <fct>   <int>  <dbl>
> 1 up        162  0.847
> 2 down      266   1.39 
> 3 non-dge 18709   97.8  

# For DESeq2 results you must specify the column names
summarize_dge(deseq_res_df, 
              fdr_col = padj,
              lfc_col = log2FoldChange,
              fdr = 0.05)
```

### Create volcano plot from differential expression results

```R
# For edgeR results default values can be used
plot_volcano(res_df) +
  ggtitle("Treatment vs Control")

# For DESeq2 results you must specify the column names
plot_volcano(deseq_res_df,
             x = log2FoldChange,
             y = padj)
```

![](man/figures/volcano1.png)

Different significance levels can be used to filter the plotted points. For example,
significance levels can be set by specifying the `fdr` and `lfc` values.

```R
plot_volcano(res_df, fdr = 0.01, lfc = log2(2))  +
  ggtitle("Treatment vs Control")
```

![](man/figures/volcano2.png)

Labels for the counts will be displayed by default. To remove them set `annotate_counts = FALSE`

```R
plot_volcano(res_df, 
             fdr = 0.01, 
             lfc = log2(2), 
             annotate_counts = FALSE)  +
  ggtitle("Treatment vs Control")
```

![](man/figures/volcano3.png)

Positions of the count labels can be adjusted by setting the `xmax_label_offset`, `xmin_label_offset` 
and `ymax_label_offset` values. Setting the values closer to 1 moves the labels away from the origin.
`xmax_label_offset` controls the 'up' gene label whereas `xmin_label_offset` controls the 'down'
genes label. `ymax_label_offset` controls the vertical position of the labels.

```R
plot_volcano(res_df,
             fdr = 0.01,
             lfc = log2(2),
             xmax_label_offset = 0.2, 
             xmin_label_offset = 0.7, 
             ymax_label_offset = 0.9)  +
  ggtitle("Treatment vs Control")
```

![](man/figures/volcano4.png)

Text labels can also be added for the DE genes by setting `label_sig = TRUE`. Caution, if there are many
DE genes this will be overplotted

```R
plot_volcano(res_df,
             fdr = 1e-4,
             lfc = log2(2),
             annotate_counts = FALSE,
             label_sig = TRUE)  +
  ggtitle("Treatment vs Control")
```

![](man/figures/volcano5.png)

### Create md plot from differential expression results

```R
# For edgeR results default values can be used
plot_md(res_df) +
  ggtitle("Treatment vs Control")

# For DESeq2 results you must specify the column names
plot_md(deseq_res_df,
        x = baseMean,
        y = log2FoldChange,
        sig_col = padj)
```

![](man/figures/md1.png)

Different significance levels can be used to filter the plotted points. For example,
significance levels can be set by specifying the `fdr` and `lfc` values.

```R
plot_md(res_df, 
        fdr = 0.01, 
        lfc = log2(2)) +
  ggtitle("Treatment vs Control")
```

![](man/figures/md2.png)

Labels for the counts will be displayed by default. To remove them set `annotate_counts = FALSE`

```R
plot_md(res_df, 
        fdr = 0.01, 
        lfc = log2(2),
        annotate_counts = FALSE) +
  ggtitle("Treatment vs Control")
```

![](man/figures/md3.png)

Positions of the count labels can be adjusted by setting the `xmax_label_offset`, `ymin_label_offset` 
and `ymax_label_offset` values. Setting the values closer to 1 moves the labels away from the origin.
`xmax_label_offset` controls the horizontal position of the labels. `ymin_label_offset` controls the 
'down' gene label. `ymax_label_offset` controls the 'up' genes label.

```R
plot_md(res_df, 
        fdr = 0.01,
        lfc = log2(2),
        xmax_label_offset = 0.6, 
        ymin_label_offset = 0.25, 
        ymax_label_offset = 0.25) +
  ggtitle("Treatment vs Control")
```

![](man/figures/md4.png)

### Heatmap with sensible defaults

We often use the same settings when making calls to `pheatmap`. This function is a wrapper around `pheatmap`
which uses sensible default values for expression data. It changes the default color scale to a diverging
blue to white to red scale, modifies the clustering parameters (row-wise euclidean, col-wise correlation) and
clustering method (complete), angles the column labels, removes border colors and rownames. 

Any of these options can be overridden by simply supplying the arguments to `quickmap` as you would `pheatmap`. 
This also allows for additional arguments to be passed to the `quickmap` function for creating row and column 
annotations. 

```R
# generate some example data and log-scale it
lcpms <- coriell::simulate_counts(n = 1000)$table %>% log1p()

# plot a heatmap of the logCPM values
quickmap(lcpms)
```

![](man/figures/quickmap1.png)

Other `pheatmap` arguments can be passed to the `quickmap` function as well.

```R
# create annotation for columns
col_df <- data.frame(treatment = rep(c("ctl", "trt"), each = 3))
rownames(col_df) <- colnames(lcpms)

# create color scheme for treatment conditions
ann_colors = list(treatment = c("ctl" = "steelblue", "trt" = "firebrick"))

# plot the heatmap, passing additional args to pheatmap
quickmap(lcpms,
         annotation_col = col_df,
         annotation_colors = ann_colors,
         main = "Treatment vs Control")
```

![](man/figures/quickmap2.png)

There are two default color scales included which can be specified by setting the `diverging_palette` argument. 
By default `diverging_palette = TRUE` which sets the color scale the same as the above heatmaps. This is useful for
scaled data. However, if you are plotting unscaled data such as normalized expression values then a continuous color
palette is more appropriate. Setting `diverging_palette = FALSE` will set the color palette to a continuous (`viridis::magma(50)`)
palette.

```R
# NOTE: different lcpms data than above
quickmap(lcpms, diverging_palette = FALSE, scale = "none")
```

![](man/figures/quickmap3.png)

### Z-score a dataframe

Z-score a dataframe by row or column

```R
# create some example data
cpms <- data.frame(a = runif(100, min = 0, max = 100),
                   b = runif(100, min = 0, max = 100),
                   c = runif(100, min = 0, max = 100),
                   d = runif(100, min = 0, max = 100))

> head(cpms)
>           a        b        c         d
> 1 93.586737 41.79316 58.59588 73.082215
> 2 70.009822 25.84383 57.03569 40.512135
> 3 60.053908 14.68176 82.66302 60.842900
> 4 95.711441 86.76281 58.07523 22.571323
> 5  2.833631 25.04612 72.85270  5.417795
> 6 24.596068 85.46398 16.12987 33.810050

# scale the data by rows
cpms_st <- zscore_df(cpms)

> head(cpms_st)
>            a           b          c          d
> 1  1.2201846 -1.13598436 -0.3716029  0.2874026
> 2  1.1247307 -1.16871823  0.4510108 -0.4070232
> 3  0.1922442 -1.39554374  0.9834448  0.2198548
> 4  0.9076295  0.63627273 -0.2336442 -1.3102580
> 5 -0.7309118 -0.04598875  1.4281295 -0.6512290
> 6 -0.4943903  1.45917080 -0.7661138 -0.1986667


# default is to scale by row, scaling by columns can also be performed by
# setting the by = "column"
cpms_st_by_col <- zscore_df(cpms, by = "column")
```

### Convert a list of sets into a binary matrix

This function is useful for comparing if a given gene is present across all
or a certain proportion of conditions.

```R
sets <- list("set1" = letters[1:4],
             "set2" = letters[3:8],
             "set3" = letters[1:5],
             "set4" = letters[4:6])

# convert the sets to a binary matrix
mat <- list_to_matrix(sets)

mat
>   set1 set2 set3 set4
> a    1    0    1    0
> b    1    0    1    0
> c    1    1    1    0
> d    1    1    1    1
> e    0    1    1    1
> f    0    1    0    1
> g    0    1    0    0
> h    0    1    0    0
```

### Get statistics for all pairwise combinations of a list of sets

Compare every set to every other set and return statistics about their intersections.

Statistics are returned in a list object. The returned list contains a named vector of the 
statistic computed. The names of the vectore indicate the pairwise comparison, i.e. "Set A : Set B"
for all combinations of sets. Use `?pairwise_intersection_stats()` for more information about
the statistics computed.

```R
sets <- list("set1" = letters[1:4],
             "set2" = letters[3:8],
             "set3" = letters[1:5],
             "set4" = letters[4:6])

# get the intersection stats -- returns a list of statistics
stats <- pairwise_intersection_stats(sets)

# individual statistics can be accessed with subsetting the list
# to get the intersection sizes of every combination of sets:
stats$intersection_size

> set1 : set2 set1 : set3 set1 : set4 set2 : set3 set2 : set4 set3 : set4 
>           2           4           1           3           3           2
```
