# coriell

Helper functions for common bioinformatics tasks. If you find yourself reusing old code over and over, let me know and we'll 
create a function and put it in this package. 

**Please test and let me know when things are broken**

## Installation

Since this package will be constantly changing be sure to install the latest version from github using:

```{r}
if (!require("devtools")) {
  install.packages("devtools")
}

devtools::install_github("coriell-research/coriell")
```

## Examples

### Perform correlation permutation test using multiple cores

Note: if the number of possible permutations of your data is less than the desired number
of permutations given by the parameter `n_perm` then an exact test on all of the real permutations
will be performed instead of random sampling. For example, if you only have 6 samples (720 possible permutations)
but set `n_perm = 1000` only 720 permutations (i.e. the exact test) will be tested.
The `permutation_correlation_test` function will display a message if this occurs.

```{r}
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

```{r}
# using the same perc_meth data.frame and ages as defined above
# get 1,000,000 random correlations from the dataset
cors <- sample_n_random_cor(df = perc_meth, 
                            y = ages$age,
                            n = 1000000,
                            cor_method = "spearman")

# simple histogram of correlation values -- drop NAs if present
hist(cors[!is.na(cors)])
```

### Plot and summarize results from edgeR analysis

```{r}
library(edgeR)
library(coriell)

# create some fake expression data
x <- data.frame(ctl1 = rnbinom(1000, size = 0.4, prob = 1e-5),
                ctl2 = rnbinom(1000, size = 0.4, prob = 1e-5),
                trt1 = rnbinom(1000, size = 0.4, prob = 1e-5),
                trt2 = rnbinom(1000, size = 0.4, prob = 1e-5),
                row.names = paste0('gene', 1:1000))

# run edger pipeline
group <- factor(c(1,1,2,2))
y <- DGEList(counts = x,group = group)
y <- calcNormFactors(y)
design <- model.matrix(~group)
y <- estimateDisp(y, design)

# To perform quasi-likelihood F-tests:
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = 2)

# coriell functions ------------------------------------------------------------
# convert the results object to a tibble
res_df <- edger_to_df(qlf)

# results can be filtered with optional arguments
edger_to_df(qlf, fdr = 0.01, lfc = 2)

# Create a volcano plot (use defaults)
plot_volcano(res_df)

# Additional plot parameters can be specified, for example adding fdr/lfc cutoffs
plot_volcano(res_df, fdr = 0.01, lfc = 2)

# Create an md plot (use defaults)
plot_md(res_df)

# Additional plot parameters can be specified, for example adding fdr/lfc cutoffs
plot_md(res_df, fdr = 0.01, lfc = 2)

# Summarize the counts into table of up/down/non-de
summarize_dge(res_df)

# Summary table can be filtered as well
summarize_dge(res_df, fdr = 0.05, lfc = 0)
```
