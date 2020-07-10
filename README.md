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

```{r}
library(coriell)
library(methylKit)

# define age dataframe -- ages matching column order
ages = data.frame(age = c(30, 80, 34, 30, 80, 40))

# simulate a methylation dataset
sim_meth <- dataSim(
  replicates = 6,
  sites = 1000,
  treatment = c(rep(1, 3), rep(0, 3)),
  covariates = ages,
  sample.ids = c(paste0("test", 1:3), paste0("ctrl", 1:3))
)

# extract the methylation as percentages and coerce to data.frame
perc_meth <- as.data.frame(percMethylation(sim_meth))

head(perc_meth)
>       test1     test2     test3     ctrl1      ctrl2     ctrl3
> 1  12.28070  44.11765  39.13043  28.57143  80.000000  32.25806
> 2  26.19048  58.10811  64.10256  53.08642  50.574713  11.11111
> 3  12.16216  10.52632  29.03226   6.47482   6.451613  34.17085
> 4   0.00000   0.00000   0.00000   0.00000   0.000000   0.00000
> 5  80.00000  58.06452  80.00000  74.78261  90.000000  50.00000
> 6 100.00000 100.00000 100.00000 100.00000 100.000000 100.00000

# permutation testing -----------------------------------------------------------
# perform permutation testing using 4 cores and 1000 permutations
res <- permutation_correlation_test(perc_meth, 
                                    y = ages$age, 
                                    n_cores = 4, 
                                    n_perm = 1000,
                                    cor_method = "spearman", 
                                    p_adjust_method = "fdr")

head(res)
>       test1     test2     test3     ctrl1      ctrl2     ctrl3         cor empirical_p       fdr
> 1  12.28070  44.11765  39.13043  28.57143  80.000000  32.25806  0.91215932       0.013 0.4116061
> 2  26.19048  58.10811  64.10256  53.08642  50.574713  11.11111  0.08827348       0.428 0.5180415
> 3  12.16216  10.52632  29.03226   6.47482   6.451613  34.17085 -0.20597146       0.345 0.5039347
> 4   0.00000   0.00000   0.00000   0.00000   0.000000   0.00000          NA          NA        NA
> 5  80.00000  58.06452  80.00000  74.78261  90.000000  50.00000 -0.04478111       0.429 0.5180415
> 6 100.00000 100.00000 100.00000 100.00000 100.000000 100.00000          NA          NA        NA
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
