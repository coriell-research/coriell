# coriell

Helper functions for common bioinformatics tasks. If you find yourself reusing old code over and over, let me know and we'll 
create a function and put it in this package. 

## Installation

First, install `devtools` with:

`install.packages("devtools")`

Then install `coriell` with:

`devtools::install_github('jcalendo/coriell')`

## Currently Available Functions

`plot_volcano` and `plot_MD` produce volcano and MD plots using `edgeR` results objects as input.

### Example

```{r}
library(coriell)
library(edgeR)


x <- read.delim("TableOfCounts.txt",row.names="Symbol")
group <- factor(c(1,1,2,2))
y <- DGEList(counts=x,group=group)
y <- calcNormFactors(y)
design <- model.matrix(~group)
y <- estimateDisp(y,design)

# To perform quasi-likelihood F-tests:
fit <- glmQLFit(y,design)
qlf <- glmQLFTest(fit,coef=2)

# call the coriell plotting functions
plot_volcano(qlf, graph_title = "Group1 vs. Group2 DGE - volcano plot")
plot_MD(qlf, graph_title = "Group1 vs. Group2 DGE- MD plot")
```

`parallel_permutation_test` performs a correlation based permutation test on a data.frame using multiple CPU cores. This version of the permutation test is at least 10X faster and more memory efficient than older code.

### Example

```{r}
library(coriell)

# create some fake methylation data
# six samples, three control and 3 treatment, 1000 rows
perc_meth <- data.frame(ctl1 = runif(n = 1000, min = 0, max = 100),
                        ctl2 = runif(n = 1000, min = 0, max = 100),
                        ctl3 = runif(n = 1000, min = 0, max = 100),
                        trt1 = runif(n = 1000, min = 0, max = 100),
                        trt2 = runif(n = 1000, min = 0, max = 100),
                        trt3 = runif(n = 1000, min = 0, max = 100))

# create a vector of sample ages matching the column order                       
ages <- c(20, 34, 65, 25, 42, 60)

# perform the permutation test
perm_res <- parallel_permutation_test(perc_meth, y = ages, n_cores = 12, n_perm = 1000)
```
