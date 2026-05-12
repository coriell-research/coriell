# coriell

This package is a collection of helper functions for common bioinformatics 
tasks. It includes plotting functions for expression data, wrappers around 
typical bioinformatics tools to ease use, some data transformation functions, 
dimensionality reduction helpers, and more. 

Check out the [package website](https://coriell-research.github.io/coriell/) 
for examples and documentation. 

Also check out the 
[Coriell bioinformatics notes](https://coriell-research.github.io/coriell-bioinformatics-notes/) 
bookdown book.

## Installation

Since this package will be constantly changing be sure to install the latest 
version from Github using:

```R
if (!require("pak")) {
  install.packages("pak")
}
pak::pkg_install("coriell-research/coriell")
```
