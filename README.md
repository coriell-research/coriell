# coriell

Helper functions for common bioinformatics tasks. If you find yourself reusing old code over and over, let me know and we'll 
create a function and put it in this package. 

## Installation

First, install `devtools` with:

`install.packages("devtools")`

Then install `coriell` with:

`devtools::install_github('jcalendo/coriell')`

## Usage

Currently available functions:

- `plot_volcano` produces volcano plots using `edgeR` results objects as input. 
- `plot_MD` produces MD plots using `edgeR` results objects as input. 
- `parallel_permutation_test` performs a correlation based permutation test on a data.frame using multiple CPU cores. This version of the permutation test is at least 10X faster and more memory efficient than older code.
