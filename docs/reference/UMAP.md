# Perform UMAP

This function provides a wrapper around `umap::umap()` that exposes the
umap defaults as function arguments.

## Usage

``` r
UMAP(x, ...)

# S3 method for class 'pca'
UMAP(
  x,
  n_neighbors = 15,
  n_components = 2,
  metric = "euclidean",
  n_epochs = 200,
  input = "data",
  init = "spectral",
  min_dist = 0.1,
  set_op_mix_ratio = 1,
  local_connectivity = 1,
  bandwidth = 1,
  alpha = 1,
  gamma = 1,
  negative_sample_rate = 5,
  a = NA,
  b = NA,
  spread = 1,
  random_state = NA,
  transform_state = NA,
  knn = NA,
  knn_repeats = 1,
  verbose = FALSE,
  umap_learn_args = NA
)

# S3 method for class 'prcomp'
UMAP(
  x,
  metadata = NULL,
  n_neighbors = 15,
  n_components = 2,
  metric = "euclidean",
  n_epochs = 200,
  input = "data",
  init = "spectral",
  min_dist = 0.1,
  set_op_mix_ratio = 1,
  local_connectivity = 1,
  bandwidth = 1,
  alpha = 1,
  gamma = 1,
  negative_sample_rate = 5,
  a = NA,
  b = NA,
  spread = 1,
  random_state = NA,
  transform_state = NA,
  knn = NA,
  knn_repeats = 1,
  verbose = FALSE,
  umap_learn_args = NA
)

# S3 method for class 'matrix'
UMAP(
  x,
  metadata = NULL,
  n_neighbors = 15,
  n_components = 2,
  metric = "euclidean",
  n_epochs = 200,
  input = "data",
  init = "spectral",
  min_dist = 0.1,
  set_op_mix_ratio = 1,
  local_connectivity = 1,
  bandwidth = 1,
  alpha = 1,
  gamma = 1,
  negative_sample_rate = 5,
  a = NA,
  b = NA,
  spread = 1,
  random_state = NA,
  transform_state = NA,
  knn = NA,
  knn_repeats = 1,
  verbose = FALSE,
  umap_learn_args = NA
)

# S3 method for class 'data.frame'
UMAP(
  x,
  metadata = NULL,
  n_neighbors = 15,
  n_components = 2,
  metric = "euclidean",
  n_epochs = 200,
  input = "data",
  init = "spectral",
  min_dist = 0.1,
  set_op_mix_ratio = 1,
  local_connectivity = 1,
  bandwidth = 1,
  alpha = 1,
  gamma = 1,
  negative_sample_rate = 5,
  a = NA,
  b = NA,
  spread = 1,
  random_state = NA,
  transform_state = NA,
  knn = NA,
  knn_repeats = 1,
  verbose = FALSE,
  umap_learn_args = NA
)

# S3 method for class 'dist'
UMAP(
  x,
  metadata = NULL,
  n_neighbors = 15,
  n_components = 2,
  metric = "euclidean",
  n_epochs = 200,
  input = "dist",
  init = "spectral",
  min_dist = 0.1,
  set_op_mix_ratio = 1,
  local_connectivity = 1,
  bandwidth = 1,
  alpha = 1,
  gamma = 1,
  negative_sample_rate = 5,
  a = NA,
  b = NA,
  spread = 1,
  random_state = NA,
  transform_state = NA,
  knn = NA,
  knn_repeats = 1,
  verbose = FALSE,
  umap_learn_args = NA
)
```

## Arguments

- x:

  PCA object, prcomp object, or numeric matrix/data.frame that can be
  converted to a numeric matrix

- n_neighbors:

  Number of nearest neighbors. Default 15

- n_components:

  Dimension of target (output) space. Default 2

- metric:

  character or function; determines how distances between data points
  are computed. When using a string, available metrics are: euclidean,
  manhattan. Other available generalized metrics are: cosine, pearson,
  pearson2. Note the triangle inequality may not be satisfied by some
  generalized metrics, hence knn search may not be optimal. When using
  metric.function as a function, the signature must be function(matrix,
  origin, target) and should compute a distance between the origin
  column and the target columns. Default "euclidean"

- n_epochs:

  Number of iterations performed during layout optimization. Default 200

- input:

  character, use either "data" or "dist"; determines whether the primary
  input argument to umap() is treated as a data matrix or as a distance
  matrix. Default "data"

- init:

  character or matrix. The default string "spectral" computes an initial
  embedding using eigenvectors of the connectivity graph matrix. An
  alternative is the string "random", which creates an initial layout
  based on random coordinates. This setting.can also be set to a matrix,
  in which case layout optimization begins from the provided
  coordinates. Default "spectral"

- min_dist:

  numeric; determines how close points appear in the final layout.
  Default 0.1

- set_op_mix_ratio:

  numeric in range \[0,1\]; determines who the knn-graph is used to
  create a fuzzy simplicial graph. Default 1

- local_connectivity:

  numeric; used during construction of fuzzy simplicial set. Default 1

- bandwidth:

  numeric; used during construction of fuzzy simplicial set. Default 1

- alpha:

  numeric; initial value of "learning rate" of layout optimization.
  Default 1

- gamma:

  numeric; determines, together with alpha, the learning rate of layout
  optimization. Default 1

- negative_sample_rate:

  integer; determines how many non-neighbor points are used per point
  and per iteration during layout optimization. Default 5

- a:

  numeric; contributes to gradient calculations during layout
  optimization. When left at NA, a suitable value will be estimated
  automatically. Default NA

- b:

  numeric; contributes to gradient calculations during layout
  optimization. When left at NA, a suitable value will be estimated
  automatically. Default NA

- spread:

  numeric; used during automatic estimation of a/b parameters. Default 1

- random_state:

  integer; seed for random number generation used during umap(). Default
  NA

- transform_state:

  nteger; seed for random number generation used during predict().
  Default NA

- knn:

  object of class umap.knn; precomputed nearest neighbors. Default NA

- knn_repeats:

  number of times to restart knn search. Default 1

- verbose:

  logical or integer; determines whether to show progress messages.
  Default FALSE

- umap_learn_args:

  vector of arguments to python package umap-learn. Default NA

- metadata:

  Optional data.frame with sample-level metadata. Used if a prcomp
  object or data.frame/matrix is supplied. Default NULL

## Value

data.frame with the UMAP embeddings. If metadata was supplied then
metadata columns are added to the results.

## Examples

``` r
# Create metadata for plotting
metadata <- data.frame(row.names = colnames(GSE161650_lc))
metadata$Group <- rep(c("DMSO", "THZ1"), each = 3)

# PCA with PCAtools
p <- PCAtools::pca(GSE161650_lc, metadata, center = TRUE, scale = TRUE)

# PCA with prcomp
pr <- prcomp(t(GSE161650_lc), center = TRUE, scale. = FALSE)

# Pre-calculated distance matrix
d <- dist(t(GSE161650_lc))

# Perform UMAP on each data type
udata <- UMAP(p, n_neighbors = 2)
#> Error in loadNamespace(x): there is no package called ‘umap’
udata2 <- UMAP(pr, metadata, n_neighbors = 2)
#> Error in loadNamespace(x): there is no package called ‘umap’
udata3 <- UMAP(d, metadata, n_neighbors = 2)
#> Error in loadNamespace(x): there is no package called ‘umap’

# Also on raw data
udata4 <- UMAP(t(GSE161650_lc), metadata, n_neighbors = 2)
#> Error in loadNamespace(x): there is no package called ‘umap’
```
