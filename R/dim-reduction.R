#' Calculate associations between variables and principal components
#'
#' Calculate associations between metadata variables and PCA rotations. This
#' function is inspired by \code{methylKit::assocComp} but is designed to work
#' with arbitrary input data.
#' @details This function returns the p-values from testing the association
#' between a given metadata column and all PC rotations. For numeric values the
#' p-value returned is computed using the \code{cor.test()} function. For
#' factor variables and variables that can be converted to factor variables the
#' function will return p-values from the \code{wilcox.test()} function or the
#' \code{kruskal.test()} function (when the metadata variable has > 2 levels).
#' @param x data.frame or matrix of values passed to \code{prcomp()}
#' @param metadata data.frame of metadata variables to test associations. rownames(metadata) must be identical to colnames(x).
#' @param N First N PCs to test associations between. Default 10.
#' @param ... Additional arguments passed to the \code{prcomp()} function.
#' @export
#' @return data.frame with rows for each metadata variable, columns for each PC, and p-values from the given test in cells.
#' @examples
#' # Specify metadata
#' metadata <- data.frame(
#'   age = c(30, 80, 34, 30, 80, 40),
#'   treatment = factor(c(rep("Treatment", 3), rep("Control", 3))),
#'   class = factor(c(rep("A", 2), rep("B", 2), rep("C", 2))),
#'   row.names = c(paste0("trt", 1:3), paste0("ctrl", 1:3))
#' )
#'
#' # Create values to perform PCA on
#' df <- data.frame(replicate(6, runif(1000, 0, 100)))
#' colnames(df) <- c(paste0("trt", 1:3), paste0("ctrl", 1:3))
#'
#' # Test for associations
#' res <- associate_components(df, metadata)
#'
#' # Show results
#' head(res)
associate_components <- function(x, metadata, N = 10, ...) {
  stopifnot("colnames of x and rownames of metadata do not match." = all(colnames(x) == rownames(metadata)))
  stopifnot("No metadata columns present" = ncol(metadata) >= 1)
  stopifnot("NA values detected in metadata" = all(sapply(metadata, anyNA) == FALSE))

  # Perform PCA
  res <- prcomp(x, ...)

  # Test only the specified number of PCs
  keep_cols <- if (N > ncol(res$rotation)) ncol(res$rotation) else N
  rotation <- res$rotation[, 1:keep_cols]

  # Test the association of each metadata column with each PC rotation
  results <- matrix(data = NA, nrow = ncol(metadata), ncol = ncol(rotation))
  for (i in seq_along(metadata)) {
    v <- metadata[, i]
    for (j in 1:ncol(rotation)) {
      if (is.factor(v) | is.character(v) | is.logical(v)) {
        v <- as.factor(v)
        d <- split(x = rotation[, j], f = v)
        if (nlevels(v) == 2) {
          results[i, j] <- wilcox.test(x = d[[1]], y = d[[2]])$p.value
        } else {
          results[i, j] <- kruskal.test(d)$p.value
        }
      } else {
        results[i, j] <- cor.test(v, rotation[, j])$p.value
      }
    }
  }
  rownames(results) <- colnames(metadata)
  colnames(results) <- colnames(rotation)
  as.data.frame(results)
}

#' Remove principal components from data
#'
#' Reconstruct data after removing principal components. This function will
#' reconstruct a truncated version of the original data matrix after removing
#' the specified principal components. The source code for this function was
#' stolen from \href{https://stats.stackexchange.com/questions/57467/how-to-perform-dimensionality-reduction-with-pca-in-r/57478#57478}{stack exchange}
#' @param x data.frame or matrix of original data.
#' @param components numeric vector of components to remove from original data. Default (1)
#' @param ... Additional arguments passed to \code{prcomp()}
#' @export
#' @return data.frame of values in the original units after removing components
#' @examples
#' # Remove first two components from dataset
#' trunc <- remove_components(USArrests, components = 1:2, scale. = TRUE, center = TRUE)
#'
#' # original data
#' head(USArrests)
#'
#' # reconstructed data
#' head(trunc)
remove_components <- function(x, components = 1, ...) {
  stopifnot("Cannot remove more components than columns in data" = max(components) < ncol(x))

  res <- prcomp(x, ...)
  keep <- setdiff(1:ncol(res$x), components)
  trunc <- res$x[, keep] %*% t(res$rotation[, keep])

  if (res$scale[[1]] != FALSE & length(res$scale) > 1) {
    trunc <- scale(trunc, center = FALSE, scale = 1 / res$scale)
  }
  if (res$center[[1]] != FALSE & length(res$center) > 1) {
    trunc <- scale(trunc, center = -1 * res$center, scale = FALSE)
  }
  data.frame(trunc)
}


#' Perform UMAP
#'
#' This function provides a wrapper around \code{umap::umap()} that exposes
#' the umap defaults as function arguments.
#' @param x PCA object, prcomp object, or numeric matrix/data.frame that can be
#' converted to a numeric matrix
#' @param metadata Optional data.frame with sample-level metadata. Used if a 
#' prcomp object or data.frame/matrix is supplied. Default NULL
#' @param n_neighbors Number of nearest neighbors. Default 15
#' @param n_components Dimension of target (output) space. Default 2
#' @param metric character or function; determines how distances between data
#' points are computed. When using a string, available metrics are: euclidean,
#' manhattan. Other available generalized metrics are: cosine, pearson,
#' pearson2. Note the triangle inequality may not be satisfied by some
#' generalized metrics, hence knn search may not be optimal. When using
#' metric.function as a function, the signature must be function(matrix,
#' origin, target) and should compute a distance between the origin column and
#' the target columns. Default "euclidean"
#' @param n_epochs Number of iterations performed during layout optimization.
#' Default 200
#' @param input character, use either "data" or "dist"; determines whether the
#' primary input argument to umap() is treated as a data matrix or as a
#' distance matrix. Default "data"
#' @param init character or matrix. The default string "spectral" computes an
#' initial embedding using eigenvectors of the connectivity graph matrix. An
#' alternative is the string "random", which creates an initial layout based on
#' random coordinates. This setting.can also be set to a matrix, in which case
#' layout optimization begins from the provided coordinates. Default "spectral"
#' @param min_dist numeric; determines how close points appear in the final
#' layout. Default 0.1
#' @param set_op_mix_ratio numeric in range [0,1]; determines who the knn-graph
#' is used to create a fuzzy simplicial graph. Default 1
#' @param local_connectivity numeric; used during construction of fuzzy
#' simplicial set. Default 1
#' @param bandwidth numeric; used during construction of fuzzy simplicial set.
#' Default 1
#' @param alpha numeric; initial value of "learning rate" of layout optimization.
#' Default 1
#' @param gamma numeric; determines, together with alpha, the learning rate of
#' layout optimization. Default 1
#' @param negative_sample_rate  integer; determines how many non-neighbor
#' points are used per point and per iteration during layout optimization.
#' Default 5
#' @param a numeric; contributes to gradient calculations during layout
#' optimization. When left at NA, a suitable value will be estimated
#' automatically. Default NA
#' @param b numeric; contributes to gradient calculations during layout
#' optimization. When left at NA, a suitable value will be estimated
#' automatically. Default NA
#' @param spread numeric; used during automatic estimation of a/b parameters.
#' Default 1
#' @param random_state integer; seed for random number generation used during
#' umap(). Default NA
#' @param transform_state nteger; seed for random number generation used
#' during predict(). Default NA
#' @param knn object of class umap.knn; precomputed nearest neighbors. Default
#' NA
#' @param knn_repeats number of times to restart knn search. Default 1
#' @param verbose logical or integer; determines whether to show progress
#' messages. Default FALSE
#' @param umap_learn_args vector of arguments to python package umap-learn.
#' Default NA
#' @export
#' @return data.frame with the UMAP embeddings. If metadata was supplied then 
#' metadata columns are added to the results.
#' @examples
#' 
#' # Create metadata for plotting
#' metadata <- data.frame(row.names = colnames(GSE161650_lc))
#' metadata$Group <- rep(c("DMSO", "THZ1"), each = 3)
#' 
#' # PCA with PCAtools
#' p <- PCAtools::pca(GSE161650_lc, metadata, center = TRUE, scale = TRUE)
#' 
#' # PCA with prcomp
#' pr <- prcomp(t(GSE161650_lc), center = TRUE, scale. = FALSE)
#' 
#' # Pre-calculated distance matrix
#' d <- dist(t(GSE161650_lc))
#' 
#' # Perform UMAP on each data type
#' udata <- UMAP(p, n_neighbors = 2)
#' udata2 <- UMAP(pr, metadata, n_neighbors = 2)
#' udata3 <- UMAP(d, metadata, n_neighbors = 2)
#' 
#' # Also on raw data
#' udata4 <- UMAP(t(GSE161650_lc), metadata, n_neighbors = 2)
#' 
UMAP <- function(x, ...) UseMethod("UMAP")

UMAP.default <- function(x) {
  stop("Object of type ", class(x), " is not supported by this function")
}

#' @rdname UMAP
#' @export
#'
UMAP.pca <- function(x, n_neighbors = 15, n_components = 2, metric = "euclidean",
                     n_epochs = 200, input = "data", init = "spectral",
                     min_dist = 0.1, set_op_mix_ratio = 1, local_connectivity = 1,
                     bandwidth = 1, alpha = 1, gamma = 1, negative_sample_rate = 5,
                     a = NA, b = NA, spread = 1, random_state = NA,
                     transform_state = NA, knn = NA, knn_repeats = 1,
                     verbose = FALSE, umap_learn_args = NA) {
  # Experiment with the code below to more elegantly replace values
  # params <- umap::umap.defaults
  # user_args <- c(as.list(environment()))
  # new_args <- modifyList(params, user_args)

  params <- umap::umap.defaults
  params$n_neighbors <- n_neighbors
  params$n_components <- n_components
  params$metric <- metric
  params$n_epochs <- n_epochs
  params$input <- input
  params$init <- init
  params$min_dist <- min_dist
  params$set_op_mix_ratio <- set_op_mix_ratio
  params$local_connectivity <- local_connectivity
  params$bandwidth <- bandwidth
  params$alpha <- alpha
  params$gamma <- gamma
  params$negative_sample_rate <- negative_sample_rate
  params$a <- a
  params$b <- b
  params$spread <- spread
  params$random_state <- random_state
  params$transform_state <- transform_state
  params$knn <- knn
  params$knn_repeats <- knn_repeats
  params$verbose <- verbose
  params$umap_learn_args <- umap_learn_args

  # Perform UMAP and add on metadata information
  u <- umap::umap(d = x$rotated, params)
  data <- u$layout
  colnames(data) <- paste0("UMAP", 1:ncol(data))
  result <- cbind(data, x$metadata)

  return(result)
}

#' @rdname UMAP
#' @export
#'
UMAP.prcomp <- function(x, metadata = NULL, n_neighbors = 15, n_components = 2, 
                        metric = "euclidean", n_epochs = 200, input = "data", 
                        init = "spectral", min_dist = 0.1, set_op_mix_ratio = 1, 
                        local_connectivity = 1, bandwidth = 1, alpha = 1, 
                        gamma = 1, negative_sample_rate = 5, a = NA, b = NA, 
                        spread = 1, random_state = NA, transform_state = NA, 
                        knn = NA, knn_repeats = 1, verbose = FALSE, 
                        umap_learn_args = NA) {

  if (!is.null(metadata)) {
    if (!all(rownames(x$x) == rownames(metadata))) {
      stop("rownames of prcomp object != rownames of metadata")
    }
  }
  
  params <- umap::umap.defaults
  params$n_neighbors <- n_neighbors
  params$n_components <- n_components
  params$metric <- metric
  params$n_epochs <- n_epochs
  params$input <- input
  params$init <- init
  params$min_dist <- min_dist
  params$set_op_mix_ratio <- set_op_mix_ratio
  params$local_connectivity <- local_connectivity
  params$bandwidth <- bandwidth
  params$alpha <- alpha
  params$gamma <- gamma
  params$negative_sample_rate <- negative_sample_rate
  params$a <- a
  params$b <- b
  params$spread <- spread
  params$random_state <- random_state
  params$transform_state <- transform_state
  params$knn <- knn
  params$knn_repeats <- knn_repeats
  params$verbose <- verbose
  params$umap_learn_args <- umap_learn_args
  
  # Perform UMAP and add on metadata information
  u <- umap::umap(d = x$x, params)
  data <- u$layout
  colnames(data) <- paste0("UMAP", 1:ncol(data))
  if (is.null(metadata)) return(data)
  result <- cbind(data, metadata)
  
  return(result)
}

#' @rdname UMAP
#' @export
#'
UMAP.matrix <- function(x, metadata = NULL, n_neighbors = 15, n_components = 2, 
                        metric = "euclidean", n_epochs = 200, input = "data", 
                        init = "spectral", min_dist = 0.1, set_op_mix_ratio = 1, 
                        local_connectivity = 1, bandwidth = 1, alpha = 1, 
                        gamma = 1, negative_sample_rate = 5, a = NA, b = NA, 
                        spread = 1, random_state = NA, transform_state = NA, 
                        knn = NA, knn_repeats = 1, verbose = FALSE, 
                        umap_learn_args = NA) {
  
  if (!is.null(metadata)) {
    if (!all(rownames(x) == rownames(metadata))) {
      stop("rownames of matrix != rownames of metadata.")
    }
  }
  
  params <- umap::umap.defaults
  params$n_neighbors <- n_neighbors
  params$n_components <- n_components
  params$metric <- metric
  params$n_epochs <- n_epochs
  params$input <- input
  params$init <- init
  params$min_dist <- min_dist
  params$set_op_mix_ratio <- set_op_mix_ratio
  params$local_connectivity <- local_connectivity
  params$bandwidth <- bandwidth
  params$alpha <- alpha
  params$gamma <- gamma
  params$negative_sample_rate <- negative_sample_rate
  params$a <- a
  params$b <- b
  params$spread <- spread
  params$random_state <- random_state
  params$transform_state <- transform_state
  params$knn <- knn
  params$knn_repeats <- knn_repeats
  params$verbose <- verbose
  params$umap_learn_args <- umap_learn_args
  
  # Perform UMAP and add on metadata information
  u <- umap::umap(d = x, params)
  data <- u$layout
  colnames(data) <- paste0("UMAP", 1:ncol(data))
  if (is.null(metadata)) return(data)
  result <- cbind(data, metadata)
  
  return(result)
}

#' @rdname UMAP
#' @export
#'
UMAP.data.frame <- function(x, metadata = NULL, n_neighbors = 15, n_components = 2, 
                        metric = "euclidean", n_epochs = 200, input = "data", 
                        init = "spectral", min_dist = 0.1, set_op_mix_ratio = 1, 
                        local_connectivity = 1, bandwidth = 1, alpha = 1, 
                        gamma = 1, negative_sample_rate = 5, a = NA, b = NA, 
                        spread = 1, random_state = NA, transform_state = NA, 
                        knn = NA, knn_repeats = 1, verbose = FALSE, 
                        umap_learn_args = NA) {
  m <- data.matrix(x)
  if (!is.null(metadata)) {
    if (!all(rownames(m) == rownames(metadata))) {
      stop("rownames of matrix != rownames of metadata. ")
    }
  }
  
  params <- umap::umap.defaults
  params$n_neighbors <- n_neighbors
  params$n_components <- n_components
  params$metric <- metric
  params$n_epochs <- n_epochs
  params$input <- input
  params$init <- init
  params$min_dist <- min_dist
  params$set_op_mix_ratio <- set_op_mix_ratio
  params$local_connectivity <- local_connectivity
  params$bandwidth <- bandwidth
  params$alpha <- alpha
  params$gamma <- gamma
  params$negative_sample_rate <- negative_sample_rate
  params$a <- a
  params$b <- b
  params$spread <- spread
  params$random_state <- random_state
  params$transform_state <- transform_state
  params$knn <- knn
  params$knn_repeats <- knn_repeats
  params$verbose <- verbose
  params$umap_learn_args <- umap_learn_args
  
  # Perform UMAP and add on metadata information
  u <- umap::umap(d = m, params)
  data <- u$layout
  colnames(data) <- paste0("UMAP", 1:ncol(data))
  if (is.null(metadata)) return(data)
  result <- cbind(data, metadata)
  
  return(result)
}

#' @rdname UMAP
#' @export
#'
UMAP.dist <- function(x, metadata = NULL, n_neighbors = 15, n_components = 2, 
                      metric = "euclidean", n_epochs = 200, input = "dist", 
                      init = "spectral", min_dist = 0.1, set_op_mix_ratio = 1, 
                      local_connectivity = 1, bandwidth = 1, alpha = 1, 
                      gamma = 1, negative_sample_rate = 5, a = NA, b = NA, 
                      spread = 1, random_state = NA, transform_state = NA, 
                      knn = NA, knn_repeats = 1, verbose = FALSE, 
                      umap_learn_args = NA) {
  
  if (!is.null(metadata)) {
    if (!all(attr(x, "Labels") == rownames(metadata))) {
      stop("Labels of distance matrix != rownames of metadata.")
    }
  }
  
  params <- umap::umap.defaults
  params$n_neighbors <- n_neighbors
  params$n_components <- n_components
  params$metric <- metric
  params$n_epochs <- n_epochs
  params$input <- input
  params$init <- init
  params$min_dist <- min_dist
  params$set_op_mix_ratio <- set_op_mix_ratio
  params$local_connectivity <- local_connectivity
  params$bandwidth <- bandwidth
  params$alpha <- alpha
  params$gamma <- gamma
  params$negative_sample_rate <- negative_sample_rate
  params$a <- a
  params$b <- b
  params$spread <- spread
  params$random_state <- random_state
  params$transform_state <- transform_state
  params$knn <- knn
  params$knn_repeats <- knn_repeats
  params$verbose <- verbose
  params$umap_learn_args <- umap_learn_args
  
  # Perform UMAP and add on metadata information
  m <- as.matrix(x)
  u <- umap::umap(d = m, params)
  data <- u$layout
  colnames(data) <- paste0("UMAP", 1:ncol(data))
  if (is.null(metadata)) return(data)
  result <- cbind(data, metadata)
  
  return(result)
}
