#' Min-Max normalize a value within a given range
#'
#' @param x numeric. Value to be transformed (normalized).
#' @param min numeric. Minimum value of range. Default NA, use the min of the supplied vector
#' @param max numeric. Maximum value of range. Default NA, use the max of the supplied vector
#' @return Normalized value (0-1) of x in the range given by min, max.
#' @export
#' @examples
#' minmax(c(0, 1, 10, 100))
#' 
#' # Use a scale outside of the range of x
#' minmax(c(0, 1, 10, 100), min = 0, max = 1000)
minmax <- function(x, min = NA, max = NA) {
  (x - min(c(x, min), na.rm = TRUE)) / 
  (max(c(x, max), na.rm = TRUE) - min(c(x, min), na.rm = TRUE))
}

#' Linear interpolation of a value
#'
#' @param x. numeric. Normalized input value to be transformed. Value must be between 0-1.
#' @param min numeric. Minimum of the desired range.
#' @param max numeric. Maximum of the desired range.
#' @return Value of x within the desired range given by min, max.
#' @export
#' @examples
#' lerp(c(0.1, 0.25, 0.5, 0.75), min = 0, max = 100)
lerp <- function(x, min, max) {
  stopifnot("x must be between 0-1 (inclusive)" = x >= 0 & x <= 1)
  return((max - min) * x + min)
}

#' Map a value in one range to a value in another
#'
#' @param x numeric. Value to be transformed
#' @param old_min numeric. Minimum value of source range.
#' @param old_max numeric. Maximum value of source range.
#' @param new_min numeric. Minimum value of desired range.
#' @param new_max numeric. Maximum value of desired range.
#' @return Value of x mapped to the desired range.
#' @export
#' @examples
#' # Map values from the range 0-10 to the range 0-1
#' map_value(c(0, 1, 2, 3), 0, 10, 0, 1)
map_value <- function(x, old_min, old_max, new_min, new_max) {
  lerp(minmax(x, old_min, old_max), new_min, new_max)
}

.clamp <- function(x, min, max) min(c(max(c(x, min)), max))

#' Limit values to a given range
#'
#' This function is useful for input validation. If the value of x is in the given
#' range then the function returns x. If the value is outside of the range then
#' the function returns either the max or the min value of the range.
#'
#' @param x numeric. Value to validate
#' @param min numeric. Minimum value of range.
#' @param max numeric. Maximum value of range.
#' @return The clamped value.
#' @export
#' @examples
#' clamp(c(-1, 121, 10, 15), min = 0, max = 100)
clamp <- Vectorize(.clamp)

#' Geometric mean of a vector
#'
#' Calculate the geometric mean of a vector
#' Stolen from \url{https://stackoverflow.com/a/25555105}.
#'
#' @param x numeric vector of non-negative values
#' @param na.rm logical. Remove NAs from calculation
#' @param zero_propagate logical. Should zeros be included in the calculation. Default FALSE.
#' @return geometric mean of the vector
#' @export
#' @examples
#' x <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 1500)
#' geometric_mean(x)
geometric_mean <- function(x, zero_propagate = FALSE, na.rm = TRUE) {
  if (any(x < 0, na.rm = TRUE)) {
    return(NaN)
  }
  if (zero_propagate) {
    if (any(x == 0, na.rm = TRUE)) {
      return(0)
    }
    exp(mean(log(x), na.rm = na.rm))
  } else {
    exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
  }
}

#' Centered Log-ratio transformation
#'
#' Calculate the centered log-ratio transformation of a vector
#'
#' @param x numeric vector
#' @param base integer base of the log function. Default 2.
#' @export
#' @examples
#' x <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 1500)
#' clr(x)
clr <- function(x, base = 2) {
  x <- log((x / geometric_mean(x)), base)
  x[!is.finite(x) | is.na(x)] <- 0.0

  return(x)
}

#' Give means of rows of matrix based on column grouping variable
#'
#' Inspired by \code{edgeR::sumTechReps} and \code{base::rowsum()}, this function
#' takes the average of the values in each group given by the group argument
#' for each row of the data matrix.
#' @param x variable (gene) by sample numeric matrix
#' @param group Factor specifying the grouping level to by averaged
#' @return variable x nLevels(group) matrix
#' @export
#' @examples
#' # Specify the Group levels
#' Group <- gl(n = 2, k = 3, labels = c("DMSO", "THZ1"))
#'
#' # Take the average of every gene by treatment group
#' by_group <- colmean(GSE161650_lc, group = Group)
#' by_group[1:5, ]
#'
colmean <- function(x, group) {
  stopifnot("x must be a numeric matrix" = is.numeric(x))
  stopifnot("group must be a factor variable" = is.factor(group))
  stopifnot("length of grouping factor must equal number of columns" = length(group) == ncol(x))

  t(rowsum(t(x), group = group, reorder = FALSE) / as.vector(table(group)))
}

#' Extract variable from an environment and remove that environment
#'
#' This function will extract all variables from the given environment and
#' assign them to the global environment and then optionally remove the
#' environment.
#' @param x Environment to extract variables from
#' @param remove Should the Environment be removed after extracting variables? Default TRUE
#' @return variables from x are assigned to the GlobalEnv after execution
#' @export
#' @examples
#' my_env <- new.env()
#' my_env$X <- 1000
#' my_env$Y <- 1:10
#'
#' # Extract X and Y to global environment and remove my_env
#' env2global(my_env)
#'
#' # Extract X and Y to global environment and keep my_env
#' env2global(my_env, remove = FALSE)
#'
env2global <- function(x, remove = TRUE) {
  if (!is.environment(x)) {
    stop("x must be an environment")
  }

  vars <- ls(x)
  for (v in vars) {
    assign(v, x[[v]], envir = .GlobalEnv)
  }

  if (remove) {
    rm(list = toString(substitute(x)), envir = .GlobalEnv)
  }
}

#' Perform simple imputation on rows of a matrix
#' 
#' Impute NA values for each row of a matrix.
#' @param x numeric matrix or data.frame that can be converted to a numeric matrix
#' @param fun Imputation function to apply to rows of the matrix. Default median
#' @export
#' @return Matrix with imputed values
#' @examples
#' # Create a matrix of values with NAs
#' X <- matrix(runif(25), 5, dimnames = list(paste0("CpG", 1:5), paste0("Sample", 1:5)))
#' X[sample.int(25, 5)] <- NA
#' X
#' 
#' # Impute missing values with row medians
#' impute(X)
#' 
#' # Impute missing values with row mins
#' impute(X, min)
impute <- function(x, fun = median) {
  if (is(x, "data.frame"))
    x <- data.matrix(x)
  
  imputed <- apply(x, 1, function(i) { 
    i[which(is.na(i))] <- do.call(fun, list(i, na.rm = TRUE))
    return(i)
  }, simplify = FALSE)
  m <- do.call(rbind, imputed)
  dimnames(m) <- dimnames(x)
  m
}

#' Remove low variance features from a matrix
#' 
#' This function removes the p lowest variance features from a matrix. The 
#' function expects features in rows and samples in columns (e.g. an expression
#' matrix).
#' @param x matrix, numeric data.frame, or SummarizedExperiment object
#' @param p proportion of low variance features to remove
#' @param assay if SummarizedExperiment, what assay to use. Default = "counts"
#' @return matrix with lowest variance features removed
#' @export
#' @examples
#' 
#' # Remove 80% lowest variance features
#' removed <- remove_var(GSE161650_lc, p = 0.8)
remove_var <- function(x, ...) UseMethod("remove_var")

remove_var.default <- function(x) {
  stop("Object of type ", class(x), " is not supported by this function")
}

#' @rdname remove_var
#' @export
remove_var.matrix <- function(x, p) {
  stopifnot("p must be between 0 and 1" = p > 0 & p < 1)
  message("Removing ", p*100, "% lowest variance features...")
  if (requireNamespace("matrixStats", quietly = TRUE)) {
    v <- matrixStats::rowVars(x, na.rm = TRUE, useNames = FALSE)
  } else if (requireNamespace("Rfast", quietly = TRUE)) {
    v <- Rfast::rowVars(x, std = FALSE, na.rm = TRUE)
  } else {
    v <- apply(x, 1, var, na.rm = TRUE)
  }
  o <- order(v, decreasing = TRUE)
  head(x[o, ], n = nrow(x) * (1 - p))
}

#' @rdname remove_var
#' @export
remove_var.data.frame <- function(x, p) {
  stopifnot("p must be between 0 and 1" = p > 0 & p < 1)
  m <- data.matrix(x)
  remove_var.matrix(m, p = p)
}

#' @rdname remove_var
#' @export
remove_var.SummarizedExperiment <- function(x, p, assay = "counts") {
  stopifnot("p must be between 0 and 1" = p > 0 & p < 1)
  m <- SummarizedExperiment::assay(x, assay)
  remove_var.matrix(m, p = p)
}
