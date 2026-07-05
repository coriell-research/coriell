#' Calculate the Dispersion Separability Criterion (DSC)
#'
#' This function computes the \href{https://ibl.mdanderson.org/public-software/tcga-batch-effects/}{Dispersion Separability Criterion (DSC)}
#' metric which was defined by the MD Anderson Bioinformatics team for computing the strength of
#' batch effects in high-throughput sequencing data from TCGA.
#'
#' @param x A numeric matrix (features in rows, samples in columns)
#' @param g A vector containing the batch assignment for each sample
#' @return DSC metric
#' @export
#' @examples
#' set.seed(123)
#'
#' # Small DSC when samples drawn from same distribution
#' x1 <- matrix(rnorm(500 * 100), nrow = 500, ncol = 100)
#' batches <- rep(c("Batch_1", "Batch_2"), each = 50)
#' dsc(x = x1, g = batches)
#'
#' # Large DSC when samples have batch effect present
#' x2 <- x1
#' x2[, 51:100] <- x2[, 51:100] + 50
#' dsc(x2, batches)
#'
dsc <- function(x, g) {
  x <- as.matrix(x)
  g <- as.factor(g)
  N <- ncol(x)
  global_mean <- rowMeans(x, na.rm = TRUE)

  trace_sb <- 0
  trace_sw <- 0

  for (batch in levels(g)) {
    batch_data <- x[, g == batch, drop = FALSE]
    n_j <- ncol(batch_data)

    if (n_j == 0) {
      next
    }

    pi_j <- n_j / N
    batch_mean <- rowMeans(batch_data, na.rm = TRUE)
    trace_sb <- trace_sb +
      pi_j * sum((batch_mean - global_mean)^2, na.rm = TRUE)
    centered_data <- sweep(batch_data, 1, batch_mean, "-")
    trace_sw <- trace_sw + (1 / N) * sum(centered_data^2, na.rm = TRUE)
  }

  db <- sqrt(trace_sb)
  dw <- sqrt(trace_sw)

  return(db / dw)
}
