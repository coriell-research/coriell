#' Flag outliers by IQR
#'
#' Flag outliers based on the interquartile range
#'
#' @param x numeric vector
#' @param threshold threshold for the upper and lower limits, A.K.A. Tukey's fences
#' @param direction return TRUE if the value is above or below the outlier cutoff. 
#' default "both", samples above and below the threshold are called outliers.
#'
#' @returns boolean vector indicating which values of the input vector are flagged as outliers
#' @export
#'
#' @examples
#' 
#' x <- c(1, 1, 2, 2, 4, 6, 11)
#' 
#' outliers_by_iqr(x)
#' 
#' # Using direction="low" disregards outliers above threshold, for example
#' outliers_by_iqr(x, direction="low")
#' 
outliers_by_iqr <- function(x, threshold = 1.5, direction = c("both", "low", "high")) {
  d <- match.arg(direction)

  iqr <- IQR(x, na.rm = TRUE)
  q1 <- quantile(x, probs = 0.25, na.rm = TRUE)
  q3 <- quantile(x, probs = 0.75, na.rm = TRUE)
  lower <- q1 - (iqr * threshold)
  upper <- q3 + (iqr * threshold)

  result <- switch(d,
    "both" = (x > upper) | (x < lower),
    "low" = x < lower,
    "high" = x > upper
  )

  return(result)
}


#' Flag outliers by MAD
#' 
#' Flag outliers based on the median absolute deviation
#'
#' @param x numeric vector
#' @param threshold threshold for number of MADs used to determine outlier. default 3
#' @param direction return TRUE if the value is above or below the outlier cutoff. 
#' default "both", samples above and below the threshold are called outliers.
#'
#' @returns boolean vector indicating which values of the input vector are flagged as outliers
#' @export
#'
#' @examples
#' 
#' x <- c(1, 1, 2, 2, 4, 6, 11)
#' 
#' outliers_by_mad(x)
#' 
#' # Using direction="low" disregards outliers above threshold, for example
#' outliers_by_mad(x, direction="low")
#' 
outliers_by_mad <- function(x, threshold = 3, direction = c("both", "low", "high")) {
  d <- match.arg(direction)

  MAD <- mad(x, na.rm = TRUE)
  cmedian <- median(x, na.rm = TRUE)
  lower <- cmedian - (threshold * MAD)
  upper <- cmedian + (threshold * MAD)

  result <- switch(d,
    "both" = (x > upper) | (x < lower),
    "low" = x < lower,
    "high" = x > upper
  )

  return(result)
}


#' Flag outliers by z-score
#' 
#' Flag outliers that exceed a given z-score value.
#'
#' @param x numeric vector
#' @param threshold threshold for number of SDs used to determine outlier. default 3
#' @param direction return TRUE if the value is above or below the outlier cutoff. 
#' default "both", samples above and below the threshold are called outliers.
#'
#' @returns boolean vector indicating which values of the input vector are flagged as outliers
#' @export
#'
#' @examples
#' x <- c(rnorm(10), 100)
#' 
#' outliers_by_z(x)
#' 
#' # Using direction="low" disregards outliers above threshold, for example
#' outliers_by_z(x, direction="low")
#' 
outliers_by_z <- function(x, threshold = 3, direction = c("both", "low", "high")) {
  d <- match.arg(direction)
  
  z <- as.numeric(scale(x, center = TRUE, scale = TRUE))
  
  result <- switch(d,
                   "both" = (z > threshold) | (z < -threshold),
                   "low" = z < -threshold,
                   "high" = z > threshold
  )
  
  return(result)
}
