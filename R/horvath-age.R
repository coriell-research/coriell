#' Transform ages using Horvath's method
#'
#' Transform age in years using the log transformation defined in Horvath (2013).
#' This function can also be used to convert transformed ages back to the
#' original scale.
#'
#' @details
#' \href{https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-10-r115}{Horvath, S. DNA methylation age of human tissues and cell types.
#' Genome Biol 14, 3156 (2013). https://doi.org/10.1186/gb-2013-14-10-r115}
#' @param x numeric vector of ages
#' @param adult_age Age of the adult of the species. From the original manuscript,
#' 20 for humans and 15 for chimpanzees.
#' @param inverse If TRUE, return transformed ages to the original scale
#' @return transformed ages
#' @export
#' @examples
#'
#' age <- 0:100
#' transformed <- horvath_age(age)
#' plot(age, transformed, xlab="Age", ylab="Transformed Age", type="l")
#'
#' # Conversion back to original scale
#' all.equal(age, horvath_age(transformed, inverse = TRUE))
horvath_age <- function(x, adult_age = 20, inverse = FALSE) {
  if (isTRUE(inverse)) {
    return(ifelse(
      x <= 0,
      exp(x) * (adult_age + 1) - 1,
      x * (adult_age + 1) + adult_age
    ))
  }

  ifelse(
    x <= adult_age,
    log1p(x) - log1p(adult_age),
    (x - adult_age) / (adult_age + 1)
  )
}
