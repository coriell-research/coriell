#' Transform ages using Horvath's method
#' 
#' Transform age in years using the log transformation defined in Horvath (2013)
#' 
#' @details
#' \href{https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-10-r115}{Horvath, S. DNA methylation age of human tissues and cell types. 
#' Genome Biol 14, 3156 (2013). https://doi.org/10.1186/gb-2013-14-10-r115}
#' @param x numeric vector of ages
#' @param adult_age Age of the adult of the species. From the original manuscript, 
#' 20 for humans and 15 for chimpanzees. 
#' @return transformed ages
#' @export
#' @examples
#' 
#' plot(\(x) coriell::horvath_age(x), from=0, to=100, xlab="Age", ylab="Transformed Age")
horvath_age <- function(x, adult_age = 20) {
  ifelse(x <= adult_age, 
         log(x + 1) - log(adult_age + 1), 
         (x - adult_age) / (adult_age + 1)
  )
}