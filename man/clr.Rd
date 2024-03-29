% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{clr}
\alias{clr}
\title{Centered Log-ratio transformation}
\usage{
clr(x, base = 2)
}
\arguments{
\item{x}{numeric vector}

\item{base}{integer base of the log function. Default 2.}
}
\description{
Calculate the centered log-ratio transformation of a vector
}
\examples{
x <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 1500)
clr(x)
}
