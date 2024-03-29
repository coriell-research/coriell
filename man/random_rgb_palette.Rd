% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/colors.R
\name{random_rgb_palette}
\alias{random_rgb_palette}
\title{Generate random RGB color palette}
\usage{
random_rgb_palette(n, alpha = 1)
}
\arguments{
\item{n}{numeric. Number of random colors to generate}

\item{alpha}{numeric. Transparency level of color palette (0-1). Default 1.0}
}
\value{
vector of hex values for RGB colors of length n
}
\description{
This function will randomly generate N colors from the rgb colorspace. It
makes no attempt at ensuring the colors are distinct or well separated.
}
