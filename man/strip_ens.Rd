% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{strip_ens}
\alias{strip_ens}
\title{Strip version IDs from ENSEMBL identifiers}
\usage{
strip_ens(x)
}
\arguments{
\item{x}{character vector containing ENSEMBL IDs}
}
\value{
character vector with trailing IDs removed
}
\description{
This function is simply a wrapper arond gsub for removing trailing IDs from 
ENSEMBL identifiers. All entries must start with 'ENS'. NAs are tolerated.
}
\examples{

ids <- c("ENSG0000001.12", "ENSMUSG00021.3", "ENST00000556.2")
strip_ens(ids)

}
