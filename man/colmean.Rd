% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{colmean}
\alias{colmean}
\title{Give means of rows of matrix based on column grouping variable}
\usage{
colmean(x, group, na.rm = FALSE)
}
\arguments{
\item{x}{variable (gene) by sample numeric matrix}

\item{group}{Factor specifying the grouping level to by averaged}

\item{na.rm}{Logical (TRUE or FALSE). Should NA (including NaN) values be discarded?}
}
\value{
variable x nLevels(group) matrix
}
\description{
Inspired by \code{edgeR::sumTechReps} and \code{base::rowsum()}, this function
takes the average of the values in each group given by the group argument
for each row of the data matrix.
}
\examples{
# Specify the Group levels
Group <- gl(n = 2, k = 3, labels = c("DMSO", "THZ1"))

# Take the average of every gene by treatment group
by_group <- colmean(GSE161650_lc, group = Group)
by_group[1:5, ]

}
