% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot-umap.R
\name{plot_umap}
\alias{plot_umap}
\title{Plot results of UMAP}
\usage{
plot_umap(
  df,
  x = "UMAP1",
  y = "UMAP2",
  colBy = NULL,
  shapeBy = NULL,
  pointSize = 3,
  pointAlpha = 1,
  hline = 0,
  vline = 0,
  hlineType = 2,
  vlineType = 2
)
}
\arguments{
\item{df}{data.frame of UMAP embeddings and metadata}

\item{x}{First UMAP component. Default "UMAP1"}

\item{y}{Second UMAP component. Default "UMAP2"}

\item{colBy}{Column name of data.frame to color points by. Default NULL}

\item{shapeBy}{Column name of data.frame to shape the points by. Default NULL}

\item{pointSize}{Size of the points. Default 3}

\item{pointAlpha}{Alpha level of the points. Default 1}

\item{hline}{y-position of horizontal line. Default 0}

\item{vline}{x-position of vertical line. Default 0}

\item{hlineType}{linetype of the horizontal line. Default 2}

\item{vlineType}{linetype of the vertical line. Default 2}
}
\value{
ggplot object
}
\description{
This function creates a simple plot of the data.frame returned by the 
\code{umap()} function.
}
\examples{
# Create metadata for plotting
metadata <- data.frame(row.names = colnames(GSE161650_lc))
metadata$Group <- rep(c("DMSO", "THZ1"), each = 3)

# PCA with PCAtools
p <- PCAtools::pca(GSE161650_lc, metadata, center = TRUE, scale = TRUE)
udata <- coriell::UMAP(p, n_neighbors = 2)
plot_umap(udata, colBy = "Group")
}
