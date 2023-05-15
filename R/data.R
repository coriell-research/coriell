#' Differential expression results from GSE161650
#'
#' A subset of differential expression results comp THZ1 (n=3) vs DMSO (n=3)
#' samples after analysis with edgeR. See Details for citation.
#'
#' Jiang B, Gao Y, Che J, Lu W et al. Discovery and resistance mechanism of a selective CDK12 degrader.
#' Nat Chem Biol 2021 Jun;17(6):675-683. PMID: 33753926
#'
#' @format
#' A data frame with 11,002 rows and 6 columns:
#' \describe{
#'   \item{feature_id}{Gene Symbol}
#'   \item{logFC}{log2 fold-change THZ1 vs DMSO}
#'   \item{unshrunk.logFC}{unshrunk log2 fold-change THZ1 vs DMSO}
#'   \item{logCPM}{Average logCPM across all samples}
#'   \item{PValue}{Uncorrected PValue}
#'   \item{FDR}{FDR Adjusted PValue}
#' }
#' @source <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161650>
"GSE161650_de"

#' Normalized log2 counts from GSE161650
#'
#' A subset of log2 expression counts containing THZ1 (n=3) vs DMSO (n=3)
#' samples after TMM normalization and log2 transformation. See Details for
#' citation
#'
#' Jiang B, Gao Y, Che J, Lu W et al. Discovery and resistance mechanism of a selective CDK12 degrader.
#' Nat Chem Biol 2021 Jun;17(6):675-683. PMID: 33753926
#' @format
#' A matrix with 11,002 rows (Gene Symbols) and 6 columns (Samples).
#' @source <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161650>
"GSE161650_lc"
