#' Perform vote counting on differential expression results
#'
#' Perform vote-counting on a list of differential expression results. This method 
#' uses user defined p-value and log fold-change cutoffs to determine the number
#' of times a gene is differentially expressed across multiple studies. A gene is
#' then considered 'meta-differentially expressed' if it is differentially expressed
#' in the user defined proportion of studies.
#' 
#' There are certain caveats to consider when using the vote counting method. 1.
#' combining FDR and logFC cutoffs after the fact destroys FDR correction. See
#' \href{https://academic.oup.com/bioinformatics/article/25/6/765/251641}{McCarthy & Smyth, 2009} for more details. 
#' 2. The direction of expression is not considered when the final vote is cast. For example, if 
#' geneA is up-regulated in 4 studies but down-regulated in 1 then the final vote for 
#' geneA will be 'up-regulated'.
#' @param exp_list list() of data.frames of differential expression results.
#' @param lfc_col Column name in expression data.frames that contains the log fold-change estimates. default = "logFC". 
#' @param pval_col Column name in expression data.frames that contains the significance values. default = "FDR".
#' @param gene_col Column name in the expression data.frames that contains the gene IDs. default = "feature_id".
#' @param lfc Numeric. Log fold-change value used as cutoff for determining differential expression in each study. default = 0. 
#' @param pval Numeric. Significance value used as cutoff for determining differential expression in each study
#' @param meta_prop Numeric. Proportion of studies a gene must be differentially expressed in to be considered as 'meta-differentially expressed'. default = 0.5.
#' @param plot Logical. Display a metavolcano ggplot2 plot of the results. default = FALSE.
#' @param all_common Logical. Use only genes present in all experiments. default = FALSE.
#' @return data.table of meta-expression results with the following columns:
#' \itemize{
#'   \item gene_col: Genes included in the final voting procedure. 
#'   \item n_de: The total number of studies in which the gene was differentially expressed (either up or down-regulated). 
#'   \item prop_de: The proportion of studies in which the gene was differentially expressed (either up or down-regulated).
#'   \item sign_consistency: The number of studies in which a gene was up-regulated minus the number of studies in which a gene was down-regulated.
#'   \item vote: The final vote for the gene. based on \code{sign(sign consistency) & abs(sign consistency) >= meta_prop * N studies}.
#' }
#' @import data.table
#' @export
meta_vote <- function(exp_list, lfc_col = "logFC", pval_col = "FDR", 
                      gene_col = "feature_id", lfc = 0, pval = 0.05, 
                      meta_prop = 0.5, all_common = FALSE, plot = FALSE) {
  N_studies <- length(exp_list)
  names(exp_list) <- if (is.null(names(exp_list))) paste0("Experiment_", 1:N_studies) else names(exp_list)
  exp_dt <- data.table::rbindlist(exp_list, idcol = "experiment")
  
  # label differential expression based on user input
  exp_dt[, direction := data.table::fcase(
    abs(get(lfc_col)) > lfc & get(pval_col) < pval & get(lfc_col) > 0, 1L,
    abs(get(lfc_col)) > lfc & get(pval_col) < pval & get(lfc_col) < 0, -1L,
    default = 0L)]
  
  # Cast wider into gene x experiment matrix
  exp_mat <- as.matrix(
    data.table::dcast(
      data = exp_dt[, .(experiment, id = get(gene_col), direction)], 
      formula = id ~ experiment, 
      value.var = "direction", 
      fill = NA), 
    rownames = "id"
  )
  
  if (all_common) {
    exp_mat <- na.omit(exp_mat)
  }
  
  # calculate voting stats
  sign_consistency <- rowSums(exp_mat, na.rm = TRUE)
  n_de <- rowSums(abs(exp_mat), na.rm = TRUE)
  prop_de <- n_de / N_studies
  vote <- data.table::fcase(
    (abs(sign_consistency) >= meta_prop * N_studies) & (sign_consistency < 0), "down",
    (abs(sign_consistency) >= meta_prop * N_studies) & (sign_consistency > 0), "up",
    default = "unperturbed")
  
  # coerce to data.table
  result_dt <- data.table::as.data.table(
    data.frame("n_de" = n_de, "prop_de" = prop_de, "sign_consistency" = sign_consistency, "vote" = vote),
    keep.rownames = gene_col)
  
  if (plot) {
    print(plot_metavolcano(result_dt, gene_col))
  }
  
  return(result_dt)
}

#' Create Meta-volcano plot
#' 
#' Specialty function for creating metavolcano plots.
#' 
#' @param meta_df data.table of meta_vote results.
#' @param label_col character string of column with gene ids
#' @export
#' @return ggplot2 object
plot_metavolcano <- function(meta_df, label_col = "feature_id") {
  ggplot2::ggplot(meta_df, ggplot2::aes_string(x = "sign_consistency", y = "n_de", color = "vote", label = label_col)) +
    ggplot2::geom_jitter(alpha = 0.8) +
    ggplot2::scale_color_manual(values = c("up" = "firebrick", "down" = "steelblue", "unperturbed" = "grey80")) +
    ggplot2::theme_light() +
    ggplot2::labs(x = "Sign Consistency",
                  y = "Number of times differentially expressed",
                  color = "Vote") +
    ggplot2::theme(legend.position = "bottom")
}

#' Perform p-value combination on differential expression results
#' 
#' This function aims to combine p-values and log fold-change values into a single
#' combined p-value and logFC value for each gene across all studies. Several methods 
#' for combining p-values are included. The useer can also pass an arbitrary function
#' to the \code{lfc_fun} argument used to summarize the logFC values across studies. 
#' 
#' The following p-value combination methods are available:
#' \itemize{
#'  \item fisher: Uses the sum of the logarithms for the p-values. Sensitive to very small p-values therefore a single significant study can lead to a very small combined p-value.
#'  \item pearson: Similar to Fisher's method but sensitive to large p-values; therefore more false negatives are obtained.
#'  \item stouffer: (not yet implemented) Recommended when weights for each study can be calculated. 
#'  \item tippet: Uses the minimum of the p-values across all studies. Recommended when the aim is to discard genes.
#'  \item wilkinson: Uses the maximum of the p-values across all studies. Recommended when the aim is to identify the most robust genes.
#' }
#'
#' @param exp_list list() of data.frames of differential expression results.
#' @param lfc_col Column name in the expression data.frames that contains log fold-change values. default = "logFC".
#' @param pval_col Column name in expression data.frames that contains the significance values. default = "FDR".
#' @param gene_col Column name in the expression data.frames that contains the gene IDs. default = "feature_id".
#' @param all_common Logical. Use only genes present in all experiments. default = FALSE.
#' @param lfc_fun Function for summarizing logFC values across studies. 
#' @param method Method for combining p-values. One of c("fisher", "pearson", "stouffer", "tippet", "wilkinson"). default "fisher".
#' @param plot Logical. Display a volcano plot using the combined values. default = FALSE
#' @param plot_lfc Numeric. LogFC cutoff used when plotting. default = 0
#' @param plot_pval Numeric. P-value cutoff used when plotting. default = 0.05.
#' @param plot_labels Logical. Display labels for significant genes on plot. default = FALSE
#' @param plot_count Logical. Display count of up/down genes on plot. default = TRUE
#' @param ... Additional arguments passed to lfc_fun. 
#' @import data.table
#' @export
#' @return data.table with columns containing the combined logFC and p-values.
meta_pcombine <- function(exp_list, lfc_col = "logFC", pval_col = "FDR", 
                          gene_col = "feature_id", all_common = FALSE, 
                          lfc_fun = mean, method = c("fisher", "pearson", "stouffer", "tippet", "wilkinson"),
                          plot = FALSE, plot_lfc = 0, plot_pval = 0.05, plot_labels = FALSE, plot_count = TRUE,
                          ...) {
  exp_dt <- data.table::rbindlist(exp_list, idcol = "experiment")
  
  pval_mat <- as.matrix(
    data.table::dcast(
      data = exp_dt[, .(experiment, id = get(gene_col), pval = get(pval_col))],
      formula = id ~ experiment,
      value.var = "pval",
      fill = NA), 
    rownames = "id")
  
  lfc_mat <- as.matrix(
    data.table::dcast(
      data = exp_dt[, .(experiment, id = get(gene_col), lfc = get(lfc_col))],
      formula = id ~ experiment,
      value.var = "lfc",
      fill = NA),
    rownames = "id"
  )
  
  if (all_common) {
    pval_mat <- na.omit(pval_mat)
    lfc_mat <- na.omit(lfc_mat)
  }
  
  # Combine p-values and logFCs
  p_method <- match.arg(method)
  P_FUN <- switch(p_method,
    fisher = meta_fisher,
    pearson = meta_pearson,
    stouffer = stop("Not yet implemented"),
    tippet = meta_tippet,
    wilkinson = meta_wilkinson
  )
  meta_p <- apply(pval_mat, 1, P_FUN)
  meta_lfc <- apply(lfc_mat, 1, lfc_fun, ... = ...)
  
  result_dt <- data.table::as.data.table(
    data.frame(meta_lfc, meta_p), 
    keep.rownames = gene_col)
  
  if (plot) {
    print(coriell::plot_volcano(
      result_dt, x = meta_lfc, y = meta_p, fdr = plot_pval, lab = get(gene_col), 
      lfc = plot_lfc, label_sig = plot_labels, annotate_counts = plot_count) +
        ggplot2::labs(x = "Combined logFC",
                      y = paste0("-log10(meta p-value) [", method, "]"),
                      subtitle = NULL,
                      caption = paste("Meta p-value cutoff:", plot_pval, "\nMeta logFC cutoff:", plot_lfc)) +
        ggplot2::theme(legend.position = "bottom")
      )
  }
  return(result_dt)
}

#' Fisher's method for combining p-values
#' 
#' @param x Numeric vector of p-values
#' @noRd
#' @return Probability of the Fisher test statistic under the chisq distribution
meta_fisher <- function(x) {
  test_stat <- -2 * (sum(log(x), na.rm  = TRUE))
  pval <- pchisq(q = test_stat, df = 2 * length(x), lower.tail = FALSE)
  return(pval)
}

#' Pearson's method for combining p-values
#' 
#' @param x Numeric vector of p-values
#' @noRd
#' @return Probability of the Pearson test statistic under the chisq distribution
meta_pearson <- function(x) {
  test_stat <- -2 * (sum(log(1 - x), na.rm = TRUE))
  pval <- pchisq(q = test_stat, df = 2 * length(x), lower.tail = FALSE)
  return(pval)
}

#' Tippet's method for combining p-values
#' 
#' @param x Numeric vector of p-values
#' @noRd
#' @return Probability of the Tippet test statistic under the beta distribution
meta_tippet <- function(x) {
  test_stat <- min(x, na.rm = TRUE)
  pval <- pbeta(q = test_stat, shape1 = 1, shape2 = length(x))
  return(pval)
}

#' Wilkinson's method for combining p-values
#' 
#' @param x Numeric vector of p-values
#' @noRd
#' @return Probability of the Wilkinson test statistic under the beta distribution
meta_wilkinson <- function(x) {
  test_stat <- max(x, na.rm = TRUE)
  pval <- pbeta(q = test_stat, shape1 = 1, shape2 = length(x))
  return(pval)
}
