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
#'   \item gene_col: Genes
#'   \item Experiment name(s): Differential expression calls for each experiment according to the user supplied cutoffs. 1 = "up-regulated", -1 = "down-regulated", and 0 = "unperturbed".
#'   \item n_de: The total number of studies in which the gene was differentially expressed.
#'   \item prop_de: The proportion of studies in which the gene was differentially expressed.
#'   \item sign_consistency: The number of studies in which a gene was up-regulated minus the number of studies in which a gene was down-regulated.
#'   \item vote: The final vote for the gene. Based on sign_consistency and meta_prop.  
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
    (prop_de >= meta_prop) & (sign_consistency < 0), "down",
    (prop_de >= meta_prop) & (sign_consistency > 0), "up",
    default = "unperturbed")
  
  # coerce to data.table
  result_dt <- data.table::as.data.table(
    data.frame(exp_mat, "n_de" = n_de, "prop_de" = prop_de, "sign_consistency" = sign_consistency, "vote" = vote),
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
