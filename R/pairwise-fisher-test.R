#' Perform pairwise fisher.tests on an input matrix
#' 
#' This function will perform a \code{fisher.test()} for all columns relative
#' to the given reference column of the input 2 x m matrix and return a 
#' data.frame of the results.
#' 
#' @details The function performs a fisher test for each test column relative to
#' the reference column of the input matrix. The input matrix must be a 2 x m 
#' matrix. It is the user's responsibility to ensure the levels of the rows align
#' with the desired odds ratio calculation. e.g. R will calculate the odds ratio
#' as \code{(x[1, "test_col"] / x[2, "test_col"]) / (x[1, "ref_col"] / x[2, "ref_col"])}. 
#' It is the user's responsibility to make sure \code{x[1, ]} and \code{x[2, ]} 
#' are in the desired order.
#' @export
#' @param x numeric 2 x m data.frame, matrix, or table
#' @param ref column name of the reference column of the input matrix
#' @return data.frame of results for each test.
#' @examples
#' # From chisq.test docs:
#' M <- as.table(rbind(c(762, 327, 468), c(484, 239, 477)))
#' dimnames(M) <- list(gender = c("F", "M"),
#'                     party = c("Democrat","Independent", "Republican"))
#' pairwise_fisher_test(M, ref = "Independent")
pairwise_fisher_test <-function(x, ref) {
  stopifnot("data must have exactly two rows" = nrow(x) == 2)
  stopifnot("ref not in colnames(x)" = ref %in% colnames(x))
  
  test_cols <- setdiff(colnames(x), ref) 
  or <- vector("numeric", length(test_cols))
  ci_low <- vector("numeric", length(test_cols))
  ci_high <- vector("numeric", length(test_cols))
  p_val <- vector("numeric", length(test_cols))
  comp <- vector("character", length(test_cols))
  
  for (i in seq_along(test_cols)) {
    test_col <- test_cols[i]
    comp[[i]] <- paste(test_col, "vs", ref)
    ft <- fisher.test(x[, c(test_col, ref)])
    
    or[[i]] <- ft$estimate
    ci_low[[i]] <- ft$conf.int[1]
    ci_high[[i]] <- ft$conf.int[2]
    p_val[[i]] <- ft$p.value
  }

  data.frame(
    comparison = comp,
    odds_ratio = or,
    ci_low = ci_low,
    ci_high = ci_high,
    p_value = p_val,
    p_value_adj = p.adjust(p_val, method = "bonferroni")
  )
}
