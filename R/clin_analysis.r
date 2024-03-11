################################################################################
###### Functions for performing factor association with clinical features ######
################################################################################

#' @name clin_association
#'
#' @description Will test if factors from a MOFA model are significantly
#' associated with clinical features. Will perform a Wilcoxon rank-sum
#' ("Mann-Whitney") test for binary features (i.e. TRUE/FALSE, yes/no etc)
#' and a kruskal rank-sum test otherwise. Make sure the samples_metadata of the
#' MOFA object contains the clin data. Make sure any missing values are true NAs
#' (not "NA" character string). NAs will be excluded.
#'
#' @inheritParams surv_association
#' @param clin Data frame with clinical data where rows are samples and columns
#' are clinical features.
#' @param clin_feat A character vector with the names of clinical features to be
#' tested.
#' @param which_fct Optional. Character vector indicating which factors to
#' perform the association for. If NULL, all factors will be tested.
#'
#' @return Results of the test
#' @export

clin_associaton <- function(factors, clin_feat, factors = NULL){
  # sanity checks
  clin <- samples_metadata(MOFA_model)
  stopifnot("Elements in clin_feat do not match the sample metadata. Check if they are spelled correctly." = clin_feat %in% colnames(clin))
  
  # extract factors
  Z <- get_factors(MOFA_model)[[1]]
  if (!is.null(factors)) {
    Z <- Z[factors]
  }
  
  # innitialise results data frame
  results_df <- data.frame()
  
  # perform test
  for (i in clin_feat){
    
    feat <- factor(clin[,i])
    if(length(levels(feat)) < 2){
      print(paste("Feature",i,"will be discarded; All observations are the same."))
      next
    }
    else if(length(levels(feat)) == 2) {
      test <- "wilcox"
      results <- apply(Z, MARGIN=2, function(x) wilcox.test(x~feat, na.action = "na.exclude"))
    } else {
      test <- "kruskal"
      results <- apply(Z, MARGIN=2, function(x) kruskal.test(x~feat, na.action = "na.exclude"))
    }
    
    # store all relevant info in a data frame
    as_df <- as.data.frame(colnames(Z))
    colnames(as_df) <- "Factor"
    as_df$feat <- i
    as_df$test <- test
    as_df$pvalue <- unlist(lapply(results, function(x) x$p.value))
    as_df$padj <- p.adjust(as_df$pval, method = "BH")
    as_df$logp <- -log10(as_df$padj)
    
    results_df <- rbind(results_df, as_df)
  }
  
  return(results_df)
}