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
#' @inheritParams surv_compare
#' @param clin Data frame with clinical data where rows are samples and columns
#' are clinical features.
#' @param clin_feat A character vector with the names of clinical features to be
#' tested.
#' @param which_fct Optional. Character vector indicating which factors to
#' perform the association for. If NULL, all factors will be tested.
#'
#' @return Results of the test
#' @export

clin_associaton <- function(factors, clin, clin_feat, which_fct = NULL,
                            p_adjust = TRUE, method = "BH", logtrans = TRUE) {

  # sanity checks
  # check feature names exist
  .check_names(clin_feat, colnames(clin), err_msg = "feature names exist in the 
  clinical data")

  # check factor names exist
  if (!is.null(which_fct)) {
    .check_names(which_fct, colnames(factors), err_msg = "Please make sure 
    factor names exist and are spelled correctly")
  }

  # select features of interest
  clin2 <- clin[, clin_feat]

  # check clinical features of interest are not uniform
  clin2 <- .check_variance(clin2)

  # innitialise results data frame
  results_df <- data.frame()

  # perform test
  for (i in clin_feat){
    
    if(length(levels(feat)) == 2) {
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

#' @title Perform a wilcoxon rank sum or kruskal-Wallis test on a factor.
#'
#' @name .perform_test
#'
#' @description
#'
#' @inheritParams surv_association
#' @param feat Vector containing the feature of interest.
#' @param feat_name Optional. Character string specifying the name of the
#' feature.
#'
#' @return Test results.
#'
#' @importFrom dplyr everything summarie across

.perform_test <- function(feat, factors, feat_name = NULL) {
  # make sure the feature vector is a factor
  feat <- as.factor(feat)

  # determine the appropriate test
  if (length(levels(feat)) == 2) {
    test <- "wilcox"
    results <- factors %>%
      summarise(across(everything(), ~ wilcox.test(. ~ feat,
                                                   na.action = "na.exclude")))
  } else {
    test <- "kruskal"
    results <- factors %>%
      summarise(across(everything(), ~ kruskal.test(. ~ feat,
                                                    na.action = "na.exclude")))
  }


}