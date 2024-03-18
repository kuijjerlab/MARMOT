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
#' @param sample_label Optional. A character string indicating the name of the
#' column with sample IDs in the clinical data. If NULL, rownames will be
#' assumed to be sample IDs.
#'
#' @returns Results of the test
#' @export
#' @importFrom magrittr %>%

clin_associaton <- function(factors, clin, clin_feat, which_fct = NULL,
                            sample_label = NULL, method = "BH") {

  # sanity checks
  # check feature names exist
  .check_names(clin_feat, colnames(clin),
               err_msg = "feature names exist in the clinical data")

  # check factor names exist
  if (!is.null(which_fct)) {
    .check_names(which_fct, colnames(factors), err_msg = "Please make sure factor names exist and are spelled correctly") # nolint
  }

  # check sample label exists
  .check_names(sample_label, colnames(clin), err_msg = "Please make sure sample_label exists and is spelled correctly.") # nolint

  # make sure samples overlap
  if (is.null(sample_label)) {
    smpl <- intersect(rownames(clin), rownames(factors))
    clin2 <- clin[smpl, ]
  } else {
    smpl <- intersect(clin[, sample_label], rownames(factors))
    clin2 <- clin[which(clin[, sample_label] %in% smpl), ]
  }

  factors2 <- factors[smpl, ]

  # select features of interest
  clin2 <- clin2[, clin_feat]

  # check clinical features of interest are not uniform
  clin2 <- .check_variance(clin2)

  # make sure missing values are NAs
  clin2 <- dplyr::mutate_all(clin2, ~ ifelse(. == "", NA, .))

  # make sure all NAs are true NAs and not NA strings
  clin2 <- dplyr::mutate_all(clin2, ~ ifelse(. == "NA", NA, .))

  # perform tests for each feature
  # still have to test that this works right
  result_list <- purrr::pmap(list(clin2, clin_feat), ~ .perform_test(..1, factors2, feat_name = ..2)) # nolint: line_length_linter.

  # convert to a data frame
  result_df <- dplyr::bind_rows(result_list)

  # add column with factor name
  result_df$Factor <- sub("\\.{3}.*", "", rownames(result_df))

  # remove row names
  rownames(result_df) <- NULL

  # add logtrans
  result_df$logp <- -log10(result_df$padj)

  return(result_df)
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
#' @returns Test results.
#'
#' @importFrom dplyr everything summarise across
#' @noRd

.perform_test <- function(feat, factors, feat_name = NULL, method = "BH") {
  
  # determine the appropriate test
  if (length(levels(feat)) == 2) {
    test <- "wilcox"
    results <- data.frame(
      test = test,
      feat_name = feat_name,
      pval = sapply(names(factors), function(x) wilcox.test(get(x) ~ feat, data = factors, na.action = "na.exclude")$p.value) # nolint: line_length_linter.
    )
  } else {
    test <- "kruskal"
    results <- data.frame(
      test = test,
      feat_name = feat_name,
      pval = sapply(names(factors), function(x) kruskal.test(get(x) ~ feat, data = factors, na.action = "na.exclude")$p.value) # nolint: line_length_linter.
    )
  }

  # add the test that was used
  results$test <- test

  # add name of tested feature
  results$feat <- feat_name

  # adjust pvals
  results$padj <- p.adjust(results$pval, method = method)

  return(results)
}