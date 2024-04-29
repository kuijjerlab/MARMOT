################################################################################
## Functions relating to testing the contribution of each view to a MOFA model #
################################################################################

#' @name fct_corr
#'
#' @description Perform pearson and spearman correlation between two sets of JDR
#' factors. Primarily intended for comparing different JDR models on the same
#' data.
#'
#' @inheritParams view_exclusion
#' @param set1 A matrix with the first set of factors.
#' @param set2 A matrix with the second set of factors.
#'
#' @returns A list containing the factor correlation matrices.
#'
#' @export

fct_corr <- function(set1, set2) {
  # sanity checks
  # check that the number of features is the same in the two sets
  if (nrow(set1) != nrow(set2)) {
    warning("Differing numbers of features (rows) between sets of factors. Results may be unreliable.") # nolint
  }

  # perform pairwise correlations of each factor in one set with the others
  corr_pears <- apply(set1, 2, function(x) {
    apply(set2, 2, function(y) {
      cor(x, y, method = "pearson")
    })
  })

  corr_spear <- apply(set1, 2, function(x) {
    apply(set2, 2, function(y) {
      cor(x, y, method = "spearman")
    })
  })

  corr_list <- list(corr_pears, corr_spear)

  return(corr_list)
}

#' @name plot_fct_corr
#'
#' @description Plots heatmaps of factor correlations.
#'
#' @param corr_list A data frame with the correlation values. Expects output
#' of \code{\link{fct_corr}}.
#' @param title Character string indicating a title for the plot.
#'
#' @returns A list of two ggplots. One for pearson and one for spearman.
#'
#' @export

plot_fct_corr <- function(corr_list, title) {
  
}