##############################################################################
## Functions to map factor weights of principle components back to features ##
##############################################################################

#' @name map_wts
#'
#' @description Function to map back to original feature weights when using
#' This function multiplies PCA loadings with JDR fctor weights to give a
#' matrix of the original features (e.g. genes) mapped to the JDR factors.
#'
#' @param fct_weights A matrix of JDR factor weights for an omic.
#' @param pca_weights A pcaMethods object for the same omic.
#'
#' @returns A matrix of features mapped to factors.
#' @export

map_wts <- function(fct_weights, pca_weights) {
  # get pca loadings
  loadings <- pca_weights@loadings

  #selecting only the PCs used in the mofa analysis
  loadings <- loadings[, seq_len(nrow(fct_weights))]

  # matrix multiplication between the PCA loadings and the factor weights
  feat_wts <- loadings %*% fct_weights

  return(feat_wts)
}
