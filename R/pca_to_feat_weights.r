##############################################################################
## Functions to map factor weights of principle components back to features ##
##############################################################################

#' @title Map weights to features.
#'
#' @name map_weights_to_features
#'
#' @description This function maps factor weights of princiople components
#' back to individual features (e.g. genes). Performing PCA on omics data prior
#' to JDR can have several advantages. However, one of the strength of JDR
#' approaches is being able to identify which individual features contribute
#' most to the factorisation. In order achieve this when performing PCA first,
#' we must multiply the matrix of PCA loadings to the matrix of factor weights.
#' This only works with JDR methods that are linear. For details, check [cite my paper] # nolint
#'
#' @inheritParams run_jdr
#' @param fct_list A list of factorissations. Expects output of
#' \code{\link{run_jdr}}, or one of the individual functions.
#' @param pca_weights A pcaMethods output corresponding to the data
#' imputted for JDR. See \code{\link{prepare_data}}.
#' @param omics Names of the omics for which this should be done.
#' If NULL, it will be done on all omics
#'
#' @return
#'
#' @examples
#'
#' @export
#'

map_weights_to_features <- function(jdr_methods = c("MOFA", "JIVE", "RGCCA", "MCIA"),
                                    fct_list, pca_weights,
                                    omics = NULL) {
  # get MOFA weights
  if (is.element("MOFA", jdr_methods)) {
    mofa_wts <- .map_mofa_wts(fct_list$MOFA, pca_weights, omics)
  }

}

#' @name .map_mofa_wts
#'
#' @description Internal function that extracts factor weights from a trained
#' MOFA object and multiplies them with PCA loadings.
#'
#' @inheritParams map_weights_to_features
#' @param mofa_object A MOFA object containing a trained MOFA model.
#'
#' @return
#'
#' @import MOFA2 #maybe just import specific functions here

.map_mofa_wts <- function(mofa_object, pca_weights, omics) {
  # get weights
  weights <- get_weights(mofa_object)

  feat_wts <- list()

  if (is.null(omics)) {
    omics <- names(weights)
  }

  #substituting into the MOFA object
  metadata <- data.frame()
  dimensions <- c()

  for (j in omics) {
    loadings <- PCA_weights[[j]]@loadings
    #selecting only the PCs used in the mofa analysis
    loadings <- loadings[,1:nrow(weights[[j]])]
    # matrix multiplication between the PCA loadings and the factor weights
    mult <- loadings %*% weights[[j]]
    feat_wts[[j]] <- mult


    #change the expectations (i.e. actual weight values)
    MOFAobject@expectations$W[[j]] <- feat_wts[[j]]

    #get new feature metadata
    view_metadata <- data.frame(feature = rownames(feat_wts[[j]]), view = rep(j, nrow(feat_wts[[j]])))
    metadata <- rbind(metadata, view_metadata)
    
    #get new dimensions
    dimensions <- c(dimensions, nrow(feat_wts[[j]]))
    
    
    #replace feature metadata
    MOFAobject@features_metadata <- metadata
    
    #replace dimensions
    MOFAobject@dimensions$D <- dimensions
  }
  
  names(MOFAobject@dimensions$D) <- omics
  return(MOFAobject)
}
