###############################################################################
## Functions to perform differential analysis on omics based on JDR factors ###
###############################################################################

#' @title Differential analysis based on JDR factors.
#'
#' @name differenitial_analysis
#'
#' @description This function performs differential analysis on omics data. It
#' defines the groups based on given JDR factors associated with a given
#' clinical feature.
#'
#' @param omic A matrix of omics data. Columns should be samples and features
#' should be rows.
#' @param factor A dataframe with the factor based on which to determine groups.
#' @param clin A dataframe of clinical data.
#' @param feat Character string with the clinical feature with which the factor
#' of interest is associated. Must be a column name in \code{clin}.
#' @param covariates Optional. A character vector with names of covariates for
#' which to correct. Only applies if limma is used.
#' @param limma Whether limma should be used for the differential analysis.
#' If FALSE, a wilcoxon signed rank test will be used instead. Default is TRUE.

differential_analysis <- function(omic, factor, clin, feat, covariates = NULL,
                                  limma = TRUE) {

}