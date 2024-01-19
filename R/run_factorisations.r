##################################
## Functions to run JDR methods ##
##################################

#' @title Run joint dimensionality reduction on multi-omics data.
#'
#' @name run_jdr
#'
#' @descripton Function for performing joint dimensionality reduction (JDR)
#' on omics data using one or more linear JDR methods.
#' Note: If only one JRD method is used, it is prefered to use the
#' individual functions:
#' \itemize{
#'  \item{\code{\link{run_mofa2}}}
#'  \item{\code{\link{run_jive}}}
#'  \item{\code{\link{run_rgcca}}}
#'  \item{\code{\link{run_mcia}}}
#' }
#'
#' @param omic_list A list of omic matrices.
#' Should be output of \code{link{prepare_data}}.
#' @param samples_overlap Logical. Whether a sample overlap was performed
#' between the omics during data preparation. See \code{\link{prepare_data}}.
#' If FALSE, omics will be filtered to only common samples for JDR methods that
#' require it. This may result in different sample sets being used with
#' different methods.
#' @param pca Logiacal. Whether PCA was performed on the data.
#' See \code{\link{prepare_data}}.
#' @param jdr_methods Character vector specifying one more more JDR methods
#' to be used. Should be from \code{c("MOFA", "JIVE", "RGCCA", "MCIA")}.
#' @param n_fct Integer. Number of factors for the factorisations. 
#' @param ... Any other parameters that can be passed to any used functions
#' from the JDR packages. See respective package documentation for further details. #nolint
#'
#' @return A list of factorisations. Each element is the factorisation
#' based on one JDR method.
#'
#' @examples
#'
#' @export
#'

run_jdr <- function(omic_list, samples_overlap = TRUE, pca = TRUE,
                    jdr_methods = c("MOFA", "JIVE", "RGCCA", "MCIA"),
                    n_fct = 5, ...) {
  # overlap samples if not already done
  if (!samples_overlap) {
    omic_fil <- .filter_omics(omic_list)
  }

  # run MOFA
  if (is.element("MOFA", jdr_methods)) {
    mofa_model <- .run_mofa2(omic_list, n_fct, ...)
  }

  # run JIVE
  if (is.element("JIVE", jdr_methods)) {
    jive_model <- .run_jive(omic_list, n_fct, ...)
  }

  # run RGCCA
  if (is.element("RGCCA", jdr_methods)) {
    if (!samples_overlap) {
      # run with omic_fil
      rgcca_model <- run_rgcca(omic_fil, n_fct, ...)
    } else {
      # run with omic_list
      rgcca_model <- run_rgcca(omic_list, n_fct, ...)
    }
  }

  # run MCIA
  if (is.element("MCIA", jdr_methods)) {
     if (!samples_overlap) {
      # run with omic_fil
    } else {
      #run with omic_list
    }
  }
}

#' @title Run JDR with MOFA2.
#'
#' @name run_mofa2
#'
#' @description Performs JDR on a list of omic matrices with MOFA2.
#'
#' @inheritParam run_jdr
#' @param ... Any other parameters that can be passed to MOFA2 functions.
#' See documentation of the MOFA2 package for details.
#'
#' @return A trained MOFA2 model.
#'
#' @export
#' @import MOFA2

run_mofa2 <- function(omic_list, n_fct, ...) {
  # make MOFA object
  mofa_object <- create_mofa(omic_list)

  # prepare MOFA object
  mofa_object <- prepare_mofa(mofa_object, ...)
  mofa_object@model # gotta figure out the best way to access the model opts and change the number of fct #nolint

  # run mofa
  mofa_model <- run_mofa(mofa_object, ...)

  return(mofa_model)
}

#' @title Run JDR with JIVE.
#'
#' @name run_jive
#'
#' @description Performs JRD on a list of omic matrices with JIVE.
#'
#' @inheritParam run_jdr
#'
#' @return A trained JIVE model.
#'
#' @export
#' @import r.jive

run_jive <- function(omic_list, n_fct, ...) {
  # run jive
  # look a little more into JIVE options and what they mean
  jive_model <- jive(omic_list, rankJ = n_fct, 
                      rankA = rep(n_fct, length(omic_list)), ...)

  return(jive_model)
}

#' @title Run JDR with RGCCA.
#'
#' @name run_rgcca
#'
#' @description Performs JDR on a list of omic matrices using RGCCA.
#'
#' @inheritParam run_jdr
#'
#' @return A trained RGCCA model.
#'
#' @export
#' @import RGCCA

run_rgcca <- function(omic_list, n_fct, ...) {
  # transpose omics
  omics_t <- lapply(omic_list, function(x) t(x))

  # run rgcca
  rgcca_model <- rgcca(omics_t, ncomp = rep(n_fct, length(omics_t)), ...)

  return(rgcca_model)
}