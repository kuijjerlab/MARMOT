################################################################################
#### Functions for performing factor association with survival and various #####
########################## other associated analyses ###########################

#' @title Calculate factor association with survival.
#'
#' @name calculate_surv_association
#'
#' @description This function takes as input computed factorisations from JDR
#' models (see \code{\link{run_jdr}}) and calculates the association to survival
#' of the factors. There are options for either univariate cox regression for
#' each factor, or multivariate cox regression for multiple factors.
#' @param factors Data frame containing factors to be associated with survival.
#' Samples should be rows and factors columns.
#' @param surv A data frame with survival information. Expects output from
#' \code{\link{prepare_surv}}.
#' @param univariate Logical, whether to run a univariate coxph for each factor
#' (TRUE) or a multivariate coxph for all together (FALSE). Default  is TRUE.
#'
#' @returns A list with either one multivariate coxph model or multiple
#' univariate coxph models.
#' @export

calculate_surv_association <- function(factors, surv, univariate = TRUE) {

  # make sure samples are the same
  samples <- surv$sampleID
  samples<-intersect(samples, rownames(factors))

  #make sure factors are named
  if (is.null(colnames(factors))) {
    colnames(factors) <- paste0("Factor", seq(1, ncol(factors)))
  }

  factors <- factors[samples, ]
  surv <- surv[which(surv$sampleID %in% samples), ]

  survival_object <- Surv(surv$days_to_event, surv$vital_status)

  # run either one multivariate or multiple univariate coxph models
  cox <- list()

  if (by_factor) {
    cox <- lapply(colnames(factors), function(colname) {
      temp <- coxph(survival_object ~ factors[, colname])
      names(temp$coefficients) <- colname
      return(temp)
    })
    names(cox) <- colnames(factors)
  } else {
    cox[[1]] <- coxph(survival_object ~ factors)
  }

  return(cox)
}