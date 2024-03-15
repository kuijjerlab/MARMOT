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
#' @inheritParams surv_association
#' @param omic A matrix of omics data. Columns should be samples and features
#' should be rows.
#' @param clin Optional. A dataframe with clinical features for which to correct
#' if using limma.
#' @param factor A dataframe with the factor based on which to determine groups.
#' @param covariates Optional. A character vector with names of covariates for
#' which to correct. Only applies if limma is used.
#' @param limma Whether limma should be used for the differential analysis.
#' If FALSE, a wilcoxon signed rank test will be used instead. Default is TRUE.
#' @param minprop Numeric between c(0,1), indicating the minimum proportion of
#' samples per group.
#'
#' @export

differential_analysis <- function(omic, factor, surv, clin = NULL,
                                  covariates = NULL, limma = TRUE,
                                  minprop = 0.1) {
  # sanity checks
  if (!is.null(covariates)) {
    # check that clinical data is provided if covariates are provided
    if (is.null(clin)) {
      stop("Please provide clinical data for covariate correction.")
    }
    # check that covariates exist
    .check_names(covariates, colnames(clin),
                 err_msg = "covariate names exist in the clinical data")

    # check that limma is selected as the method
    if (!limma) {
      stop("Limma must be used for covariate correction.")
    }
  }

  # check that sample names for survival and clinical are the same
  

  # only take samples for which there is clin info
  samples <- surv[, 1]
  samples <- intersect(samples, rownames(factor))
  factor <- factor[samples, ]
  surv <- surv[which(surv[, 1] %in% samples), ]
  clin <- clin[which(clin[, 1] %in% samples), ]

  # define groups
  df <- .fct_cutpoint(factor = factor, surv)



}

#' @name .fct_cutpoint
#' @description Uses maxstat to find the optimal survival cutpoint or the median
#' for non-survival clinical features based on a given JDR factor.
#'
#' @inheritParams differential_analysis
#' @param factor Factor based on which to split the cohort. Expects a matrix.
#' @returns Data frame containing the information about the two survival groups.

.fct_cutpoint <- function(factor, surv, minprop = 0.1) {
  # cut the data
  time <- surv$time_to_event
  event <- surv$vital_status
  df <- data.frame(sample = samples, time = time, event = event,
                   factor = factor)
  cut <- survminer::surv_cutpoint(df, variables = "factor", minprop = minprop)
  df$FactorValue <- survminer::surv_categorize(cut)$Z
  df$FactorCluster <- df$Z > cut$cutpoint$cutpoint

  # determine which is the good survival group and which
  # is the poor survival group
  # not sure if this is the best way to do it
  df_split <- split(df, df$FactorCluster)
  if (max(df_split[[1]]$time) > max(df_split[[2]]$time)) {
    df_split[[1]]$group <- "long"
    df_split[[2]]$group <- "short"
  } else if (max(df_split[[1]]$time) < max(df_split[[2]]$time)) {
    df_split[[1]]$group <- "short"
    df_split[[2]]$group <- "long"
  }

  df2 <- rbind(df_split[[1]], df_split[[2]])
  df2 <- df2[order(df2$group), ]
  df2$group <- factor(df2$group)

  return(df2)
}