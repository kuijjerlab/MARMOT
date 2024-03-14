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
#' @param survival Logical. Whether the feature of interest is survival. This
#' will determine how the groups are divided.
#'
#' @export

differential_analysis <- function(omic, factor, clin, feat, covariates = NULL,
                                  limma = TRUE, survival = TRUE) {
  # sanity checks
  # check that feature name exists in clin data frame
  .check_names(feat, colnames(clin), err_msg = "feature name exists in the clinical data") # nolint

  # check that covariates exist
  if (!is.null(covariates)) {
    .check_names(covariates, colnames(clin). err_msg = "covariate names exist in the clinical data") # nolint
  }
}

#' @name .fct_cutpoint
#' @description Uses maxstat to find the optimal survival cutpoint or the median
#' for non-survival clinical features based on a given JDR factor.
#' @param factor Factor based on which to split the cohort. Expects a matrix.
#' @param surv Data frame containing the survival data
#' @param minprop Numeric between c(0,1), indicating the minimum proportion of samples per group.
#' @returns Data frame containing the information about the two survival groups.

.fct_cutpoint <- function(factor, surv, minprop = 0.1) {

  # only take samples for which there is survival info
  samples <- surv$sample_id
  samples <- intersect(samples, rownames(factor))
  factor <- factor[samples, ]
  surv <- surv[which(surv$sample_id %in% samples), ]

  # cut the data
  time <- surv$time_to_event
  event <- surv$vital_status
  df <- data.frame(sample = samples, time = time, event = event, Z = factor)
  cut <- surv_cutpoint(df, variables = "Z", minprop = minprop)
  df$FactorValue <- surv_categorize(cut)$Z
  df$FactorCluster <- df$Z > cut$cutpoint$cutpoint
  
  # determine which is the good survival group and which is the poor survival group
  df.split <- split(df, df$FactorCluster)
  if(max(df.split[[1]]$time) > max(df.split[[2]]$time)){
    df.split[[1]]$group <- "long"
    df.split[[2]]$group <- "short"
  }else if(max(df.split[[1]]$time) < max(df.split[[2]]$time)){
    df.split[[1]]$group <- "short"
    df.split[[2]]$group <- "long"
  }
  
  # determine which group has the high factor values and which has low factor values
  # if(max(df.split[[1]]$Z) < min(df.split[[2]]$Z)){
  #   df.split[[1]]$FactorValue <- "low"
  #   df.split[[2]]$FactorValue <- "high"
  # }else if(min(df.split[[1]]$Z) > max(df.split[[2]]$Z)){
  #   df.split[[1]]$FactorValue <- "high"
  #   df.split[[2]]$FactorValue <- "low"
  # }
  
  df2 <- rbind(df.split[[1]], df.split[[2]])
  df2 <- df2[order(df2$group),]
  df2$group <- factor(df2$group, levels = c("long", "short"))

  return(df2)
}