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
#' @inheritParams clin_association
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
                                  minprop = 0.1, sample_label = NULL) {
  # sanity checks
  if (!is.null(covariates)) {
    # check that clinical data is provided if covariates are provided
    if (is.null(clin)) {
      stop("Please provide clinical data for covariate correction.")
    }
    # check that covariates exist
    .check_names(covariates, colnames(clin),
                 err_msg = "covariate names exist in the clinical data")

    # check that sample label exists
    .check_names(sample_label, colnames(clin),
                 err_msg = "sample_label exists in the clinical data")

    # check that limma is selected as the method
    if (!limma) {
      stop("Limma must be used for covariate correction.")
    }
  }

  # check that sample names for survival and clinical are the same
  .check_names(rownames(clin), colnames(surv), partial = TRUE
               err_msg = "sample names are the same in clinical and survival data") # nolint

  # overlap samples
  # Define the sets
  sets <- list(
    surv_samples = surv$sample_id,
    factor_samples = rownames(factor),
    clin_samples = clin[, sample_label],
    omic_samples = colnames(omic)
  )

  # Compute the intersection of sets
  samples <- Reduce(intersect, sets)

  # subset data frames
  factor <- factor[samples, ]
  surv <- surv[which(surv$sample_id %in% samples), ]
  clin <- clin[which(clin[, sample_label] %in% samples), ]
  omic <- omic[, samples]

  # define groups
  df <- .fct_cutpoint(factor = factor, surv)

  # grab covariates
  cov <- clin[, covariates]

  cov <- cov[order(match(cov[, sample_label], df$sample)),]

  # add to df
  df <- cbind(df, cov)

  # remove any rows that contain NAs
  df <- na.omit(df)

  # create design matrix
  if (limma) {
    toptable <- .run_limma(omic, df, covariates)
  }

  
}

#' @name .fct_cutpoint
#' @description Uses maxstat to find the optimal survival cutpoint or the median
#' for non-survival clinical features based on a given JDR factor.
#'
#' @inheritParams differential_analysis
#' @param factor Factor based on which to split the cohort. Expects a matrix.
#' @returns Data frame containing the information about the two survival groups.

.fct_cutpoint <- function(factor, surv, minprop = 0.1) {
  # get samples
  samples <- surv[, 1]

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

#' @title Run limma on omic.
#'
#' @name .run_limma
#'
#' @description Run limma on an omic using groups defined by a JRD factor.
#'
#' @inheritParams differential_analysis
#' @param df Data frame. Output of .fct_cutpoint
#'
#' @returns A limma toptable.
#' @noRd

.run_limma <- function(omic, df, covariates) {
    # create formula with covariates
    formula <- as.formula(paste("~ 0 + group +", paste(covariates,
                                                       collapse = " + ")))
    design <- limma::model.matrix(formula, data = df)
    
    # specify contrasts
    contrasts <- limma::makeContrasts(short_vs_long = groupshort - grouplong,
                                levels = design)
    fit <- limma::eBayes(lmFit(omic, design))
    fit2 <- limma::contrasts.fit(fit, contrasts)
    fit2 <- limma::eBayes(fit2)
    toptable <- limma::topTable(fit2, coef="short_vs_long", number=Inf)
    toptable <- toptable[order(row.names(toptable)), ]

    return(toptable)
}