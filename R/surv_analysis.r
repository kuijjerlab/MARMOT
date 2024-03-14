################################################################################
#### Functions for performing factor association with survival and various #####
########################## other associated analyses ###########################

#' @title Calculate factor association with survival.
#'
#' @name surv_association
#'
#' @description This function takes as input computed factorisations from JDR
#' models (see \code{\link{run_jdr}}) and calculates the association to survival
#' of the factors. There are options for either univariate cox regression for
#' each factor, or multivariate cox regression for multiple factors.
#'
#' @param factors Data frame containing factors to be associated with survival.
#' Samples should be rows and factors columns.
#' @param surv A data frame with survival information. Expects output from
#' \code{\link{prepare_surv}}.
#' @param univariate Logical, whether to run a univariate coxph for each
#' factor (TRUE) or a multivariate coxph for all together (FALSE).
#' Default is TRUE.
#'
#' @returns A list with either one multivariate coxph model or multiple
#' univariate coxph models.
#'
#' @export

surv_association <- function(factors, surv, univariate = TRUE) {

  # make sure samples are the same
  samples <- surv$sample_id
  samples <- intersect(samples, rownames(factors))

  # make sure factors are named
  if (is.null(colnames(factors))) {
    colnames(factors) <- paste0("Factor", seq(1, ncol(factors)))
  }

  # subset to only samples that have survival data
  factors <- factors[samples, ]
  surv <- surv[which(surv$sampleID %in% samples), ]

  # create survival object
  survival_object <- survival::Surv(surv$time_to_event, surv$vital_status)

  # run either one multivariate or multiple univariate coxph models
  cox <- list()

  if (univariate) {
    cox <- lapply(colnames(factors), function(colname) {
      temp <- survival::coxph(survival_object ~ factors[, colname])
      names(temp$coefficients) <- colname
      return(temp)
    })
    names(cox) <- colnames(factors)
  } else {
    cox[[1]] <- survival::coxph(survival_object ~ factors)
  }

  return(cox)
}

#' @name surv_compare
#'
#' @description This function takes one or more coxph outputs and creates a
#' data frame that summarises the results and makes them easy to compare across
#' multiple models.
#'
#' @inheritParams surv_association
#' @param model_labels Laabels for the coxph models input. This should be in the
#' same order as the models passed to the "models" variable. Will be used to
#' distinguish between models for comparison.
#' @param models Trained JDR models. Must be one of the supported JDR methods.
#' See \code{link{run_jdr}}.
#' @param p_adjust Whether the pvalues should be corrected for multiple testing.
#' @param method Method for pvalue correction. Can be c("holm", "hochberg",
#' "hommel", "bonferroni", "BH", "BY", "fdr", "none").
#' @param logtrans Logical. Whether to include a -log10() transformation of the
#' p-values. Recommended for plotting. If p_adjust = TRUE, the adjusted p-values
#' will be transformed, otherwise, the raw p-values will be used.
#'
#' @returns A dataframe with a summary of the information obtained from
#' \code{\link{calculate_surv_association}} for one or multiple JDR models.
#'
#' @export

surv_compare <- function(models, model_labels, univariate = TRUE,
                         p_adjust = TRUE, method = "BH", logtrans = TRUE) {
  #initialise surv data frame
  surv_df <- data.frame(factor = character(),
                        hazard_ratio = numeric(),
                        pval = numeric(),
                        upper = numeric(),
                        lower = numeric(),
                        label = character())

  # innitialise pvals
  pvals2 <- c()

  for (i in seq_along(models)) {

    cox <- models[[i]]

    if (univariate) {
      # Extract information for each factor in the model
      model_data <- lapply(cox, function(model) {
        s <- summary(model)
        cfs <- coef(s)
        fct <- rownames(cfs)
        hr <- cfs[, "exp(coef)"]
        pvals <- cfs[, "Pr(>|z|)"]
        u <- s[["conf.int"]][, "upper .95"]
        l <- s[["conf.int"]][, "lower .95"]

        temp <- data.frame(
          factor = fct,
          hazard_ratio = as.numeric(hr),
          pval = as.numeric(pvals),
          upper = as.numeric(u),
          lower = as.numeric(l),
          label = model_labels[i]
        )
      })

      df <- dplyr::bind_rows(model_data)

      if (p_adjust) {
        df$padj <- p.adjust(df$pval)
      }

      surv_df <- rbind(surv_df, df)
    } else {
      s <- summary(cox)
      cfs <- coef(s)
      # get factor names and rename for convenience
      fct <- rownames(cfs)
      fct <- gsub("[, 1]", "", rownames(fct))

      # get hazard ratio
      hr <- cfs[, "exp(coef)"]

      # get pval
      pvals <- cfs[, "Pr(>|z|)"]

      # get confidence intervals
      u <- s[["conf.int"]][, "upper .95"]
      l <- s[["conf.int"]][, "lower .95"]

      #get model label
      lbl <- rep(model_labels[i], length(fct))

      # get pvalues for all factors for later adjustment
      if (p_adjust) {
        padj <- p.adjust(pvals)
        temp <- cbind(factor = fct,
                      hazard_ratio = as.numeric(hr),
                      pval = as.numeric(pvals),
                      padj = padj,
                      upper = as.numeric(u),
                      lower = as.numeric(l),
                      label = lbl)
      } else {
        # combine them in a df
        temp <- cbind(factor = fct, 
                      hazard_ratio = as.numeric(hr),
                      pval = as.numeric(pvals),
                      upper = as.numeric(u),
                      lower = as.numeric(l),
                      label = lbl)
      }

      surv_df <- rbind(surv_df, temp)
    }
  }

  # log transform p values
  if (logtrans) {
    if (p_adjust) {
      surv_df$logp <- -log10(surv_df$padj)
    } else {
      surv_df$logp <- -log10(surv_df$pval)
    }
  }

  return(surv_df)
}