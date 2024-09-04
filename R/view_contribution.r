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
#' @param labels Option character vector with labels for the two sets of factors.
#' @param as_data_frame Logical, indicating whether the output should be a
#' data frame instead of a list of matrices.
#' @param abs Logical. Whether absolute correlation should be reported.
#'
#' @returns A list containing the factor correlation matrices.
#'
#' @export

fct_corr <- function(set1, set2, labels = NULL, as_data_frame = TRUE,
                     abs = TRUE) {
  # sanity checks
  # check that the number of features is the same in the two sets
  if (nrow(set1) != nrow(set2)) {
    warning("Differing numbers of samples (rows) between sets of factors. Results may be unreliable.") # nolint
  }

  # add labels
  if (!is.null(labels)) {
    colnames(set1) <- paste(colnames(set1), labels[1], sep = " ")
    colnames(set2) <- paste(colnames(set2), labels[2], sep = " ")
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

  # concatenate the two matrices into a list
  corr_res <- list(corr_pears, corr_spear)

  # take absolute value
  if (abs) {
    corr_res <- lapply(corr_res, abs)
  }

  # ensure outputs are matrices
  if (!is.matrix(corr_pears)) {
    corr_res <- lapply(corr_res, as.matrix)
    # transpose matrices in the edge case of one of the sets
    # only having one column
    if (ncol(corr_res[[1]]) == 1) {
      corr_res <- lapply(corr_res, t)
    }
  }

  # ensure labels are correct
  corr_res <- lapply(corr_res, function(x) {
    rownames(x) <- colnames(set2)
    return(x)
  })

  corr_res <- lapply(corr_res, function(x) {
    colnames(x) <- colnames(set1)
    return(x)
  })

  names(corr_res) <- c("pearson", "spearman")

  return(corr_res)
}

#' @name format_fct_corr
#'
#' @description Function to format the factor correlations results for plotting.
#'
#' @inheritParams fct_corr
#' @param corr_res A list of correlation matrices. Expects output of
#' \code{\link{fct_corr}}.
#'
#' @returns A long data frame ready for plotting.
#'
#' @export
#'
format_fct_corr <- function(corr_res) {
  # sanity checks
  # check that input is a list
  if (!is.list(corr_res)) {
    stop("Invalid input format. This function expects a list of two matrices as input.") # nolint
  }

  # format matrices as data frames
  corr1 <- reshape2::melt(corr_res[[1]])
  corr1$method <- names(corr_res)[1]
  corr2 <- reshape2::melt(corr_res[[2]])
  corr2$method <- names(corr_res)[2]

  # concatenate the two methods
  corr <- rbind(corr1, corr2)

  return(corr)
}