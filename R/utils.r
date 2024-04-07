################################################################################
###################### Various small utility functions #########################
################################################################################

#' @title Check if all elements of a character vector exist as elements in
#' another character vector.
#'
#' @name .check_colnames
#'
#' @description Check if all elements of a character vector exist as elements in
#' another character vector.
#'
#' @param vector_1 Vector to be checked.
#' @param vector_2 Reference vector.
#' @param err_msg Optional. A custom error message so I can give some meaningful
#' errors.
#' @param partial Logical, if to look for a full match or only partial.
#'
#' @noRd

.check_names <- function(vector_1, vector_2, err_msg = NULL, partial = FALSE) {
  names_exist <- sapply(vector_1, function(x) all(x %in% vector_2))
  if (!partial) {
    if (!all(names_exist)) {
      stop(paste("Elements could not be found. Please make sure",
                 err_msg, "and are spelled correctly.", sep = " "))
    }
  } else {
    if (any(!names_exist)) {
      stop(paste("Elements could not be found. Please make sure",
                 err_msg, "and are spelled correctly.", sep = " "))
    }
  }
}

#' @title Check if each column in a df has some variance.
#'
#' @name .check_variance
#'
#' @description Check if each column in a df has some variance and remove
#' columns that don't.
#'
#' @param df Data frame to check.
#'
#' @returns A filtered data frame.
#'
#' @importFrom dplyr summarise_all filter everything
#' @importFrom tidyr pivot_longer
#' @importFrom magrittr %>%
#'
#' @noRd

.check_variance <- function(df) {
  constant_columns <- df %>%
    summarise_all(~ if(is.numeric(.)) var(.) == 0 else n_distinct(.) == 1) %>%
    pivot_longer(everything(), names_to = "Column", values_to = "Constant") %>%
    filter(Constant)

  if (nrow(constant_columns) > 0) {
    warning("The folowing features have no varaince and have been removed:",
            paste(constant_columns$Column, collapse = ", "))
  }

  result <- df[, setdiff(colnames(df), constant_columns$Column)]

  return(result)
}

#' @title reduces a set of vectors to common elements
#'
#' @name .overlap_sets
#'
#' @inheritParams .check_names # not sure yet if to implement err_msg here
#' @param sets A list of vectors to overlap.
#'
#' @returns A vector with all the common elements.
#' @noRd

.overlap_sets <- function(sets, err_msg = NULL) {
  # check that at least some elements are common to all vectors
  .check_overlaps(sets, err_msg)

  # overlap vectors
  result <- Reduce(intersect, sets)
  return(result)
}

#' @name .check_overlaps
#'
#' @description Function to check that a list of vectors have at least some
#' overlapping elements.
#'
#' @inheritParams .overlap_sets
#' @inheritParams .check_names
#'
#' @noRd

.check_overlaps <- function(sets, err_msg = NULL) {
  overlap <- Reduce(intersect, sets)

  if (length(overlap) == 0) {
    # have to come up with a better error message at some point
    # and maybe implement the custom err_msg
    stop(err_msg)
  }
}

# For now, this can live here.
#' @name .normalise_indeg
#' @description This function takes in a list of omic matrices that includes
#' indegree and quantile normalises the indegree
#' @param data A list of omic matrices (including indegree)
#' @param indeg_label Character stings indicating the label of the list item
#' containing the indegree. Default is "indegree".
#' @returns A list of omic matrices with normalised indegree.
#' @noRd

.normalise_indeg <- function(data, indeg_label = "indegree") {
  indeg <- data[[indeg_label]]
  indeg_n <- preprocessCore::normalize.quantiles(indeg, copy = FALSE)

  data_n <- data
  data_n[[indeg_label]] <- indeg_n

  return(data_n)
}

#' @name .numeric_mat
#' @description Function to check if a matrix is numeric.
#' @param A matrix to test.
#' @noRd
#'
.numeric_mat <- function(mat) {
  all(sapply(as.vector(mat), function(elem) {
    is.numeric(elem) || (is.character(elem) && !is.na(as.numeric(elem))) ||
      is.na(elem)
  }))
}
