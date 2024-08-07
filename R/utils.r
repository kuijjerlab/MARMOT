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
    if (!any(names_exist)) {
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
    summarise_all(~ if (is.numeric(.)) var(.) == 0 else n_distinct(.) == 1) %>%
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
#' @param mat A matrix to test.
#' @noRd

.numeric_mat <- function(mat) {
  all(sapply(as.vector(mat), function(elem) {
    is.numeric(elem) || (is.character(elem) && !is.na(as.numeric(elem))) ||
      is.na(elem)
  }))
}

#' @name .pos_omics
#' @description Function that takes an omic matrix as input and transforms it so
#' that all values are positive and scaled between (0,1).
#' @param omic An omic matrix with rows as features and columns as samples.
#' @noRd

.pos_omics <- function(omic) {
  # ensure positivity
  if (min(omic) < 0) {
    omic_pos <- omic + abs(min(omic))
  } else {
    omic_pos <- omic
  }

  #scale between (0,1)
  omic_pos <- omic_pos / max(omic_pos)

  return(omic_pos)
}

#' @name .pad_matrix
#' @description Function to pad a matrix with NAs up to a specified number of
#' columns.
#' @param mat A matrix.
#' @param n_cols Number of columns.
#' @param mat_colnames A vector of column names for the padded matrix.
#' @noRd

.pad_matrix <- function(mat, n_cols, mat_colnames) {
  # If the matrix has fewer columns than max_cols, pad it with NA columns
  if (ncol(mat) < n_cols) {

    current_colnames <- colnames(mat)
    common_colnames <- intersect(mat_colnames, current_colnames)

    # Reorder columns of the current matrix to match the order in max_colnames
    mat <- mat[, common_colnames, drop = FALSE]

    # Create a matrix with NA columns to fill in
    na_cols <- matrix(NA, nrow = nrow(mat), ncol = n_cols - ncol(mat))

    # Combine the original matrix with the NA columns
    mat <- cbind(mat, na_cols)

    # Ensure the columns are in the same order as max_colnames
    colnames(mat) <- mat_colnames
    mat <- mat[, mat_colnames, drop = FALSE]
  }
  return(mat)
}

#' @name .pad_mat_wrapper
#' @description A wrapper that applies .pad_matrix to a list of matrices,
#' preserving column names and order.
#' @param mat_list A list of matrices
#' @noRd

.pad_mat_wrapper <- function(mat_list) {
  # Find the maximum number of columns
  max_cols <- max(sapply(mat_list, ncol))

  # Get all column names, preserving order
  max_colnames <- colnames(mat_list[[which.max(sapply(mat_list, ncol))]])

  # pad matrices
  mat_pad <- lapply(mat_list, function(mat) .pad_matrix(mat, max_cols, max_colnames)) #nolint

  return(mat_pad)
}

#' Convert Interval to Median Value
#'
#' This function converts an interval string to its median value.
#' The interval string can be in the form of ">=X" or "X-Y".
#'
#' @param interval A character string representing an interval.
#'        The interval can be in the format ">=X" or "X-Y".
#'
#' @return A numeric value representing the median of the interval.
#'         If the interval is ">=X", it returns X. If the interval is "X-Y", it returns the median of X and Y.
#'         If the input is `NA` or an unrecognized format, it returns `NA`.
#'
#' @examples
#' convert_interval_to_median(">=70") # Returns 70
#' convert_interval_to_median("50-69") # Returns 59.5
#' convert_interval_to_median(NA) # Returns NA
#' convert_interval_to_median("1-17") # Returns 9

.convert_interval_to_median <- function(interval) {
  if (is.na(interval)) {
    return(NA)
  } else if (grepl(">=", interval)) {
    # Handle the '>=X' case, assuming median to be X
    return(as.numeric(sub(">=", "", interval)))
  } else if (grepl("-", interval)) {
    # Handle the 'X-Y' case
    bounds <- as.numeric(unlist(strsplit(interval, "-")))
    return(median(bounds))
  } else {
    return(NA)
  }
}