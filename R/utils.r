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
#' error messages depending on the context in which I use this.

.check_names <- function(vector_1, vector_2, err_msg = NULL) {
  names_exist <- sapply(vector_1, function(x) all(x %in% vector_2))
  if (!all(names_exist)) {
    stop(paste("Elements could not be found. Please make sure",
               err_msg, "and are spelled correctly.", sep = " "))
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
#' @return A filtered data frame.
#'
#' @importFrom dplyr summarise_all pivot_longer filter everything
#' @importFrom magrittr %>%

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

# For now, this can live here. 
#' @name .normalise_indeg
#' @description This function takes in a list of omic matrices that includes
#' indegree and quantile normalises the indegree
#' @param data A list of omic matrices (including indegree)
#' @param indeg_label Character stings indicating the label of the list item
#' containing the indegree. Default is "indegree".
#' @return A list of omic matrices with normalised indegree.

.normalise_indeg <- function(data, indeg_label = "indegree") {
  indeg <- data[[indeg_label]]
  indeg_n <- normalize.quantiles(indeg, copy = FALSE)

  data_n <- data
  data_n[[indeg_label]] <- indeg_n

  return(data_n)
}