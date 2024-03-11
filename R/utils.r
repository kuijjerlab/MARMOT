################################################################################
###################### Various small utility functions #########################
################################################################################

#' @title Check if elements of a character vector exist as column names in a
#' data frame.
#'
#' @name .check_colnames
#'
#' @description Check if elements of a character vector exist as column names
#' in a data frame.
#'
#' @param df Dataframe for which column names should be checked.
#' @param name_vector A character vector containing column names to be checked.
#'

.check_colnames <- function(df, name_vector) {
  names_exist <- sapply(name_vector, function(x) all(x %in% colnames(df)))
  if (!all(names_exist)) {
    stop("Feature names could not be found. Please make sure the column names 
    provided exist in your data frame and are spelled correctly")
  }
}