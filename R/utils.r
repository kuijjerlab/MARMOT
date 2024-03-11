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