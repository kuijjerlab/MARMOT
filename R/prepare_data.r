#' Prepare data for JDR.
#'
#' This function prepares omics data for joint dimensionality reduction. 
#' @param omics A string vector containing filenames for omics matrices; accepts text or .RData files;
#'               matrices must have samples as columns and rows as features.
#'               The omics data must come from the same samples.
#' @param names Character vector with the names to be used for the omics;
#'               if left NULL, the omics will just be numbered.
#' @param sep Character vector containing the separator used if text files are provided.
#' @param overlap_samples Logical. Whether to ensure only samples with data in all omics are kept.
#'                      Some JDR methods require this.

prepare_data <- function(omics, names = NULL, sep = NULL, overlap_samples = TRUE){

}

#' Create a list of omic matrices for JDR.
#'
#' This function loads omics data and binds them into a list.
#'
#' @inheritParams prepare_data
#' @return A list of omics matrices for JDR.
#' @examples 

.create_omics_list <- function(omics, names = NULL, sep = NULL, overlap_samples = TRUE) {
  #initialise list of omics
  omic <- list()

  #sanity checks and loading input files
  for (i in seq_along(omics)) {
    if (grepl(x = omics[i], pattern = ".RData" | ".Rda", ignore.case = TRUE)) {
      omic[[i]] <- as.matrix(get(load(omics[i])))
    } else if (grepl(x = omics[i], pattern = ".txt" | ".tsv" | ".csv")) {
      if (is.null(sep)) {
        stop("Please provide a separator for your text files.")
      } else if (length(sep) > 1) {
        stop("Please only provide one separator for your text files.")
      } else {
        omic[[i]] <- as.matrix(read.table(omics[i], head = TRUE, sep = sep))
      }
    } else {
      stop("Invalid file format. Please make sure the data is in a text or .RData file format")
    }
  }

  # check if omics share at least one sample with each-other
  smpl_overlap <- length(Reduce(intersect, lapply(omic, colnames))) > 0
  if (smpl_overlap == FALSE) {
    stop("No common samples between the omics. Please ensure omics share at least some samples.")
  }

  # get all sample names
  all_samples <- lapply(omic, colnames)

  # ensure only common samples are kept
  if (overlap_samples) {
    # Find common samples
    common_samples <- Reduce(intersect, all_samples)

    # Subset to only common samples
    omic <- lapply(omic, function(x) x[, common_samples])
  } else {
    # ensure samples are all in the same order
    omic <- lapply(omic, function(x) { x[, all_samples, drop = FALSE]})
  }

  # name the omics
  if (!is.null(names)) {
    if (length(names) == length(omic)) {
      names(omic) <- names
    } else {
      stop("Please make sure the length of 'names' is the same as the number of omics")
    }
  }
  return(omic)
}
