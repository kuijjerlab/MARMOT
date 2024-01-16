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
#' @param PCA Logical. Whether PCA should be performed on the omics.
#' @param thresh NULL or numeric between 0-1. Threshold for the R2_cum to be used for selecting PCs; default 0.85.
#' Only needed if 'PCA' is TRUE. 
#' @param noPCs Numeric, indicating the number of principle components to keep. 
#' If 'thresh' is not NULL, it will ensure at least this many PCs are selected if fewer are needed to reach the threshold. 
#' @returns A list of omics ready for JDR. 
#' @examples

prepare_data <- function(omics, names = NULL, sep = NULL, overlap_samples = TRUE, PCA = TRUE, thresh = 0.85, noPCs = 20) {
  # create omics list
  omic <- .create_omics_list(omics, names, sep, overlap_samples)

  if (PCA){
   omic <- lapply(omic, .omic_pca(omic, thresh, noPCs))
  }
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

#' Perform PCA on a list of omics.
#' This function takes a list of omics and performs PCA on them. Expects the output of 
#' @inheritParams prepare_data
#' @param omic A matrix of omics data. 
#' @returns Returns a list of omics PCs that account for a specific variability threshold.

.omics_pca <- function(omic, thresh = 0.85, noPCs = NULL) {
  # data must be transposed for the PCA
  dat <- t(as.matrix(omic))
    
  # removing any features that are 0 in all samples
  dat <- dat[, which(colSums(dat[]) !=0)]
    
  # scaling data and running pca
  dat <- scale(dat)
  dat_pca <- pca(dat, nPcs = nrow(dat))
    
    # determine if to use R2_cum threshold or number of PCs 
    if(is.null(noPCs)){
      # make sure that each omic has at least 'noPCs' components 
      PCs <- which(dat_pca@R2cum <= thresh)
      if(length(PCs) < noPCs) {
        PCs <- 1:noPCs
      }
      omic_pca <- t(dat_pca@scores[, PCs]) 
      
    } else {
      PCs <- 1:noPCs
      omic_pca <- t(dat_pca@scores[, PCs])
    }
    
    pca_omics[[i]] <- dat_pca
  
  
  names(pca_omics) <- names(omic)
  names(MOFA_data) <- names(omic)
  
  can_pca <- list(pca_omics, MOFA_data)
  
  return(can_pca)
  
}
