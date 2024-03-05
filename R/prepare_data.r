#############################################################################
## Functions to manipulate omics  and clinical data and prepare it for JDR ##
################# and further analysis with or without PCA ##################
#############################################################################

#' @title Prepare data for JDR.
#'
#' @name prepare_data
#'
#' @description This function prepares omics data for joint dimensionality
#' reduction.
#'
#' @param omics A string vector containing filenames for omics matrices;
#' accepts text or .RData files; matrices must have samples as columns
#' and rows as features.The omics data must come from the same samples.
#' @param names Character vector with the names to be used for the omics;
#' if left NULL, the omics will just be numbered.
#' @param sep Character vector containing the separator used if text
#' files are provided.
#' @param overlap_samples Logical. Whether to ensure only samples with data
#' in all omics are kept. Some JDR methods require this.
#' @param pca Logical. Whether PCA should be performed on the omics.
#' @param thresh NULL or numeric between 0-1. Threshold for the R2_cum
#' to be used for selecting PCs; Only needed if 'PCA' is TRUE.
#' @param n_pcs Numeric, indicating the number of principle components to keep.
#' If 'thresh' is not NULL, it will ensure at least this many PCs are selected
#' if fewer are needed to reach the threshold.
#' @param save_pca Logical. Whether to save a copy of the unfiltered PCA results
#' as an RData file.
#' @param file_name Character string specifying the file name for the
#' unfiltered PCA results. If left NULL, a generic name will be assigned.
#' @param scale Logical. Whether data should be scaled prior to performing PCA.
#'
#' @return A list of omics ready for JDR.
#' @examples
#' @export

prepare_data <- function(omics, names = NULL, sep = NULL,
                         overlap_samples = TRUE, pca = TRUE,
                         thresh = 0.85, n_pcs = 20, save_pca = TRUE,
                         file_name = NULL, scale = TRUE) {
  # create omics list
  omic_list <- .create_omics_list(omics, names, sep, overlap_samples)

  if (pca) {
    omics_pca <- lapply(omic_list, function(omic) .omics_pca(omic, scale))
    if (save_pca) {
      if (is.null(file_name)) {
        file_name <- "omics_list_PCA.RData"
      }
      save(omics_pca, file = file_name)
    }

    # filtering PCA results
    omic_fil <- lapply(omics_pca, function(omic) .filter_pcs(omic, thresh, n_pcs)) #nolint

    return(omic_fil)
  } else {
    return(omic_list)
  }
}

#' @title Create a list of omic matrices for JDR.
#'
#' @name .create_omics_list
#'
#' @description Internal function that loads omics data and
#' binds them into a list.
#'
#' @inheritParams prepare_data
#' @return A list of omics matrices for JDR.

.create_omics_list <- function(omics, names, sep, overlap_samples) {
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
      stop("Invalid file format. Please make sure the data is in a text or .RData file format") #nolint
    }
  }

  # check if omics share at least one sample with each-other
  smpl_overlap <- length(Reduce(intersect, lapply(omic, colnames))) > 0
  if (smpl_overlap == FALSE) {
    stop("No common samples between the omics. Please ensure omics share at least some samples.") #nolint
  }

  # get all sample names
  all_samples <- lapply(omic, colnames)

  # ensure only common samples are kept
  if (overlap_samples) {
    # Subset to only common samples
    omic <- .filter_omics(omic)
  } else {
    # ensure samples are all in the same order
    omic <- lapply(omic, function(x) { x[, all_samples, drop = FALSE]})
  }

  # name the omics
  if (!is.null(names)) {
    if (length(names) == length(omic)) {
      names(omic) <- names
    } else {
      stop("Please make sure the length of 'names' is the same as the number of omics") #nolint
    }
  }
  return(omic)
}

#' @title Filter omics to only overlapping samples.
#'
#' @name .filter_omics
#'
#' @description Internal function to filter omics to only overlapping samples.
#' This is required by some JDR methods.
#'
#' @param omic_list A list of omic matrices.
#'
#' @return A list of omics that are filtered to only common samples.

.filter_omics <- function(omic_list) {
  # Get all sample names
  all_samples <- lapply(omic_list, colnames)

  # Find common samples
  common_samples <- Reduce(intersect, all_samples)

  # Filter to only common samples
  omic_fil <- lapply(omic_list, function(x) omic_list[, common_samples])

  return(omic_fil)
}

#' @title Perform PCA on a list of omics.
#'
#' @name .omics_pca
#'
#' @description Internal function that takes a list of omics and
#' performs PCA on them.
#'
#' @inheritParams prepare_data
#' @param omic A matrix of omics data.
#' @return PCA results for the omic matrix.
#' @importFrom pcaMethods pca

.omics_pca <- function(omic, scale) {
  # data must be transposed for the PCA
  dat <- t(as.matrix(omic))

  # removing any features that are 0 in all samples
  dat <- dat[, which(colSums(dat[]) != 0)]

  # scaling data
  if (scale) {
    dat <- scale(dat)
  }
  # running pca
  dat_pca <- pca(dat, nPcs = nrow(dat))

  return(dat_pca)
}

#' @title Filter PCA results based on a R2 threshold.
#'
#' @name .filter_pcs
#'
#' @description Internal function that filters a list of omics
#' based on a cummulative R2 threshold; alternatively,
#' it can select a specified number of principle components.
#'
#' @inheritParams prepare_data
#' @param omic_pca PCA results for one omic.
#' @return A matric of omic principle components filtered based on a
#' cummulative R2 threshold.

.filter_pcs <- function(omic_pca, thresh, n_pcs) {
  if (!is.null(thresh)) {
    # filter by R2 cum threshold
    PCs <- which(omic_pca@R2cum <= thresh)

    # make sure that each omic has at least 'n_pcs' principle components
    if (length(PCs) < n_pcs) {
      PCs <- 1:n_pcs
    }
    omic_fil <- t(omic_pca@scores[, PCs])

  } else {
    PCs <- 1:n_pcs
    omic_fil <- t(omic_pca@scores[, PCs])
  }

  return(omic_fil)
}

#' @title Extract survival information from clinical data.
#'
#' @name prepare_surv
#'
#' @description Function that extracts survival information from a larger
#' dataframe of clinical features and formats it for downstream analysis.
#'
#' @inbheritParams prepare_data
#' @param clinical Character string containing path to a file with clinical
#' information.The expected format of clinical data is with rows as samples
#' and columns as features.
#' @param feature_names Named list containing the feature names for
#' the survival features of interest. Two feature names can be
#' provided per element. If so, the first feature will be used unless a sample
#' has NA values for that feature, in which case it will check for values in the
#' second feature. See examples for details. The list must be named as follows:
#' list(sample_id = "sample_id", vital_status = "vital_status",
#' time_to_event = "time_to_event")
#' @param keep_nas Logical. If samples with NA values in any survival feature
#' should be kept.
#'
#' @return A data frame containing survival information.
#' @examples
#' @export

prepare_surv <- function(clinical, feature_names,
                         sep = "\t", keep_nas = FALSE) {
  # load data
  clin <- read.table(clinical, sep = sep, head = TRUE)

  # sanity checks
  # check the features are labelled correctly
  if (names(feature_names != c("sample_id", "vital_status", "time_to_event"))) {
    stop("List element names unknown. Please make sure names(feature_names) are
    as detailed in the documentation and are spelled correctly.")
  }

  # check if provided feature names exist in clinical data
  names_exist <- sapply(feature_names, function(x) all(x %in% colnames(clin)))
  stopifnot(names_exist, message = "Feature names could not be found. Please 
  make sure feature names exist in clinical data and are spelled correctly.")

  # grab specified features
  temp <- clin[, unlist(feature_names)]

  # merge columns that were provided for the same feature
  temp <- lapply()

  feat <- which(sapply(feature_names, function(x) length(x) > 1))
  if (length(feat > 0)) {
    temp <- lapply()
    temp <- ifelse(is.na(temp$days_to_death) & !is.na(temp$days_to_last_followup),
                                 temp$days_to_last_followup,
                                 temp$days_to_death)
  }


  # making sure the relevant columns are all named the same
  # not elegant, but you work with what you got
  if (TCGA) {
    clin <- clin %>%
      rename(
        sampleID = colnames(clin)[grep("submitter_id.samples", colnames(clin))],
        days_to_death = colnames(clin)[grep("days_to_death.diagnoses", colnames(clin))],
        vital_status = colnames(clin)[grep("vital_status.diagnoses", colnames(clin))],
        days_to_last_followup = colnames(clin)[grep("days_to_last_follow_up.diagnoses", colnames(clin))],
        days_to_last_followup = colnames(clin)[grep("days_to_lastfollowup", colnames(clin))]
      )

    #get the relevant columns
    temp <- data.frame(cbind(sampleID = clin$sampleID, 
                             vital_status = clin$vital_status, 
                             days_to_death = clin$days_to_death, 
                             days_to_last_followup = clin$days_to_last_followup, 
                             days_to_event = NA))

    #merge days_to_death & days_to_last_followup
    temp$days_to_death <- as.numeric(temp$days_to_death)
    temp$days_to_last_followup <- as.numeric(temp$days_to_last_followup)
    temp$days_to_event <- ifelse(is.na(temp$days_to_death) & !is.na(temp$days_to_last_followup),
                                 temp$days_to_last_followup,
                                 temp$days_to_death)

    #change sample names so they match with those in data
    temp$sampleID <- gsub("-", "\\.", temp$sampleID)

    #replace vital status with logical
    temp$vital_status <- gsub("DECEASED", TRUE, temp$vital_status, ignore.case = TRUE)
    temp$vital_status <- gsub("ALIVE", FALSE, temp$vital_status, ignore.case = TRUE)
    temp$vital_status <- gsub("dead", TRUE, temp$vital_status, ignore.case = TRUE)
    temp$vital_status <- gsub("LIVING", FALSE, temp$vital_status, ignore.case = TRUE)

    surv <- data.frame(sampleID = temp$sampleID,
                       vital_status = temp$vital_status,
                       days_to_event = temp$days_to_event)
  }else {
    surv <- data.frame(sampleID = clin$sampleID,
                       vital_status = clin$vital_status,
                       days_to_event = clin$time_to_event)
  }

  surv$days_to_event <- as.numeric(surv$days_to_event)
  surv$vital_status <- as.logical(surv$vital_status)

  #remove any remaining NAs
  surv <- surv[!is.na(surv$days_to_event), ]
  surv <- surv[!is.na(surv$vital_status), ]

  return(surv)
}

#' @title Merge multiple columns of surv data frame
#'
#' @name .merge_feat
#'
#' @description Merges the columns of a data frame if more than one feature
#' is provided for a single feature.
#'
#' @inheritParams prepare_surv
#' @param cols Column names
#' @param surv Data frame of extracted survival features.
#'
#' @return A data frame with merged columns.

.merge_surv <- function(cols, surv) {
  if (length(cols) > 1) {
    col <- names(cols)
    surv[, col] <- ifelse(is.na(surv[, cols[1]]) & !is.na(surv[, cols[2]]),
                          surv[, cols[2]], surv[, cols[1]])
  }
}