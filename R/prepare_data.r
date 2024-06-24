################################################################################
#### Functions to manipulate omics  and clinical data and prepare it for JDR ###
################# and further analysis with or without PCA #####################
################################################################################

#' @title Prepare data for JDR.
#'
#' @name prepare_data
#'
#' @description This function prepares omics data for joint dimensionality
#' reduction.
#'
#' @param omics A string vector containing filenames for omics matrices;
#' accepts .RData files or any format accepted by fread; matrices must have
#' samples as columns and rows as features. The omics data must come from
#' the same samples, but some samples may have some missing values in some
#' omics.
#' @param names Character vector with the names to be used for the omics;
#' if left NULL, the omics will just be numbered.
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
#' @param scale_data Logical. Whether data should be scaled prior to performing
#' PCA.
#' @param logs Logical. Whether a log file containing information about the
#' input parameters should be saved.
#' @param log_name String indicating name for log file. If left NULL,
#' a generic name will be used.
#'
#' @returns A list of omics ready for JDR.
#' @examples
#' @export

prepare_data <- function(omics, names = NULL, overlap_samples = TRUE,
                         pca = TRUE, thresh = 0.85, n_pcs = 20, save_pca = TRUE,
                         file_name = NULL, scale_data = TRUE, logs = FALSE,
                         log_name = NULL) {
  # create log file
  if (logs) {
    if (is.null(log_name)){
      log_file <- file("prepare_data_log.txt", open = "a")
    } else {
      log_file <- file(log_name, open = "a")
    }

    timestamp <- format(Sys.time, "%d-%m-%Y %H:%M:%S")

    # Get the function's environment
    env <- environment()

    # Get the names and values of the function's arguments
    args <- as.list(env)[-1]
    arg_names <- names(args)
    arg_values <- lapply(args, function(arg) {
      if (is.character(arg) || is.numeric(arg) || is.logical(arg)) {
        as.character(arg)
      } else {
        paste("<", class(arg), "object>")
      }
    })

    # Write the timestamp and input parameters to the log file
    writeLines(c(paste("Timestamp:", timestamp),
                 "Input parameters:",
                 paste(paste(arg_names, arg_values, sep = ": "),
                       collapse = "\n"), ""), con = log_file)

    # Close the log file
    close(log_file)
  }

  # create omics list
  omic_list <- .create_omics_list(omics, names, overlap_samples)

  if (pca) {
    omics_pca <- lapply(omic_list, function(omic) .omics_pca(omic, scale_data))
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

#' @title Load omic data from file.
#'
#' @name .load_data
#'
#' @description Load data from file.
#'
#' @param file Character string indicating file name.
#'
#' @returns An omic matrix
#' @noRd

.load_data <- function(file) {
  # check and load file
  if (grepl(x = file, pattern = "\\.RData$|\\.Rda$", ignore.case = TRUE)) {
    return(as.matrix(get(load(file))))
  } else {
    # reading the data without the first column
    dat <- as.matrix(data.table::fread(file, drop = 1))

    # reading the first column separately and setting as row names
    rows <- data.table::fread(file, select = 1, colClasses = "character")[, 1]
    rownames(dat) <- rows

    return(dat)
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
#' @returns A list of omics matrices for JDR.
#' @noRd

.create_omics_list <- function(omics, names = NULL, overlap_samples = TRUE) {
  # read in omics and create a list
  omic <- lapply(omics, .load_data)

  # set omic names
  if (is.null(names)) {
    names <- paste0(rep("omic_", length(omic)), seq_along(length(omic)))
  }

  # test that the data is numeric or coercible to numeric
  is_num <- all(sapply(omic, .numeric_mat))
  if (!is_num) {
    stop("Data is not numeric. Please make sure all your omics are numeric.")
  }

  # set everything to lowercase
  # test if sample names are unique after lowercase

  # check if omics share at least one sample with each-other
  smpl_overlap <- length(Reduce(intersect, lapply(omic, colnames)))
  if (smpl_overlap == 0) {
    stop("No common samples between the omics. Please ensure omics share at least some samples.") #nolint
  }

  # warning if too few samples
  if (smpl_overlap < 15) {
    warning("Fewer than 15 samples common to all omics. This may cause unpredictable results.") # nolint
  }

  # get all sample names
  all_samples <- colnames(omic[[1]])

  # ensure only common samples are kept
  if (overlap_samples) {
    # Subset to only common samples
    omic <- .filter_omics(omic)
  } else {
    # pad matrices
    .pad_mat_wrapper(omic)
  }

  # name the omics
  if (length(names) == length(omic)) {
    names(omic) <- names
  } else {
    stop("Please make sure the length of 'names' is the same as the number of omics") # nolint
  }
  return(omic)
}

#' @title Filter omics to only overlapping samples.
#'
#' @name .filter_omics
#'
#' @description Function to filter omics to only overlapping samples.
#' This is required by some JDR methods.
#'
#' @param omic_list A list of omic matrices.
#'
#' @returns A list of omics that are filtered to only common samples.
#' @noRd

.filter_omics <- function(omic_list) {
  # Get all sample names
  all_samples <- lapply(omic_list, colnames)

  # Find common samples
  common_samples <- Reduce(intersect, all_samples)

  # Filter to only common samples
  omic_fil <- lapply(omic_list, function(x) x[, common_samples])

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
#' @returns PCA results for the omic matrix.
#' @noRd

.omics_pca <- function(omic, scale_data) {
  # data must be transposed for the PCA
  dat <- t(as.matrix(omic))

  # removing any features that are 0 in all samples
  # change to the .check_variance function
  dat <- dat[, which(colSums(dat[]) != 0)]

  # scaling data
  if (scale_data) {
    dat <- scale(dat)
  }
  # running pca

  dat_pca <- pcaMethods::pca(dat, nPcs = nrow(dat))

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
#'
#' @returns A matric of omic principle components filtered based on a
#' cummulative R2 threshold.
#' @noRd

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
#' @inheritParams prepare_data
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
#' @param vital_status_values A named list of length two indicating how vital
#' status is recorded in the clinical data. This will be converted to logical,
#' where alive = FALSE and dead = TRUE. Each element of the list can have
#' multiple entries. Default is list(alive = c(FALSE, "alive", "living", 0),
#' dead = c(TRUE, "dead", "deceased", 1))
#' @param keep_nas Logical. If samples with NA values in any survival feature
#' should be kept. Not recommended for most downstream analysis.
#'
#' @returns A data frame containing survival information.
#' @examples
#' @export

prepare_surv <- function(clinical, feature_names,
                         vital_status_values = list(alive = c(FALSE, "alive", "living", 0), # nolint
                                                    dead = c(TRUE, "dead", "deceased", 1)), # nolint
                         keep_nas = FALSE, logs = FALSE) {
  # load data
  clin <- data.table::fread(clinical)

  # sanity checks
  # check the features are labelled correctly
  .check_names(names(feature_names),
               c("sample_id", "vital_status", "time_to_event"),
               err_msg = "names(feature_names) are as detailed in the 
               documentation")

  # check if provided feature names exist in clinical data
  .check_names(feature_names, colnames(clin), err_msg = "feature names exist in 
  the clinical data")

  # grab specified features
  surv <- clin[, unlist(feature_names)]

  # merge columns that were provided for the same feature and rename columns
  result_list <- lapply(feature_names, .merge_surv, surv)
  surv_merged <- do.call(cbind, result_list)
  colnames(surv_merged) <- names(feature_names)

  #replace vital status with logical
  surv_merged$vital_status <- ifelse(!is.na(surv_merged$vital_status),
                                     toupper(surv_merged$vital_status) %in%
                                       toupper(vital_status_values[[2]]), NA)

  # make sure time is numeric and vital status is logical
  surv_merged$time_to_event <- as.numeric(surv_merged$time_to_event)
  surv_merged$vital_status <- as.logical(surv_merged$vital_status)

  #remove any remaining NAs
  if (!keep_nas) {
    surv_merged <- surv_merged[complete.cases(surv_merged), ]
  }

  # make sure it's a data frame and not data table
  surv_merged <- data.frame(surv_merged)

  return(surv_merged)
}

#' @title Merge multiple columns of surv data frame
#'
#' @name .merge_feat
#'
#' @description Merges the columns of a data frame if more than one feature
#' is provided for a single feature.
#'
#' @param cols Character vector containing column names.
#' @param surv Data frame of extracted survival features.
#'
#' @returns A data frame with merged columns.
#' @noRd

.merge_surv <- function(cols, surv) {
  # check if column vector has more than one element
  if (length(cols) > 1) {
    # get feature name
    col <- names(cols)
    # innitialise data frame
    temp <- data.frame()
    # if first feature is NA, replace with value from second
    temp <- ifelse(is.na(surv[, cols[1]]) & !is.na(surv[, cols[2]]),
                   surv[, cols[2]], surv[, cols[1]])
    colnames(temp) <- col
  } else {
    temp <- surv[, cols, drop = F]
  }

  return(temp)
}