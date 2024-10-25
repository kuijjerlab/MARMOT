################################################################################
####################### Functions to run JDR methods ###########################
################################################################################

#' @title Run joint dimensionality reduction on multi-omics data.
#'
#' @name run_jdr
#'
#' @description Function for performing joint dimensionality reduction (JDR)
#' on omics data using one or more linear JDR methods.
#' Note: If only one JRD method is used, it is prefered to use the
#' individual functions:
#' \itemize{
#'  \item{\code{\link{run_mofa2}}}
#'  \item{\code{\link{run_jive}}}
#'  \item{\code{\link{run_rgcca}}}
#'  \item{\code{\link{run_mcia}}}
#' }
#'
#' @param omic_list A list of omic matrices or path to .RData file containing
#' said list. Should be output of \code{link{prepare_data}}.
#' @param samples_overlap Logical. Whether a sample overlap was performed
#' between the omics during data preparation. See \code{\link{prepare_data}}.
#' If FALSE, omics will be filtered to only common samples for JDR methods that
#' require it. This may result in different sample sets being used with
#' different methods.
#' @param pca Logiacal. Whether PCA was performed on the data.
#' See \code{\link{prepare_data}}.
#' @param jdr_methods Character vector specifying one more more JDR methods
#' to be used. Should be from \code{c("MOFA", "JIVE", "RGCCA", "MCIA")}.
#' @param n_fct Integer. Number of factors for the factorisations.
#' @param seed Seed for factorisation.
#' @param convergence Only needed if MOFA2 is used. Determines the convergence
#' mode. Can be one of c("slow", "medium", "fast").
#' @param use_basilisk Logical. Whether basilisk should be used to automatically
#' set up conda environment. Only needed when MOFA2 is used.
#'
#' @return A list of factorisations. Each element is the factorisation
#' based on one JDR method.
#'
#' @examples
#'
#' @export
#'

run_jdr <- function(omic_list, samples_overlap = TRUE, pca = TRUE,
                    jdr_methods = c("MOFA", "JIVE", "RGCCA", "MCIA"),
                    n_fct = 5, seed = 42, convergence = "slow",
                    use_basilisk = TRUE) {

  # check if file path was provided instead of list and load file if so
  if (is.character(omic_list)) {
    omics <- get(load(omic_list))
  } else if (is.list(omics_list)) {
    omics <- omic_list
  } else {
    stop("omic_list must be either a list of matrices or a path to an .RData file.")
  }

  # overlap samples if not already done
  if (!samples_overlap) {
    omic_fil <- .filter_omics(omics)
  }

  # initialise factorisation list
  fct_list <- list()

  # run MOFA
  if (is.element("MOFA", jdr_methods)) {
    print("Performing JDR with MOFA2...")
    mofa_model <- run_mofa2(omic_list, n_fct, seed = seed,
                            convergence = convergence,
                            use_basilisk = use_basilisk)

    fct_list$MOFA <- MOFA2::get_factors(mofa_model)
    }

  # run JIVE
  if (is.element("JIVE", jdr_methods)) {
    print("Performing JDR with JIVE...")
    jive_model <- run_jive(omic_list, n_fct)

    # I don't actually know what this does
    # just copying cantini's code
    # will have to look into it more
    rankJV <- jive_model$rankJ
    J <- numeric(0)
    for(j in 1:length(omics)){
      J <- rbind(J, jive_model$joint[[j]])
    }
    svd.o <- svd(J)
    jV <- svd.o$v %*% diag(svd.o$d)
    colnames(jV) <- paste0("Factor", ncol(jV))
    rownames(jV) <- colnames(omic_list[[1]])

    fct_list$JIVE <- jV[, 1:rankJV]
  }

  # run RGCCA
  if (is.element("RGCCA", jdr_methods)) {
    print("Performing JDR with RGCCA...")
    if (!samples_overlap) {
      # run with omic_fil
      rgcca_model <- run_rgcca(omic_fil, n_fct)
      fct_list$RGCCA <- rgcca_model
    } else {
      # run with omic_list
      rgcca_model <- run_rgcca(omic_list, n_fct)
      fct_list$RGCCA <- as.matrix(rgcca_model$Y[[1]])
    }
  }

  # run MCIA
  if (is.element("MCIA", jdr_methods)) {
    print("Performing JDR with MCIA...")
    if (!samples_overlap) {
      # run with omic_fil
      mcia_model <- run_mcia(omic_fil, n_fct)
      fct_list$MCIA <- mcia_model
    } else {
      # run with omic_list
      mcia_model <- run_mcia(omic_list, n_fct)
      fct_list$MCIA <- as.matrix(mcia_model$mcoa$SynVar)
    }
  }

  return(fct_list)
}

#' @title Run JDR with MOFA2.
#'
#' @name run_mofa2
#'
#' @description Performs JDR on a list of omic matrices with MOFA2.
#'
#' @inheritParams run_jdr
#'
#' @return A trained MOFA2 model.
#'
#' @export
#' @importFrom MOFA2 run_mofa create_mofa prepare_mofa

run_mofa2 <- function(omic_list, n_fct, seed, convergence, use_basilisk) {
  # make MOFA object
  mofa_object <- create_mofa(omic_list)

  # prepare MOFA object
  mofa_object <- prepare_mofa(mofa_object)
  mofa_object@model_options$num_factors <- n_fct
  mofa_object@training_options$seed <- seed
  mofa_object@training_options$convergence_mode <- convergence

  # run mofa
  mofa_model <- run_mofa(mofa_object, use_basilisk = use_basilisk)

  return(mofa_model)
}

#' @title Run JDR with JIVE.
#'
#' @name run_jive
#'
#' @description Performs JRD on a list of omic matrices with JIVE.
#'
#' @inheritParams run_jdr
#'
#' @return A trained JIVE model.
#'
#' @export
#' @importFrom r.jive jive

run_jive <- function(omic_list, n_fct) {
  # run jive
  # look a little more into JIVE options and what they mean
  jive_model <- jive(omic_list, rankJ = n_fct,
                     rankA = rep(n_fct, length(omic_list)),
                     method = "given", conv = "default",
                     maxiter = 100)

  return(jive_model)
}

#' @title Run JDR with RGCCA.
#'
#' @name run_rgcca
#'
#' @description Performs JDR on a list of omic matrices using RGCCA.
#'
#' @inheritParams run_jdr
#'
#' @return A trained RGCCA model.
#'
#' @export
#' @importFrom RGCCA rgcca

run_rgcca <- function(omic_list, n_fct) {
  # transpose omics
  omics_t <- lapply(omic_list, function(x) t(x))

  # run rgcca
  rgcca_model <- rgcca(omics_t, ncomp = rep(n_fct, length(omics_t)))

  return(rgcca_model)
}

#' @title Run JDR with MCIA.
#'
#' @name run_mcia
#'
#' @description Performs JDR on a list of omic matrices using MCIA.
#'
#' @inheritParams run_jdr
#'
#' @return A trained MCIA model.
#'
#' @export
#' @importFrom omicade4 mcia

run_mcia <- function(omic_list, n_fct) {
  # make omics positive & scale between (0,1)
  # not sure why this is required specifically for this method
  # but that it what cantini did (but not for the other methods)
  # as far as I can tell, there is no specification anywhere in the docu
  # saying positive input is needed
  # so I dunno
  omics_pos <- lapply(omic_list, .pos_omics)

  # removing any features with uniform values across all samples
  # because MCIA has a stroke otherwise
  # actually, it only has a stroke if the uniform value is also the minimum value
  # i.e. 0
  for (i in seq_along(omic_list)) {

  }

  mcia_model <- mcia(omics_pos, cia.nf = n_fct)

  return(mcia_model)
}