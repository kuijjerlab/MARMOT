################################################################################
####### Functions to run GSEA based on differential analysis on based on #######
############################### JDR factors ####################################

#' @title Run GSEA based on differential analysis results
#'
#' @name perform_gsea
#'
#' @description This function performs Gene Set Enrichment Analysis (GSEA) on
#' results of differential analysis based on JDR factors. See
#' \code{\link{differential_analysis}} for details.
#'
#' @inheritParams differential_analysis
#' @param diff_results A data frame containing either a limma toptable or the
#' output of \code{link{differential_analysis}}.
#' @param gene_set A .gmt file containing gene sets of interest or a list of
#' gene sets to test.
#' @param ... Any other parameters accepted by fgsea.
#' @seealso \code{\link{fgsea::fgsea}}
#'
#' @returns Results of GSEA
#' @export

perform_gsea <- function(diff_results, limma = TRUE, gene_set, save_file = TRUE,
                         file_name = NULL, ...) {
  # get gene set
  if (is.character(gene_set)) {
    gene_set <- fgsea::gmtPathways(gene_set)
  }

  # sanity checks
  # check that omic and gene_set use the same annotation
  gs_names <- unlist(unique(gene_set))
  omic_names <- rownames(diff_results)
  .check_names(omic_names, gs_names, partial = TRUE,
               err_msg = "the gene sets and the omic data use the same annotation") # nolint

  # create rank
  rnk <- .create_rank(diff_results, limma)

  # get gene set
  if (is.character(gene_set)) {
    gene_set <- fgsea::gmtPathways(gene_set)
  }

  gsea_res <- fgsea::fgsea(pathways = gene_set, rnk, ...)


  if (save_file) {
    if (is.null(file_name)) {
      file_name <- "gsea_analysis_results.RData"
    }
    save(rnk, gsea_res, file = file_name)
  }
  return(gsea_res)
}

#' @name select_stable_path
#' @description Function that takes multiple gsea results and selects
#' pathways that are common to all of them.
#' @param gsea_res A list of dataframe outputs of \code{\link{fgsea::fgsea}}
#'
#' @returns A list with the subsetted data frames.

select_stable_path <- function(gsea_res) {
  # find common pathways
  common_path <- purrr::reduce(gsea_res, ~ dplyr::inner_join(.x, .y))$pathway

  # subset to only common pathways
  gsea_subset <- purrr::map(gsea_res, ~ filter(.x, pathway %in% common_path))

  return(gsea_subset)
}

#' @name .create_rank
#' @description Function to create ranks either from limma or wilcox tests for
#' GSEA. For limma, it simply extracts the test statistic. For wilcox, it
#' multiplies the -log10 p-value with the log2 FC to get directionality.
#'
#' @inheritParams perform_gsea
#'
#' @returns A dataframe of gene ranks.
#' @noRd

.create_rank <- function(diff_results, limma) {
  if (limma) {
    # order by absolute value of statistic
    diff_results <- diff_results[order(abs(diff_results$t),
                                       decreasing = TRUE), ]
    rnk <- diff_results$t
  } else {
    diff_results <- diff_results[order(abs(diff_results$logFC),
                                       decreasing = TRUE), ]
    rnk <- diff_results$W * diff_results$logFC
  }

  # name rnk
  names(rnk) <- rownames(diff_results)

  return(rnk)
}

#' @name .get_leading_edge
#'
#' @description Function to extract the leading edge from GSEA results and output it in a
#' readable format
#'
#' @param gsea_res A list of GSEA results.

.get_gsea <- function(gsea_res) {
  

  genes <- unlist(strsplit(str, "/"))
  genes <- gsub("/", "", genes)  # Remove "/"
  return(genes)
}