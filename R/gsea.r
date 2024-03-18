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
  # create rank
  rnk <- .create_rank(diff_results, limma)

  # get gene set
  gmt <- grepl(".gmt", gene_set, ignore.case = TRUE)
  if (gmt) {
    gset <- fgsea::gmtPathways(gene_set)
    gsea_res <- fgsea::fgsea(pathways = gset, rnk, ...)
  } else {
    gsea_res <- fgsea::fgsea(pathways = gene_set, rnk, ...)
  }

  if (save_file) {
    if (is.null(file_name)) {
      file_name <- "gsea_analysis_results.RData"
    }
    save(list = c(rnk, gsea_res), file = file_name)
  }
  return(gsea_res)
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
    rnk <- diff_results$logp * diff_results$logFC
  }

  return(rnk)
}