###

#' @name run_gsva
#'
#' @description Function to run GSVA on samples grouped by JDR factors.
#'
#' @inheritParam run_gsea

run_gsva <- function(omic, gene_set, method = "zscore") {
  # get gene set
  if (is.character(gene_set)) {
    gene_set <- fgsea::gmtPathways(gene_set)
  }

  # sanity check
  # check that omic and gene_set use the same annotation
  gs_names <- unlist(unique(gene_set))
  omic_names <- rownames(omic)
  .check_names(omic_names, gs_names, partial = TRUE,
               err_msg = "the gene sets and the omic data use the same annotation") # nolint

  gsva_res <- GSVA::gsva(omic, gene_set, method = method)

  return(gsva_res)

}