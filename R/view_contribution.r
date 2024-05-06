################################################################################
## Functions relating to testing the contribution of each view to a MOFA model #
################################################################################

#' @name fct_corr
#'
#' @description Perform pearson and spearman correlation between two sets of JDR
#' factors. Primarily intended for comparing different JDR models on the same
#' data.
#'
#' @inheritParams view_exclusion
#' @param set1 A matrix with the first set of factors.
#' @param set2 A matrix with the second set of factors.
#' @param labels Option character vector with labels for the two sets of factors.
#' @param as_data_frame Logical, indicating whether the output should be a
#' data frame instead of a list of matrices.
#'
#' @returns A list containing the factor correlation matrices.
#'
#' @export

fct_corr <- function(set1, set2, labels = NULL, as_data_frame = TRUE) {
  # sanity checks
  # check that the number of features is the same in the two sets
  if (nrow(set1) != nrow(set2)) {
    warning("Differing numbers of samples (rows) between sets of factors. Results may be unreliable.") # nolint
  }

  # add labels
  if (!is.null(labels)) {
    colnames(set1) <- paste(colnames(set1), labels[1], sep = "_")
    colnames(set2) <- paste(colnames(set2), labels[2], sep = "_")
  }

  # perform pairwise correlations of each factor in one set with the others
  corr_pears <- apply(set1, 2, function(x) {
    apply(set2, 2, function(y) {
      cor(x, y, method = "pearson")
    })
  })

  corr_spear <- apply(set1, 2, function(x) {
    apply(set2, 2, function(y) {
      cor(x, y, method = "spearman")
    })
  })

  if (as_data_frame) {
    corr_pears <- reshape2::melt(corr_pears)
    corr_spear <- reshape2::melt(corr_spear)
    corr_pears$method <- "pearson"
    corr_spear$method <- "spearman"
    corr_res <- rbind(corr_pears, corr_spear)
  } else {
    corr_res <- list(corr_pears, corr_spear)
    names(corr_res) <- c("pearson", "spearman")
  }

  return(corr_res)
}

#' @name plot_fct_corr
#'
#' @description Plots heatmaps of factor correlations.
#'
#' @inheritParams plot_data_dim
#' @param corr_res A list or data frame with the correlation values.
#' Expects output of \code{\link{fct_corr}}.
#' @param grid Logical. If true one  grid plt will be output with a panel for
#' spearman and one for pearson. If FALSE, a list will be output with each
#' element being one of the plots.
#' @param title Optional. Character string indicating title for the grid.
#' Only used if grid = TRUE.
#' @param ... Any other ggplot parameters.
#'
#' @returns A list of two ggplots. One for pearson and one for spearman.
#'
#' @export
#' @import ggplot2

plot_fct_corr <- function(corr_res, grid = TRUE, colours = NULL, title = NULL, ...) {
  # set colours
  if (is.null(colours)) {
    col <- RColorBrewer::brewer.pal(name = "Dark2", n = 8)
    col <- col[c(3, 4)]
  } else {
    if (length(colours) != 1) {
      stop(paste0(length(colours), " colours were specified, when 2 were expected. ",
      "Please make sure you specify the correct number of colours."))
    }
    col <- colours
  }

  if (class(corr_res) == list) {
    # Convert matrix to a data frame
    pears <- reshape2::melt(corr_res[["pearson"]])
    spear <- reshape2::melt(corr_res[["spearman"]])
  }


  # Plot heatmap
  p <- ggplot(pears, aes(x = Var2, y = Var1, fill = value, label = round(value, 2))) +
    geom_tile() +
    geom_text(color = "black") +
    scale_fill_gradient2(low = col[1], mid = "white", high = col[2], midpoint = 0) +
    labs(title = "Pearson", x = NULL, y = NULL, fill = "r") +
    theme_bw()

  s <- ggplot(spear, aes(x = Var2, y = Var1, fill = value, label = round(value, 2))) +
    geom_tile() +
    geom_text(color = "black") +
    scale_fill_gradient2(low = col[1], mid = "white", high = col[2], midpoint = 0) +
    labs(title = "Spearman", x = NULL, y = NULL, fill = "r") +
    theme_bw()

  if (grid) {
    q <- cowplot::plot_grid(p, s, labels = title)
  } else {
    q <- list(pearson = p, spearman = s)
    names(q) <- c("pearson", "spearman")
  }

  return(q)
}