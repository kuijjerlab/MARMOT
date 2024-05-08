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

  # concatenate the two matrices into a list
  corr_res <- list(corr_pears, corr_spear)
  names(corr_res) <- c("pearson", "spearman")

  return(corr_res)
}

#' @name format_fct_corr
#'
#' @description Function to format the factor correlations results for plotting.
#'
#' @inheritParams fct_corr
#' @param corr_res A list of correlation matrices. Expects output of
#' \code{\link{fct_corr}}.
#'
#' @returns A long data frame ready for plotting.
#'
#' @export
#'
format_fct_corr <- function(corr_res) {
  # sanity checks
  # check that input is a list
  if (!is.list(corr_res)) {
    stop("Invalid input format. This function expects a list of two matrices as input.") # nolint
  }

  # format matrices as data frames
  corr1 <- reshape2::melt(corr_res[[1]])
  corr1$method <- names(corr_res)[1]
  corr2 <- reshape2::melt(corr_res[[2]])
  corr2$method <- names(corr_res)[2]

  # concatenate the two methods
  corr <- rbind(corr1, corr2)

  return(corr)
}

#' @name plot_fct_corr
#'
#' @description Plots heatmaps of factor correlations.
#'
#' @inheritParams plot_data_dim
#' @param corr_res A data frame with the correlation values.
#' Expects output of \code{\link{fct_corr}} with \code{as_dat_frame = TRUE}.
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

plot_fct_corr <- function(corr_res, grid = TRUE, colours = NULL,
                          title = NULL, ...) {
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

  corr  <- corr_res

  # Plot heatmap
  p <- ggplot(corr, aes(x = Var2, y = Var1, fill = pearson, label = round(pearson, 2))) +
    geom_tile() +
    geom_text(color = "black") +
    scale_fill_gradient2(low = col[1], mid = "white", high = col[2], midpoint = 0) +
    facet_wrap(~cancer, nrow = 3) +
    labs(title = "Pearson", x = NULL, y = NULL, fill = "r") +
    theme_bw()

  s <- ggplot(corr, aes(x = Var2, y = Var1, fill = spearman, label = round(spearman, 2))) +
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