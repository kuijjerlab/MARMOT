################################################################################
#################### Functions for various custom plots ########################
################################################################################

#' @name plot_data_dim
#' @description Barplot that shows the dimensionalities for a set of multi-omics
#' data. Expects the output of \code{\link{prepare_data}} as input.
#'
#' @param data Character vector containing the path to two .Rda files. The data
#' files are expected to be the output of \code{\link{prepare_data}}.
#' @param data_labels Optional. Character label for the data. This will be used
#' for the title of the plot.
#' @param which_omics Optional. Character vector specifying the omics to be
#' plotted. If not specified, all the omics will be plotted.
#' @param colours Optional. Character vector with colour IDs to use for plotting.
#' If NULL, the "Dark2" pallette from RColorBrewer will be used. If using custom
#' colours, please make sure the number of colours specified matches the number
#' of colours needed.
#'
#' @returns Returns a ggplot.
#' @import ggplot2
#' @export

plot_data_dim <- function(data, data_label = NULL, which_omics = NULL,
                          colours = NULL) {
  # loading data
  dat <- get(load(data[1]))

  #setting colours
  if (is.null(colours)) {
    col <- RColorBrewer::palette("Dark2")
    col <- col[-1]
  } else {
    col <- colours
  }

  #filtering data
  if (!is.null(which_omics)) {
    dat <- dat[[which_omics]]
  }

  #creating data frames for plotting
  df <- data.frame(
    omic = names(dat),
    features = sapply(dat, function(mat) dim(mat)[1]),
    samples = sapply(dat, function(mat) dim(mat)[2]),
    label = data_label
  )

  #set levels
  df$omic <- factor(df$omic, levels = df$omic)

  p <- ggplot(df, aes(x = omic, y = features, fill = omic)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = paste0("n = ",samples)), angle =-90, vjust = -0.2, size = 4) +  # Annotate with number of columns
    labs(title = df$label, y = "Number of features") +
    scale_fill_manual(values = col[seq_len(nrow(df))]) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Rotate x-axis labels for readability
    coord_flip()

  return(p)

}

#' @title
#'
#' @name plot_variance_bar
#'
#' @description
#'
#' @param 
#'
#' @returns A ggplot.
#' @import ggplot2
#' @export

plot_variance_bar <- function() {
    ggplot(data = var_df, aes(x = cancer, y = var_perc, fill = omic))+
      geom_bar(stat="identity")+
      scale_fill_manual(labels = labels[[j]], values = col[2:6])+
      ylab("Variance")+ xlab("")+
      theme_classic()+
      theme(text = element_text(size=30),
            axis.text.x = element_text(angle = 45, hjust = 1))
}

#' @name plot_clin_association
#' @description This takes the output of clin_association and plots a heatmap of
#' the association of each factor with each clinical feature.
#'
#' @inheritParams plot_data_dim
#' @param clin_assoc A data frame containing the association results.
#' Expects the output of clin_association
#'
#' @returns A heatmap.
#' @export
#' @import ggplot2

plot_clin_association <- function(clin_assoc, colours = NULL) {
  # sanity checks
  # check that log p value was calculated
  if (!is.element("logp", colnames(clin_assoc))) {
    stop("-log10(p-value) required for plotting. Please make sure to run 'clin_asociation' with logtrans = TRUE") # nolint
  }

  # get colours
  if (is.null(colours)) {
    col <- RColorBrewer::palette("Dark2")
  } else {
    col <- colours
  }

  p <- ggplot(data = clin_assoc, aes(x = Factor, y = feat, fill = logp)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = col[3], name = "-log10(FDR)") +
    labs(x = NULL, y = "clinical feature") +
    theme_classic()

  return(p)
}

#' @name surv_compare_dotplot
#'
#' @description Creates a dotplot that compares factors from two JDR models in
#' terms of their association with patient survival.
#'
#' @inheritParams plot_data_dim
#' @param surv_df A dataframe. Expects the output of \code{\link{surv_compare}}.
#' @param models_to_compare A vector with the models to compare. Should be two
#' of the model labels provided in \code{\link{surv_compare}}.
#'
#' @returns A ggplot.
#' @export
#' @import ggplot2 ggbeeswarm
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate case_when

surv_compare_dotplot <- function(surv_df, models_to_compare, colours = NULL,
                                 add_intercepts) { # implement the intercepts arguments
  # sanity checks
  # check that specified models exist
  .check_names(models_to_compare, surv_df$labels, err_msg = "elements of 'models_to_compare' exist in your data frame ") # nolint

  # check logtrans exists
  if (!is.element("logp", colnames(surv_df))) {
    stop("-log10(p-value) required for plotting. Please make sure to run 'surv_compare' with logtrans = TRUE") # nolint
  }

  # select models to compare
  surv_df <- surv_df[which(surv_df$label %in% models_to_compare), ]

  # set colour palette
  if (is.null(colours)) {
    col <- RColorBrewer::palette("Dark2")
  } else {
    col <- colours
  }

  surv_df <- surv_df %>% mutate(label = case_when(logp < 1.30103 ~ "ns",
                                label == models_to_compare[1] ~ models_to_compare[1],
                                label == models_to_compare[2] ~ models_to_compare[2]))

  # setting colours
  if (is.null(colours)) {
    cols <- c(col[6], col[8], "grey")
  } else {
    cols <- col
  }


  p <- ggplot(surv_df, aes(x = cancer, y = logp, fill = label)) +
    geom_beeswarm(shape = 21, cex = 0.7, size = 7) +
    scale_y_continuous(breaks = seq(0, 10, 1), limits = c(0, 10)) +
    scale_fill_manual(values = cols) +
    geom_hline(yintercept = 2, linetype = "dashed", col = "red") +
    geom_hline(yintercept = 1.30103, linetype = "dashed", col = "blue") +
    labs(x = NULL, y = expression("-log"[10]*"(FDR)")) +
    theme_classic() +
    theme(legend.position = c(0.9, 0.85),
          legend.background = element_rect(color = "black"),
          legend.title = element_blank(),
          text = element_text(size = 20),
          axis.text.x = element_text(angle = 45, hjust = 1))

  return(p)

}