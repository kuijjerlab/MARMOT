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
#' If NULL, the "Dark2" pallette from RColorBrewer will be used (recommended).
#' If using custom colours, please make sure the number of colours specified
#' matches the number of colours needed.
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
    col <- RColorBrewer::brewer.pal(name = "Dark2", n = 8)
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

  # get colours
  if (is.null(colours)) {
    col <- RColorBrewer::brewer.pal(name = "Dark2", n = 8)
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
  .check_names(models_to_compare, surv_df$label, err_msg = "elements of 'models_to_compare' exist in your data frame ") # nolint

  # select models to compare
  surv_df <- surv_df[which(surv_df$label %in% models_to_compare), ]

  # set colour palette
  if (is.null(colours)) {
    col <- RColorBrewer::brewer.pal(name = "Dark2", n = 8)
  } else {
    if (length(colours) != 3) {
      stop(paste0(length(colours), " colours were specified, when 3 were expected. ",
      "Please make sure you specify the correct number of colours."))
    }
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

#' @name surv_compare_tile
#'
#' @description Creates a heatmap that compares factors from two JDR models in
#' terms of their association with patient survival.
#'
#' @inheritParams surv_compare_dotplot
#'
#' @returns A ggplot.
#' @export
#' @import ggplot2 ggbeeswarm
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarise

surv_compare_tile <- function(surv_df, models_to_compare, colours = NULL,
                                 add_intercepts) { # implement the intercepts arguments
  # sanity checks
  # check that specified models exist
  .check_names(models_to_compare, surv_df$label, err_msg = "elements of 'models_to_compare' exist in your data frame ") # nolint

  # select models to compare
  surv_df <- surv_df[which(surv_df$label %in% models_to_compare), ]

  # set colour palette
  if (is.null(colours)) {
    col <- RColorBrewer::brewer.pal(name = "Dark2", n = 8)
    col <- col[c(3, 4)]
  } else {
    if (length(colours) != 2) {
      stop(paste0(length(colours), " colours were specified, when 2 were expected. ",
      "Please make sure you specify the correct number of colours."))
    }
    col <- colours
  }

  # setting colours
  if (is.null(colours)) {
    cols <- c(col[1], col[2])
  } else {
    cols <- col
  }

  # Calculate maximum logp value for each factor group within each cancer group
  max_logp <- surv_df %>%
  group_by(cancer, factor) %>%
  summarise(max_logp = max(logp))

  # Merge maximum logp values back into the dataframe
  surv_df <- merge(surv_df, max_logp, by = c("factor", "cancer"),
                   suffixes = c("", "_max"))

  # Create a column to specify color based on whether logp matches
  # the maximum within its factor group
  surv_df$logp_value <- ifelse(surv_df$logp == surv_df$max_logp, "high", "low")

  # Plot heatmap with logp values annotated
  p <- ggplot(surv_df, aes(x = label, y = factor, fill = logp_value)) +
    geom_tile(width = 0.95, height = 0.95,) +
    geom_text(aes(label = round(logp, 2)), color = "black") +
    scale_fill_manual(values = c(col[2], col[1]), guide = "legend") +
    facet_wrap(~cancer, nrow = 3) +
    labs(x = NULL, y = NULL) +
    coord_equal() +
    theme_classic()

  return(p)

}

#' @name gsea_dotplots
#' @description Takes as input gsea results on surv associated MOFA factors
#' and the survival association pvalues for the factors and creates a dotplot.
#'
#' @inheritParams plot_data_dim
#' @param gsea_results Results of GSEA analysis. Expects a data frame where each
#' row is a pathway. See \code{\link{fgsea::fgsea}} for more details.
#' @param surv_df Output of surv_compare. Should be filtered beforehand for
#' significant factors if required.
#' @param title Optional. Character string to be used as plot title.
#' @param n_path Number of pathways to plot. Will select the top nPath pathways
#' sorted by pvalue. If NULL, all pathways will be plotted.
#' @param thresh Pvalue (in -log10 scale) based on which to select pathways for
#' plotting. If NULL, all pathways will be plotted. Only use if \code{n_path}
#' is not used.
#'
#' @param ... Any other ggplot2 parameters.
#'
#' @import ggplot2
#' @returns A ggplot object.
#' @export

gsea_dotplots <- function(gsea_results, surv_df, gene_set = NULL, title = NULL,
                          n_path = 20, thresh = NULL, colours = NULL, ...) {
  # sanity checks
  # check that either n_path or thresh are set
  if(!is.null(n_path) && !is.null(thresh)) {
    stop("Both 'n_path' and 'thresh' have value. Please set one or the other to NULL.")
  }

  # set colours
  if (is.null(colours)) {
    col <- RColorBrewer::brewer.pal(name = "Dark2", n = 8)
    col <- col[c(3, 4)]
  } else {
       if (length(colours) != 2) {
      stop(paste0(length(colours), " colours were specified, when 2 were expected. ",
      "Please make sure you specify the correct number of colours."))
    }
    col <- colours
  }

  df <- data.frame(gsea_res)
  df$logp <- -log10(gsea_results$padj)

  #order pathway fct by the gsea signif
  df$pathway <- factor(df$pathway, levels = unique(df$pathway[order(df$logp)]))

  # only keep top 20 or ones above threshhold
  if (!is.null(n_path)) {
    df <- df[1:n_path, ]
  } else {
    df <- df[which(df$logp >= thresh), ]
  }

  p <- ggplot(data = df, aes(x = logp, y = pathway, color = NES)) +
   geom_point(data = df, aes(size = abs(NES) * 25)) +
   scale_size(range = c(20, 150), name = "abs(NES)", guide = "none") +
   scale_color_gradientn(colours = col, name = "NES") +
   labs(y = NULL, x = expression("-log"[10] * "(FDR)"), title = title) +
   theme_bw() +
   theme(text = element_text(size = 100),
        legend.key.size = unit(5, "cm"),
        legend.text = element_text(size = 100),
        axis.text.y = element_text(size = 120),
        axis.text.x = element_text(size = 100),
        panel.grid.major = element_line(color = "black"))

  return(p)
}

#' @name volcano_plot
#' @description Plots a volcano plot from a limma model
#' @param limma TopTable of limma output
#' @param labels If true, the significant features will be labeled. Not
#' recommended if there are a lot of significant features
#' @param round_to Integer. What should the x axis be rounded to. e.g. 
#' nearest 5s, 10s, 100s etc. Default rounds to the nearest 10s.
#' @param signif_thresh P-value threshold for significance. Default 0.05.
#'
#' @returns A volcano plot
#' @export
#' @import ggplot2

volcano_plot <- function(limma, labels = FALSE, round_to=10, signif_thresh = 0.05){
  #set axis parameters (breaks & labels)
  breaks <- c(seq(plyr::round_any(min(limma$logFC),round_to), plyr::round_any(max(limma$logFC), round_to), 
                  plyr::round_any(max(abs(limma$logFC)),round_to)/4))
  limits <- c(plyr::round_any(min(limma$logFC),round_to), 
              plyr::round_any(max(limma$logFC),round_to))

  #set fold change limits
  #set to be 1/10 of the max FC value
  FClimit <- plyr::round_any(max(abs(limma$logFC)),round_to)/10

  col <- palette("Dark2")
  #setting rules for colouring & shading
  limma <- limma %>% dplyr::mutate(gene_type = case_when(logFC >= FClimit & adj.P.Val <= signif_thresh ~ "up",
                                                  logFC <= -FClimit & adj.P.Val <= signif_thresh ~ "down",
                                                  TRUE ~ "ns"))

  cols <- c("up" = col[4], "down" = col[3], "ns" = "grey") 
  sizes <- c("up" = 2, "down" = 2, "ns" = 1) 
  alphas <- c("up" = 1, "down" = 1, "ns" = 0.5)

  #pdf(file="volcano_test.pdf")
  p <- ggplot(data=limma, aes(x = logFC,y = -log10(adj.P.Val), fill = gene_type,    
                         size = gene_type, alpha = gene_type)) + 
    geom_point(shape = 21, colour = "black")+
    geom_hline(yintercept = -log10(signif_thresh),
               linetype = "dashed") + 
    geom_vline(xintercept = c(-FClimit, FClimit),
               linetype = "dashed") +
    scale_x_continuous(breaks = breaks, # Modify x-axis tick intervals    
                      limits = limits)+
    theme_classic()+
    scale_fill_manual(values = cols) + # Modify point colour
    scale_size_manual(values = sizes) + # Modify point size
    scale_alpha_manual(values = alphas)
  #dev.off()

  #determining features to be labeled
  if(labels){
    limma$genes <- rownames(limma)
    top_feat<- limma[which(limma$adj.P.Val <= signif_thresh),]
    top_feat <- top_feat[which(abs(top_feat$logFC) >= FClimit),]

    sig_genes <- limma %>%
      filter(genes %in% top_feat$genes)

    p <- p + geom_label_repel(data = sig_genes, # Add labels last to appear as the top layer  
                       aes(label = genes),
                       force = 2,
                       nudge_y = 0.25)
  }

  return(p)
}