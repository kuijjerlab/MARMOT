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
  .check_names(models_to_compare, surv_df$label, err_msg = "elements of 'models_to_compare' exist in your data frame ") # nolint

  # select models to compare
  surv_df <- surv_df[which(surv_df$label %in% models_to_compare), ]

  # set colour palette
  if (is.null(colours)) {
    col <- RColorBrewer::brewer.pal(name = "Dark2", n = 8)
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

#' @name gsea_dotplots
#' @description Takes as input gsea results on surv associated MOFA factors
#' and the survival association pvalues for the factors and creates a dotplot.
#'
#' @inheritParams plot_data_dim
#' @param gsea_results Results of fgsea analysis. Can be a list with multiple
#' fgsea runs or a single data frame with the fgsea format. - not implemented
#' for now cos lazy. Just works with the lists now.
#' @param surv_df Output of surv_compare. Should be filtered beforehand for
#' significant factors if required.
#' @param gene_set What gene set was used for the gsea. For labelling.
#' @param net_metric What net metric was included in the model.
#' @param title Additional string to be added to the start of the plot title.
#' Optional.
#' @param n_path Number of pathways to plot. Will select the top nPath pathways
#' sorted by pvalue. If NULL, all pathways will be plotted.
#' @param thresh Pvalue (in -log10 scale) based on which to select pathways for
#' plotting. If NULL, all pathways will be plotted. Only use if nPath is not
#' used.
#' @param ... Any other ggplot2 parameters.
#'
#' @returns A ggplot object.
#' @export

gsea_dotplots <- function(gsea_results, surv_df, gene_set, net_metric,
                          title = NULL, n_path = 20, thresh = NULL,
                          colours = NULL, ...) {
  # set colours
  if (is.null(colours)) {
    col <- palette("Dark2")
    col <- col[c(3, 4)]
  } else {
    col <- colours
  }


  #df <- data_frame(pathway = character(), pval = numeric(), padj = numeric(),
  #                 log2err = numeric, factor = character(), pathway_db = character(),
  #                net_metric = character(), logp = numeric(),
  #                 logp_surv = numeric())
  
  #surv <- surv_df[which(surv_df$label == net_metric),]
  
  
  fct <- names(gsea_results)
    #surv_fct <- surv[which(surv$factor == fct),]
    #df <- gsea_results[[i]] #if using fgsea package
    df <- gsea_results[[i]]@result
    df$factor <- fct
    df$gene_set <- gene_set
    df$net_metric <- net_metric
    #df$logp <- as.numeric(-log10(df$padj)) #if using fgsea package
    df$logp <- as.numeric(-log10(df$p.adjust))
    #df$logp_surv <- as.numeric(surv_fct$logp)
    
    #order pathway fct by the gsea signif
    # df$pathway <- factor(df$pathway, levels = unique(df$pathway[order(df$logp)])) #if using fgsea package
    df$ID <- factor(df$ID, levels = unique(df$ID[order(df$logp)])) 
    if(!is.null(n_path)){
      df <- df[order(df$logp, decreasing = T),]
      df <- df[1:n_path,]
    }else if(!is.null(thresh)){
      df <- df[which(df$logp >= thresh),]
    }
    
    #p <- ggplot(data = df, aes(x = logp, y = pathway, color=logp))+ #if using fgsea package
    p <- ggplot(data = df, aes(x = logp, y = ID, color = NES)) +  
      geom_point(data=df, aes(size = abs(NES)*25))+
      scale_size(range = c(20, 150),name = "abs(NES)", guide = "none") +
      scale_color_gradientn(colours = col, name = "NES") +
      labs(y = NULL, x = expression("-log"[10]*"(FDR)"), title = paste(title, gene_set, net_metric, fct, sep = " ")) +
      theme_bw()+
      theme(text = element_text(size = 100),
            legend.key.size = unit(5, 'cm'),
            legend.text = element_text(size = 100),
            axis.text.y = element_text(size = 120),
            axis.text.x = element_text(size = 100),
            panel.grid.major = element_line(color = "black"))
    
    if(i == 1){
      q <- p
    }else{
      q <- plot_grid(q, p)
    }
    
    # temp2 <- temp[[j]]
    # temp2 <- temp2[,1:4]
    # temp2$factor <- fct
    # temp2$pathway_db <- name[1]
    # temp2$net_metric <- name[2]
    # temp2$logp <- -log10(temp2$padj)
    # temp2$logp_surv <- surv$logp
    # df <- rbind(df, temp2)
  
  
  
  #legend scale
  #max_surv <- max(df$logp_surv, na.rm = T)
  #max_gsea <- max(df$logp, na.rm = T)
  
  return(q)
}