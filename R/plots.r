################################################################################
#################### Functions for various custom plots ########################
################################################################################

#' @name plot_data_dim
#' @description Barplot that shows the dimensionalities for a set of multi-omics
#' data. Expects the output of \code{\link{prepare_data}} as input.
#'
#' @param data Character vector containing the path to one or two .Rda files.
#' The data files are expected to be the output of \code{\link{prepare_data}}.
#' @param data_labels Character vector with labels for the data. This will be used
#' for the title of the plot. Also used for comparison. Order should be the same
#' as \code{data}.
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

plot_data_dim <- function(data, data_labels, which_omics = NULL,
                          colours = NULL, compare = TRUE) {
  # sanity
  if (compare) {
    if (length(data) !=2) {
      stop("You must provide two data files if comparison = TRUE")
    } else {
      # loading data
      dat <- get(load(data[1]))
      dat_2 <- get(load(data[2]))
    }
  } else {
    if (length(data) !=1) {
      stop("You mut provide one data file if comparison = FALSE")
    } else {
      dat <- get(load(data))
    }
  }

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
    if (compare) {
      dat_2 <- dat_2[[which_omics]]
    }
  }

  #creating data frame for plotting
  df <- data.frame(
    omic = names(dat),
    features = sapply(dat, function(mat) dim(mat)[1]),
    samples = sapply(dat, function(mat) dim(mat)[2]),
    label = data_labels[1]
  )

  # add second data info
  if (compare) {
    df_2 <- data.frame(
    omic = names(dat_2),
    features = sapply(dat_2, function(mat) dim(mat)[1]),
    samples = sapply(dat_2, function(mat) dim(mat)[2]),
    label = data_labels[2]
    )

    # merge
    df <- rbind(df, df_2)
  }

  # oder by number of feat
  #df <- df %>% arrange(features)

  #set levels
  df$omic <- factor(df$omic, levels = unique(df$omic))

  p <- ggplot(df, aes(x = omic, y = features, fill = omic)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = paste0("n = ",samples)), angle =-90, vjust = -0.2, size = 4) +  # Annotate with number of columns
    labs(title = df$label, y = "Number of features") +
    scale_fill_manual(values = col[seq_len(nrow(df))]) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Rotate x-axis labels for readability
    coord_flip() +
    facet_wrap(~label, ncol = 2, scales = "free_x")

  return(p)

}

#' @name plot_data_distribution
#' @description Plot distribution of the data.

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
      geom_violin(stat="identity")+
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

  p <- ggplot(data = clin_assoc, aes(x = Factor, y = feat, fill = logp,
              label = round(logp, 2))) +
    geom_tile() +
    geom_text(color = "black") +
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
    scale_y_continuous(breaks = seq(0, 5, 1), limits = c(0, 5)) +
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
#' @param file_name Optional. Character string with a file name. If used, the
#' plot will be saved in a file with the specified name, as well as returned
#' as an object.
#'
#' @param ... Any other ggplot2 parameters.
#'
#' @import ggplot2
#' @returns A ggplot object.
#' @export

gsea_dotplots <- function(gsea_results, surv_df, gene_set = NULL, title = NULL,
                          n_path = 20, thresh = NULL, colours = NULL,
                          file_name = NULL, ...) {
  # sanity checks
  # check that either n_path or thresh are set
  if(!is.null(n_path) && !is.null(thresh)) {
    stop("Both 'n_path' and 'thresh' have value. Please set one or the other to NULL.")
  }

  # set colours
  if (is.null(colours)) {
    col <- RColorBrewer::brewer.pal(name = "Dark2", n = 8)
    col <- col[c(3, 4)]
    col <- append(col, "white", after = 1)
  } else {
       if (length(colours) != 3) {
      stop(paste0(length(colours), " colours were specified, when 3 were expected. ",
      "Please make sure you specify the correct number of colours."))
    }
    col <- colours
  }

  df <- data.frame(gsea_results)
  df$logp <- -log10(gsea_results$padj)

  #order pathway fct by the gsea signif
  df$pathway <- factor(df$pathway, levels = unique(df$pathway[order(df$logp)]))

  # only keep top 20 or ones above threshhold
  if (!is.null(n_path)) {
    df <- df[1:n_path, ]
  } else {
    df <- df[which(df$logp >= thresh), ]
  }

  p <- ggplot(data = df, aes(x = logp, y = pathway, color = NES, shape = factor)) +
   geom_point(data = df, aes(size = abs(NES) * 25)) +
   scale_size(range = c(20, 150), name = "abs(NES)", guide = "none") +
   scale_color_gradient2(low = col[1], mid = col[2], high = col[3],midpoint = 0,
                         name = "NES") +
   labs(y = NULL, x = expression("-log"[10] * "(FDR)"), title = title) +
   theme_bw() +
   theme(text = element_text(size = 100),
        legend.key.size = unit(5, "cm"),
        legend.text = element_text(size = 100),
        axis.text.y = element_text(size = 120),
        axis.text.x = element_text(size = 100),
        panel.grid.major = element_line(color = "black")) +
   guides(
     shape = guide_legend(override.aes = list(size = 30))  # Increase size of shapes in legend
   )
  
  if (!is.null(file_name)) {
    ggsave(p, file = file_name, ...)
  }

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

volcano_plot <- function(limma, labels = FALSE, round_to = 10, signif_thresh = 0.05) {
  #set axis parameters (breaks & labels)
  breaks <- c(seq(plyr::round_any(min(limma$logFC),round_to), plyr::round_any(max(limma$logFC), round_to), 
                  plyr::round_any(max(abs(limma$logFC)), round_to) / 4))
  limits <- c(plyr::round_any(min(limma$logFC), round_to),
              plyr::round_any(max(limma$logFC), round_to))

  #set fold change limits
  #set to be 1/10 of the max FC value
  FClimit <- plyr::round_any(max(abs(limma$logFC)), round_to) / 10

  col <- palette("Dark2")
  #setting rules for colouring & shading
  limma <- limma %>% dplyr::mutate(gene_type = case_when(logFC >= FClimit & adj.P.Val <= signif_thresh ~ "up",
                                                  logFC <= -FClimit & adj.P.Val <= signif_thresh ~ "down",
                                                  TRUE ~ "ns"))

  cols <- c("up" = col[4], "down" = col[3], "ns" = "grey") 
  sizes <- c("up" = 2, "down" = 2, "ns" = 1) 
  alphas <- c("up" = 1, "down" = 1, "ns" = 0.5)

  #pdf(file="volcano_test.pdf")
  p <- ggplot(data=limma, aes(x = logFC, y = -log10(adj.P.Val), fill = gene_type,
                         size = gene_type, alpha = gene_type)) +
    geom_point(shape = 21, colour = "black")+
    geom_hline(yintercept = -log10(signif_thresh),
               linetype = "dashed") +
    geom_vline(xintercept = c(-FClimit, FClimit),
               linetype = "dashed") +
    scale_x_continuous(breaks = breaks, # Modify x-axis tick intervals
                      limits = limits) +
    theme_classic() +
    scale_fill_manual(values = cols) + # Modify point colour
    scale_size_manual(values = sizes) + # Modify point size
    scale_alpha_manual(values = alphas)
  #dev.off()

  #determining features to be labeled
  if (labels) {
    limma$genes <- rownames(limma)
    top_feat <- limma[which(limma$adj.P.Val <= signif_thresh), ]
    top_feat <- top_feat[which(abs(top_feat$logFC) >= FClimit), ]

    sig_genes <- limma %>%
      filter(genes %in% top_feat$genes)

    p <- p + ggrepel::geom_label_repel(data = sig_genes, # Add labels last to appear as the top layer  
                       aes(label = genes),
                       force = 2,
                       nudge_y = 0.25)
  }

  return(p)
}

#' @name plot_feat_wts
#'
#' @description Function to plot feature weights.
#'
#' @inheritParams plot_data_dim
#' @param feat_wts Matrix with omic feature weights.
#' @param fct Character vector indicating the names of the factors to be plotted.
#' If \code{NULL}, all factors will be plotted.
#' @param n_feat Number of top features to label. Is overridden by \code{thresh}
#' Should a positive integer.
#' @param manual_lab Character vector of feature names to manually label.
#' Will be checked against omic rownames.
#' @param scale Logical. Whether to scale the feature weight between c(-1,1).
#' @param thresh Numeric. Weight threshold above which to label features.
#' Absolute value will be considered.
#' @param plot_distribution If \code{TRUE}, this will show a distribution of the
#' feature weights, annotating the top features. If \code{FALSE}, it will
#' plot only the top features.
#' @param ... Any other parameters of \code{ggplot} functions.
#'
#' @returns A ggplot object.
#' @export
#' @import dplyr ggplot2, tidyr
#' @importFrom magrittr %>%

plot_feat_wts <- function(feat_wts, fct = NULL, n_feat = 10, manual_lab = NULL,
                          scale = TRUE, file_name = NULL, thresh = NULL,
                          plot_distribution = FALSE, colours = NULL, ...) {
  # check that manual labels exist in the data
  if (!is.null(manual_lab)) {
    .check_names(manual_lab, rownames(feat_wts),
                 err_msg = "features marked for manual labelling exist in the omic data") #nolint
  }

  # set colours
  if (is.null(colours)) {
    col <- RColorBrewer::brewer.pal(name = "Dark2", n = 8)
    col <- col[c(3, 4)]
  } else {
       if (length(colours) != 3) {
      stop(paste0(length(colours), err_msg = " colours were specified, when 2 were expected. ",
      "Please make sure you specify the correct number of colours."))
    }
    col <- colours
  }

  # reshape for plotting
  df <- reshape2::melt(feat_wts)
  colnames(df) <- c("feature", "factor", "value")

  # get factors of interest
  if (!is.null(fct)) {
    df <- df %>% filter(factor %in% fct)
  }

  # scale by weight with highest absolute value
  if (scale) {
    df <- df %>%
           group_by(factor) %>%
           mutate(value = value / max(abs(value), na.rm = TRUE)) %>%
           ungroup()
  }

  # add labelling groups
  df$to_label <- FALSE

  # Define group of features to color according to the loading
  if (is.null(manual_lab) && n_feat > 0) {
    features <- df %>%
                group_by(factor) %>%
                top_n(n =  n_feat, abs(value)) %>%
                ungroup()
    features <- features$feature
    df <- df %>% 
            group_by(factor) %>%
            mutate(to_label = ifelse(feature %in% features, TRUE, to_label)) %>%
            ungroup()
  } else if (!is.null(manual_lab) && n_feat > 0) {
      features <- df %>%
                group_by(factor) %>%
                top_n(n = n_feat, abs(value)) %>%
                filter(feature %in% manual_lab) %>%
                ungroup()
      features <- features$feature
      df <- df %>% 
            group_by(factor) %>%
            mutate(to_label = ifelse(feature %in% features, TRUE, to_label)) %>%
            ungroup()
    } else if (!is.null(manual_lab) && n_feat == 0) {
      features <- df %>%
                group_by(factor) %>%
                filter(feature %in% manual_lab) %>%
                ungroup()
      features <- features$feature
      df <- df %>% 
            group_by(factor) %>%
            mutate(to_label = ifelse(feature %in% features, TRUE, to_label)) %>%
            ungroup()
    }

  # add anything above threshhold if specified
  if (!is.null(thresh)) {
    df <- df %>%
    group_by(factor) %>%
    mutate(to_label = ifelse(abs(value) >= thresh, TRUE, to_label)) %>%
    ungroup()
  }

  # Sort features by weight
  df <- df %>%
  group_by(factor) %>%
  arrange(feature, value) %>%
  ungroup()

  # store sign
  df$sign <- ifelse(df$value > 0, "+", "-")

  df$feature_id <- paste(df$feature, df$factor, sep = "_")
  df$feature_id <- factor(df$feature_id, levels = unique(df$feature_id[order(df$value)]))

  if (plot_distribution) {
    p <- ggplot(df, aes(x = value, y = feature_id), ...) +
      geom_point() +
      scale_y_discrete(expand = c(0.03, 0.03)) +
      labs(x = "Weight", y = "Rank", size = 10) +
      theme_minimal() +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank()
      )

    # Add labels
    if (n_feat > 0 || length(unique(df$to_label)) > 0) {
      p <- p + ggrepel::geom_text_repel(
        data = df[df$to_label != FALSE, ],
        aes(label = feature, col = to_label),
        size = 3, segment.alpha = 0.25, segment.color = "black",
        segment.size = 0.3, show.legend = FALSE, max.overlaps = Inf)
    }

    # configure axes
    if (scale) {
      p <- p + 
        coord_cartesian(xlim = c(-1, 1)) +
        scale_x_continuous(breaks = c(-1, 0, 1)) +
        expand_limits(x = c(-1, 1))
    }

    # Define dot size
    p <- p + scale_size_manual(values=c(2/2, 2*2)) + guides(size = "none")
  } else {
    # filter to only selected features
    df <- df %>% filter(to_label == TRUE)

    p <- ggplot(df, aes(x = reorder(feature, value), y = value, fill = sign)) +
      geom_bar(stat = "identity", show.legend = FALSE) +
      scale_fill_manual(values = col) +
      ylab("Weight") +
      xlab("") +
      coord_flip() +
      theme_bw() +
      theme(
        axis.title.x = element_text(color = 'black'),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=rel(1.1), hjust = 1, color = 'black'),
        axis.text.x = element_text(color = 'black'),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(color = "black"),
        legend.key = element_rect(fill = "transparent"),
        plot.title = element_text(size = rel(1.3), hjust = 0.5),

        # facets
        #strip.text = element_text(size = rel(1.2)),
        #panel.background = element_blank(),
        #panel.spacing = unit(1, "lines"),

        # gridlines
        panel.grid.major.y = element_blank(),
        ) +
        facet_wrap(~factor, nrow = 1, scales = "free")
  }

  # Facet if multiple factors
  if (ncol(feat_wts) > 1) {
    p <- p + facet_wrap(~factor, nrow = 1, scales = "free")
  }

  if (!is.null(file_name)) {
    ggsave(p, file = file_name, ...)
  }

  return(p)
}

#' @name top_surv_factors_km
#' @description Plots Kaplan-Meier curves for the significantly associated
#' survival factors (identified in surv_compare)
#' @param surv Survival information as data frame.
#' @param title Character string. Title of the plot.
#' @param conf.int Whether or not to show confidence intervals. Will be passed to ggsurvplot.
#' @param factor Vector containing the factor for which to do this.
#' Only one factor at a time, call mutiple times for multiple factors.
#' @param minprops Minimum proportion of samples per group for the split into
#' survival groups.
#' @returns A list containing KM plots for each significantly survival associated factor
#' @export
#' @import survminer ggplot2 survival
surv_factor_km <- function(surv, factor, title, conf_int = FALSE,
                           minprops) {

  # overlap samples
  # Define the sets
  sets <- list(
    surv_samples = surv$sample_id,
    factor_samples = names(factor)
  )
  # Compute the intersection of sets
  samples <- .overlap_sets(sets, 
                           err_msg = "There is no overlap between the samples in the survival and factor data. Please ensure at least some samples are common.")

  factor <- factor[samples]
  surv <- surv[which(surv$sample_id %in% samples), ]

  #get cutpoint
  df_list <- list()
  fit_list <- list()
  legend_labels <- c()
  for (minprop in minprops) {
    df <- .fct_cutpoint(factor, surv, minprop = minprop)
    df$FactorValue <- paste0(df$FactorValue, "_", minprop)
    fit <- survival::survfit(survival::Surv(time, event) ~ FactorValue, df)
    df_list <- append(df_list, list(df))
    fit_list <- append(fit_list, list(fit))
    legend_labels <- c(legend_labels, unique(df$FactorValue))
  }

  names(df_list) <- paste0("minprop_", minprops)
  names(fit_list) <- paste0("minprop_", minprops)

  #colors to use
  col <- palette("Dark2")
  high_cols <- colorRampPalette(c(col[1], "white"))(length(fit_list)+1)
  low_cols <- colorRampPalette(c(col[2], "white"))(length(fit_list)+1)
  cols <- unlist(mapply(c, low_cols, high_cols, SIMPLIFY = FALSE))
  cols <- head(cols, -2)
  names(cols) <- NULL

  km <- survminer::ggsurvplot_combine(fit_list, data = df_list,
    conf.int = conf_int,
    pval = TRUE,
    pval.method = TRUE,
    pval.coord = c(max(df$time) / 4, 0),
    fun = function(y) y * 100,
    palette = cols,
    strata = "group",
    legend = "top",
    xlab = "Time to last follow-up", ylab = "Survival probability (%)",
    legend.labs = legend_labels,
    title = title
  )

  km <- km$plot +
      theme(
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        plot.title = element_text(size = 25)
  )

  return(km)
}

#' Plot Data Distributions with Violin Plots
#'
#' This function creates violin plots to visualize the distributions of data 
#' from one or two lists of omics data. It allows for custom coloring and 
#' labeling of the different datasets.
#'
#' @param omic_list A list of matrices representing the first set of omics data
#' to be plotted.
#' @param omic_list2 An optional second list of matrices representing a second
#' set of omics data to be plotted.
#' @param labels A character vector of length 2 specifying labels for the two 
#' omics lists (if \code{omic_list2} is provided). If only one list is provided, 
#' only one label is required.
#' @param colours A character vector of length 2 specifying the fill colors for 
#' the violin plots. If not provided, default colors from the "Dark2" palette 
#' will be used.
#' @param title A character string specifying the title of the plot.
#'
#' @return A ggplot object representing the violin plot of the data distributions.
#'
#' @details This function can handle one or two lists of omics data, converting 
#' them into a long format data frame for visualization with ggplot2. If two 
#' lists are provided, they are combined into a single data frame before plotting.
#'
#' @importFrom ggplot2 ggplot aes geom_violin scale_fill_manual labs theme_minimal 
#' theme
#' @importFrom RColorBrewer brewer.pal
#'
#' @examples
#' omic_list1 <- list(matrix1 = matrix(rnorm(100), nrow = 10, ncol = 10))
#' omic_list2 <- list(matrix2 = matrix(rnorm(100), nrow = 10, ncol = 10))
#' plot_data_distributions(omic_list1, omic_list2, labels = c("Set 1", "Set 2"))
#'
#' @export

plot_data_distributions <- function(omic_list, omic_list2 = NULL, labels = NULL,
                                    colours = NULL, title = NULL) {
  #sanity
  if (!is.null(omic_list2) && length(labels) != 2) {
    stop("Please provide labels for both lists of omics.")
  } 

  # set colours
  if (is.null(colours)) {
    col <- RColorBrewer::brewer.pal(name = "Dark2", n = 8)
    col <- col[c(8, 6)]
  } else {
       if (length(colours) != 2) {
      stop(paste0(length(colours), err_msg = " colours were specified, when 2 were expected. ",
      "Please make sure you specify the correct number of colours."))
    }
    col <- colours
  }

  # create df for plotting
  df <- .convert_list_to_df(omic_list, label = labels[1])

  if (!is.null(omic_list2)) {
    df_2 <- .convert_list_to_df(omic_list2, label = labels[2])
    df <- rbind(df, df_2)
  }

  p <- ggplot(df, aes(x = omic, y = value, fill = pca)) +
    geom_violin() +
    scale_fill_manual(values = col) +
    labs(x = NULL, y = NULL, title = title) +
    theme_minimal() +
    theme(axis.text = element_text(size = 25),
      legend.position = "right")
  
  return(p)
}