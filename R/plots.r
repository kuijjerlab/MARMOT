################################################################################
#################### Functions for various custom plots ########################
################################################################################

# @name plot_data_dim
# @description Barplot that shows the dimensionalities for a set of multi-omics data. Expects the output of "create_omics_list" 
# @param data Expects a character vector containing the path to two .Rda files. The data files are expected to be the output of "create_omics_list"
# @param data_labels Character label for the data. This will be used for the title of the plot. 
# @param Optional. Character vector specifying the omics to be plotted. If not specified, all the omics will be plotted. 
# @cols Optional. Character vector with colour IDs to use for the barplots. If NULL, the "Dark2" pallette from RColorBrewer will be used
# @returns Returns a ggplot.

plot_data_dim <- function(data, data_label, which_omics = NULL, cols = NULL){
  #loading data
  dat <- get(load(data[1]))
  
  #setting colours
  if(is.null(cols)){
    col <- palette("Dark2")
    col <- col[-1]
  }else{
    col <- cols
  }
  
  #filtering data
  if(!is.null(which_omics)){
    dat1 <- dat[[which_omics]]
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
    scale_fill_manual(values=col[1:nrow(df)])+
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +# Rotate x-axis labels for readability
    coord_flip()
  
  return(p)
  
}
