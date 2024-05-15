
# Usage:
# dendogram: function name
# dis_mat: path to dissimilarity matrix
# k: number of clusters (default=3)

################ original dendrogram ################

dendrogram <- function(dis_mat,k = 3,map) {
  
  #################### perform hierarchical clustering ####################
  
  # perform hierarchical clustering
  hc <- hclust(dis_mat)
  # get clusters from dendogram
  clusters <- cutree(hc, k = k)
  # convert to dendrogram
  dend <- hc  %>%
    # convert to dendrogram (hclust object)
    as.dendrogram()
  # sort to same order as dendrogram
  clusters <- clusters[order.dendrogram(dend)]
  # get number of clusters (k) from dendogram
  k <- length(unique(clusters))
  
  #################### create color pallete ####################
  
  # select palette from RColorBrewer package
  seg_pal_name <- "Dark2"
  label_pal_name <- "Set2"
  # add more colors to segments palette :
  if (k > 8) {
    # get segment color palette
    seg_palette <- brewer.pal(8, seg_pal_name)
    # expand number of colors to match total clusters
    seg_palette <- colorRampPalette(seg_palette)(k)
    # if k < 8 use the default color palette
  } else {
    # get segment color palette
    seg_palette <- brewer.pal(k, seg_pal_name)
  }
  # add more colors to labels palette :
  if (k > 8) {
    # get label color palette
    label_palette <- brewer.pal(8, label_pal_name)
    # expand number of colors to match total clusters
    label_palette <- colorRampPalette(label_palette)(k)
    # if k < 8 use the default color palette
  } else {
    # get label color palette
    label_palette <- brewer.pal(k, label_pal_name)
  }
  # shuffle sequence 1-k to prevent rainbow scheme
  set.seed(123)
  seg_palette <- seg_palette[sample(1:k)]
  set.seed(123)
  label_palette <- label_palette[sample(1:k)]
  # create a cluster to color mapping
  key_value <- tibble(cls = unique(clusters), pal = label_palette[1:k])
  # convert clusters to tibble
  cls_df <- tibble(cls = clusters)
  # join to get palette to add color to clusters
  cls_df <- cls_df %>% left_join(key_value, by = "cls")
  # get label palette as vector
  label_palette2 <- cls_df$pal
  
  #################### Apply color to dendogram ####################
  
  dend2 <- dend %>%
    # apply color to branches
    branches_attr_by_clusters(clusters, values = seg_palette) %>%
    # apply color to labels
    dendextend::set("labels_colors", label_palette2)
  
  ############### Convert to ggplot friendly table ##################
  
  # convert to ggplot friendly format
  ggd1 <- as.ggdend(dend2)
  # extract labels
  dend_labs <- ggd1$labels
  # map to symbols
  if(!is.null(map)) {
    names(dend_labs)[which(names(dend_labs) == "label")] <- "ensembl"
    map_genes <- suppressMessages(read_tsv(map, col_names = TRUE, show_col_types = FALSE))
    names(map_genes) <- c("ensembl","label")
    if(nrow(map_genes) > nrow(dend_labs)) {
      extra_genes <- map_genes$ensembl[!(map_genes$ensembl %in% dend_labs$ensembl)]
      cat(paste0('\n\nWARNING: ',length(extra_genes),' genes from users map.txt file were not found in multiplex or were not in geneset.txt file; genes removed from map.txt file:\n'))
      cat(paste0(extra_genes,collapse = "\n"))
    }
    if(nrow(map_genes) < nrow(dend_labs)) {
      missing_genes <- dend_labs$ensembl[!(dend_labs$ensembl %in% map_genes$ensembl)]
      cat(paste0('\n\nWARNING: ',length(missing_genes),' genes from dendrogram visualization that were missing in users map.txt file:\n'))
      cat(paste0(missing_genes,collapse = "\n"))
    }
    dend_labs <- merge(dend_labs,map_genes,by = "ensembl",all.x = TRUE)
    dend_labs <- dend_labs %>% dplyr::select(ensembl,label,x,y,col,cex)
    dend_labs <- dend_labs[order(dend_labs$x,decreasing = FALSE),]
  }
  # extract segments
  dend_seg <- ggd1$segments %>%
    # if color is NA set to darkgrey
    mutate(col = ifelse(is.na(col),"darkgrey",col))
  
  ######################### Pad labels ##############################
  
  dend_labs2 <- dend_labs  %>%
    # get max number of characters in label
    mutate(max_nchar = max(nchar(as.character(label)),na.rm = TRUE))  %>%
    # pad label to make all labels same length
    mutate(label = str_pad(label,max_nchar,pad = " ",side = "right"))
  
  #################### Plot dendogram ##############################
  
  # get total number of labels
  n_labels <- nrow(dend_labs)
  # get label length
  label_length <- unique(dend_labs2$max_nchar) / 40 * -1
  # build base ggplot
  p <- ggplot() +
    # add segments for dendrogram
    geom_segment(
      data = dend_seg,
      aes(x = x, y = y, xend = xend, yend = yend, color = col), 
      linewidth = 1
    ) +
    # add labels for dendrogram
    geom_label(
      data = dend_labs2,
      aes(x, y, label = label, fill = col),
      color = 'black',
      hjust = 0,
      family = 'mono',
      fontface = "bold",
      label.padding = unit(0.125, "lines")
    ) +
    # reverse y axis and add axis limits
    scale_y_reverse(limit = c(max(dend_seg$y + 0.01), label_length)) +
    # remove plot borders with expand and set limits to max labels
    scale_x_continuous(expand = c(0, 0),limits = c(-1,n_labels + 1)) +
    # pass segment colors
    scale_color_identity() +
    # pass label colors
    scale_fill_identity() +
    # flip the dendrogram 90 degrees
    coord_flip() +
    # apply dengrogram theme
    theme_dendro() +
    # remove legend
    theme(legend.position = "none")
  
  #################### Return dendogram ##############################
  
  result <- list(dend_labs,dend2,p)
  names(result) <- c("dendrogram_labels","dendrogram","plot")
  return(result)
  
}

