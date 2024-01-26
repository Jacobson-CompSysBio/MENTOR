
# Usage:
# subclustered_dendrogram: function name
# dis_mat: path to dissimilarity matrix
# k: number of clusters to start with (default=3)
# k_increment: number of clusters to add at each iteration (default=5)
# max_size: maximum number of observations in each cluster (default=40)

################ subclustered dendrogram ################

subclustered_dendrogram <- function(
      dis_mat,
      k = 3,
      k_increment = 5,
      max_size = 40,
      map
  ) {
  
  # Perform hierarchical clustering
  hc <- hclust(dis_mat)
  
  ################### define cutree_plus function ####################
  
  # enhance cutree function to return tibble
  cutree_plus <- function(hc, k = 3) {
    # get clusters
    clusters <- cutree(hc, k = k)
    # convert to tibble
    clusters <- tibble(
      label = names(clusters),
      cluster = clusters
    )
    # return clusters
    return(clusters)
  }
  
  #################### perform sub clustering ####################
  
  # Generate first clustering
  df <- cutree_plus(hc, k)  %>%
    group_by(cluster)  %>%
    # count the number of observations in each cluster
    mutate(n = n())  %>%
    # ungroup
    ungroup()
  
  # Get maximum cluster size
  cluster_size <- max(df$n)
  
  # extract all clusters with size <= max_size
  df_converged <- df  %>%
    dplyr::filter(n <= max_size)  %>%
    # convert cluster to character
    mutate(cluster = as.character(cluster))
  
  # append df_converged to df_list
  df_list <- list(df_converged)
  
  # Start while loop for subclustering
  while (cluster_size > max_size) {
    # Increase the number of clusters by k_increment
    k <- k + k_increment
    
    # Create clusters
    df2 <- cutree_plus(hc, k)  %>%
      # group by cluster
      group_by(cluster)  %>%
      # count the number of observations in each cluster
      mutate(n = n())  %>%
      # ungroup
      ungroup()
    
    # Merge df and df2
    df <- df  %>%
      # filter clusters with size > max_size
      dplyr::filter(n > max_size)  %>%
      # Join df and df2 by label
      inner_join(df2, by = c("label"))  %>%
      # convert cluster.x and cluster.y to character
      mutate(cluster.x = as.character(cluster.x),
             cluster.y = as.character(cluster.y))  %>%
      # reset cluster.y to start at 1
      group_by(cluster.x) %>%
      mutate(cluster.y = as.numeric(factor(cluster.y)))  %>%
      # ungroup
      ungroup()  %>%
      mutate(cluster = ifelse(n.x != n.y, 
                              paste(cluster.x, cluster.y, sep = "."),
                              cluster.x))  %>%
      # Select label and cluster columns
      # rename n.y to n and drop -n.x
      dplyr::select(label, cluster, n = n.y, -n.x)
    
    # Extract all converging clusters (with size <= max_size)
    df_converged <- df  %>%
      dplyr::filter(n <= max_size)
    
    # append converged clusters to df_list
    df_list <- append(df_list, list(df_converged))
    
    # Update cluster_size
    cluster_size <- max(df$n)
    
  }
  
  # Combine all dataframes in df_list
  sub_clusters <- bind_rows(df_list)  %>%
    # ungroup
    ungroup()  %>%
    arrange(cluster)
  
  # give an unique number to each subcluster
  cluster_to_cluster_id <- sub_clusters  %>%
    distinct(cluster)  %>%
    mutate(cluster_id = seq_along(cluster))
  
  sub_clusters <- sub_clusters %>%
    # add cluster_id
    left_join(cluster_to_cluster_id, by = "cluster")
  
  # get cluster names order in hc
  cluster_names <- cutree(hc, k = 10)  %>% names()
  
  # reoder subclusters
  clusters <- sub_clusters  %>%
    arrange(match(label, cluster_names))  %>%
    pull(cluster_id)
  
  # Convert to dendrogram
  dend <- hc  %>%
    # convert to dendrogram (hclust object)
    as.dendrogram()
  
  # we need to sort them to the order of the dendrogram:
  clusters <- clusters[order.dendrogram(dend)]
  clusters_numbers <- unique(clusters) - (0 %in% clusters)
  n_clusters <- length(clusters_numbers)
  
  #################### create color pallete ####################
  
  # select palette from RColorBrewer package
  seg_pal_name <- "Dark2"
  label_pal_name <- "Set2"
  
  # get final k
  k <- n_clusters
  
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
    branches_attr_by_clusters(clusters, values = seg_palette)  %>%
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
    if(nrow(dend_labs) != nrow(map_genes)|!(all(dend_labs$ensembl %in% map_genes$ensembl))) {
        print("ERROR: please ensure your mapping file has the same number of genes in your dissimilarity matrix")
        quit(save = "no")
    }
    dend_labs <- merge(dend_labs,map_genes,by = "ensembl",all.x = TRUE)
    dend_labs <- dend_labs %>% dplyr::select(label,x,y,col,cex)
    dend_labs <- dend_labs[order(dend_labs$x,decreasing = FALSE),]
  }
  # extract segments
  dend_seg <- ggd1$segments %>%
    # if color is NA set to darkgrey
    mutate(col = ifelse(is.na(col), "darkgrey", col))
  
  ######################### Pad labels ##############################
  
  # Groups
  groups <- sub_clusters  %>%
    dplyr::select(label, group = cluster)
  
  dend_labs2 <- dend_labs  %>%
    # join with groups to add group names
    left_join(groups)  %>%
    # get the max number of characters in the symbol column
    mutate(max_nchar = max(nchar(label), na.rm = TRUE))  %>%
    # pad the symbol column with spaces to make all symbols the same length
    # This is important to make the labels line up in the dendogram
    mutate(label = str_pad(label, max_nchar, pad = " ", side = "right"))  
  
  #################### Plot dendogram ##############################
  
  # get total number of labels
  n_labels <- nrow(dend_labs)
  
  lab_max_char <- unique(dend_labs2$max_nchar)
  # get label length
  label_length <- (lab_max_char + lab_max_char) * -0.02
  
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
      fontface = "bold"
    ) +
    # reverse y axis and add axis limits
    # This is to visualize labels at the right
    # side and to allow enough space for them
    scale_y_reverse(limit = c(1, label_length)) +
    # remove plot borders with expand and set limits to max number of labels
    scale_x_continuous(expand = c(0, 0), limits = c(-1, n_labels + 1)) +
    # Add segment colors
    scale_color_identity() +
    # Add label colors
    scale_fill_identity() +
    # flip the dendrogram in 90 degrees
    coord_flip() +
    # apply dengrogram theme
    theme_dendro() +
    # remove legend
    theme(legend.position = "none")
  
  #################### Plot dendogram ##############################
  
  result <- list(dend_labs,p)
  names(result) <- c("dendrogram_labels","dendrogram")
  return(result)
  
}
