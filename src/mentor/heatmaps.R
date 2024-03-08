
# Usage:
# heatmap: function name
# log2fc: path to log2fc file

################ heatmap function ################

heatmap <- function(heatmap,dend_labs,reordercols,p_cutoff,squish_bounds) {
  
  # read in logfc table (must be a tsv with columns: label, log2fc)
  heat_labs <- suppressMessages(read_tsv(heatmap,col_names = TRUE, show_col_types = FALSE))
  # changing first column name to "label" to match dend_labs
  names(heat_labs) <- c("label","value","source")
  if(any(duplicated(heat_labs))) {
    cat("\nWARNING: duplicated rows in heatmap table; make sure all rows are unique!")
  }
  if(reordercols) {
    cat("\n\nre-ordering the columns of the heatmap table by clustering")
    # do the rearrangement clustering
    heat_labs_no_dups <- heat_labs %>% group_by(source,label) %>% top_n(1,abs(value)) %>% distinct()
    df <- spread(heat_labs_no_dups, source, value)
    df[is.na(df)] <- 0
    columnar_data <- t(df[, 2:ncol(df)])
    distance_mat <- dist(columnar_data, method = 'euclidean')
    clusters <- hclust(distance_mat)
    ordered_labels <- clusters$order
    row_order <- row.names(columnar_data)[clusters$order]
    heat_labs <- heat_labs %>% arrange(factor(source, levels=row_order))
  }
  # create order column on dend_labs
  dend_labs$row_order <- 1:nrow(dend_labs)
  # loop through different sources and create data frame for heatmap sources
  heatmap_list <- list()
  for(i in 1:length(unique(heat_labs$source))) {
    # if dealing with first source set x = 0
    if(i == 1) {
      x <- 0
    }
    # filter to specific source
    heat_labs_source <- heat_labs %>% dplyr::filter(source == unique(heat_labs$source)[i])
    # merge with heat_labs
    heat_labs_source <- merge(heat_labs_source,dend_labs,by = "label",all.y = TRUE)
    # extract na rows
    na_rows <- heat_labs_source[is.na(heat_labs_source$source),]
    # check for duplicates and take max abs value and
    heat_labs_source <- heat_labs_source %>% group_by(label) %>% top_n(1,abs(value)) %>% bind_rows(.,na_rows)
    # reorder heat labels
    heat_labs_source <- heat_labs_source[order(heat_labs_source$row_order,decreasing = FALSE),]
    # reset the source column
    heat_labs_source$source <- unique(heat_labs$source)[i]
    # adding x coordinates for heatmap
    heat_labs_source$x <- x
    # adding y coordinates for heatmap
    heat_labs_source$y <- 1:(nrow(heat_labs_source))
    # add new data to list
    heatmap_list[[i]] <- heat_labs_source
    # increment x
    x <- x + 1
  }
  # collapse list
  heat_labs <- do.call("rbind",heatmap_list)
  # add a factor column for sources where values are all 1 or NA
  heat_labs$factor <- do.call("c",lapply(unique(heat_labs$x),function(column){
    # filter to current source
    heats <- heat_labs[heat_labs$x == column,]
    # check if all values for a given source are all 1 or NA
    if(all(heats$value == 1|is.na(heats$value))) {
      # if all values are 1 or NA then set NA values to 0
      value <- heats$value
      return(case_when(value == 1 ~ 1,is.na(value) ~ 0))
    } else {
      # else set all values to 0
      return(rep(0,length(heats$value)))
    }
  }))
  # change levels of factor column to "absence" or "presence"
  heat_labs$factor <- factor(heat_labs$factor,levels = c(0,1),labels = c("absent","present"))
  # add a p-value column for sources where all values are between 0 and 1 inclusive
  #heat_labs$pvalue <- do.call("c",lapply(unique(heat_labs$x),function(column){
  #  # filter to current source
  #  heats <- heat_labs[heat_labs$x == column,]
  #  # check if all values fall between 0 and 1 inclusive and that all values for factor are "absent"
  #  if(all((heats$value >= 0 & heats$value <= 1)|is.na(heats$value)) & !any(heats$factor == "present")) {
  #    # if all values fall between 0 and 1 inclusive and all values for factor are "absent"
  #    x <- heats$value
  #    # set yes p-values less than 0.05 as yes and greater than or equal to 0.05 as no
  #    return(case_when(
  #      is.na(x) ~ "na",
  #      x < p_cutoff ~ "yes",
  #      x >= p_cutoff ~ "no"
  #    ))
  #  } else {
  #    # otherwise set all p-values to "na"
  #    return(rep("na",length(heats$value)))
  #  }
  #}))
  # create heatmap
  heat <- ggplot(data = heat_labs,aes(x,y)) + 
    # set all tiles to be grey initially
    geom_tile(aes(fill = !!sym("value")),colour = "grey50",fill = "grey50",show.legend = FALSE) + 
    # set the labels for the x-axis
    scale_x_continuous(
      # labels for x-axis
      labels = unique(heat_labs$source)[!is.na(unique(heat_labs$source))],
      # breaks for x-axis labels
      breaks = unique(heat_labs$x)
    ) +
    # adjusting y-scale
    scale_y_continuous(expand = c(0, 0),limits = c(-1,nrow(dend_labs) + 1)) +
    # adjusting the theme of the ggplot
    theme(
      # removing axis lines
      axis.line = element_blank(),
      # removing axis ticks
      axis.ticks = element_blank(),
      # removing axis text
      axis.text.y = element_blank(),
      # removing axis titles
      axis.title = element_blank(),
      # change angle of x-axis title
      axis.text.x = element_text(angle = 90,vjust = 0.15,hjust = 1,size = 14),
      # removing background panel
      panel.background = element_blank(),
      # adjust legend title size and position
      # legend.title = element_text(size = 8,vjust = 2)
    )
  # if there are any values that contain a decimal and fall below 0 or greater than 1
  if(any(grepl("\\.",as.character(heat_labs$value))) & any(abs(heat_labs$value) >= 1|heat_labs$value < 0)) {
    # add new scale for log2fc values
    heat <- heat + 
      # add new scale
      new_scale_fill() + 
      # add geom_tile fills for log2fc values
      geom_tile(aes(x,y,fill = value), show.legend = TRUE)
    # check if squish_bounds is null
    if(!is.null(squish_bounds)) {
      # get upper and lower bound squish values
      squish_bounds <- as.numeric(do.call("c",strsplit(squish_bounds,",")))
      # draw new scale with squish bounds parameters
      heat <- heat + scale_fill_gradient2(
        # set the value legend to log2fc (assume these values represent log2fc?)
        name = TeX("$\\log_{2}(FC)$"),
        # set low value to blue
        low = '#0017FF',
        # set 0 to white
        mid = "white",
        # set high value to red
        high = '#FF2D00',
        # set midpoint to 0
        midpoint = 0,
        # set limits based on squish bounds
        limits = c(squish_bounds[1],squish_bounds[2]),
        # squish the scales
        oob = scales::squish,
        # set NA values to grey
        na.value = "grey50"
      )
    } else {
      heat <- heat + scale_fill_gradient2(
        # set the value legend to log2fc (assume these values represent log2fc?)
        name = TeX("$\\log_{2}(FC)$"),
        # set low value to blue
        low = '#0017FF',
        # set 0 to white
        mid = "white",
        # set high value to red
        high = '#FF2D00',
        # set midpoint to 0
        midpoint = 0,
        # set NA values to grey
        na.value = "grey50"
      )
    }
    heat <- heat + theme(
        # adjust legend title size and position
        legend.title = element_text(size = 8,vjust = 2)
      )
  }
  # if there are any factor values set to "present" (when p-values not present)
  if(any(heat_labs$factor == "present")) {
    # add new scale for factor values
    heat <- heat + 
      # add new scale
      new_scale_fill() + 
      # add geom_tile fills for factor values
      geom_tile(aes(x,y,fill = factor),show.legend = TRUE) + 
      scale_fill_manual(
        # set factor values legend title to ""
        name = "",
        # set absent to transparent and present to black
        values = c("absent" = "transparent","present" = "black"),
        # override the legend specs
        guide = guide_legend(
            override.aes = list(
              # override the legend fills to grey and black 
              fill = c("grey50","black")
            )
          )
      )
  }
  # if there are any p-value values set to "yes"
  #if(any(heat_labs$pvalue == "yes")) {
  #  # add new scale for p-values
  #  heat <- heat +
  #    # add new scale
  #    new_scale_fill() +
  #    # add geom_tile fills for p-values
  #    geom_tile(aes(x,y,fill = pvalue)) +
  #    # adjust legend title size and position
  #    theme(legend.title = element_text(size = 8,vjust = 2)) +
  #    scale_fill_manual(
  #      # set the p-value legend to adjusted p-value
  #      name = "p-value (adj)",
  #      # set "yes" to black "no" to grey and "na" to transparent
  #      values = c("yes" = "black","no" = "grey50","na" = "transparent"), 
  #      # set labels of legend to "< 0.05" and ">= "0,05"
  #      labels = c(paste0("< ",p_cutoff),paste0(">= ",p_cutoff),""),
  #      # override the legend specs
  #      guide = guide_legend(
  #        override.aes = list(
  #          # override the legend fills to black grey and white (removes "na")
  #          fill = c("black","grey50","white")))
  #    )
  #}
  
  return(heat)
  
}
