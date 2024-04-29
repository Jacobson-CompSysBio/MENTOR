
# function to create polar dendrogram
circlize <- function(dend_list,heatmap,heatmap_factor,clusters,cluster_label_size,split,col_fun1,labels_size,group_colors,max_height,track_height,highlight_index,highlight_color) {
  
    # clear circos
    circos.clear()
    # setting up global parameters for circos plot
    circos.par(
        "canvas.xlim" = c(-1, 1),
        "canvas.ylim" = c(-1, 1),
        gap.degree = c(rep(5,nrow(clusters) - 1),8)
    )
    # create track height
    track_height <- as.numeric(do.call("c",strsplit(track_height,",")))
    # check if continuous heatmap is present
    if(!is.null(heatmap)) {
        # add continuous heatmap track to circos plot
        circos.heatmap(
            heatmap,
            col = col_fun1,
            split = split,
            na.col = "grey",
            cluster = FALSE,
            rownames.side = "outside",
            rownames.cex = labels_size,
            cell_width = rep(1,nrow(heatmap)),
            cell.border = NA,
            track.height = track_height[1]
        )
    }
    # adding cluster numbers to outside of circos plot
    max_label <- max(nchar(rownames(heatmap)))
    circos.track(
        track.index = get.current.track.index(),
        panel.fun = function(x, y) {
            circos.text(
                CELL_META$xcenter,
                CELL_META$cell.ylim[2] + 3, # (1.75 - (1/max_label))
                paste0("C",CELL_META$sector.index),
                facing = "inside",
                cex = cluster_label_size,
                col = clusters[as.numeric(CELL_META$sector.index) + 1,]$color,
                adj = c(0.5, 0), niceFacing = TRUE
            )    
        },bg.border = NA
    )
    # check again if continuous heatmap is present
    if(!is.null(heatmap)) {
        # add color blocks for continuous heatmap (group labels)
        circos.track(
            track.index = get.current.track.index(),
            panel.fun = function(x, y) {
                if(CELL_META$sector.numeric.index == nrow(clusters)) {
                    for(i in 1:ncol(heatmap)) {
                        circos.rect(
                            CELL_META$cell.xlim[2] + convert_x(5.5,"mm"),i - 1, # might have to make 5.5 and 15.5 dynamic
                            CELL_META$cell.xlim[2] + convert_x(15.5,"mm"),i,
                            col = group_colors[i], border = NA
                        )
                    }
                }
            },bg.border = NA,track.height = track_height[1]
        )
        group_colors <- group_colors[(ncol(heatmap) + 1):length(group_colors)]
    }
    # check if factor heatmap is present
    if(!is.null(heatmap_factor)) {
        # add factor heatmap tack to circos plot
        circos.heatmap(
            heatmap_factor,
            col = c("1" = "black"),
            split = split,
            na.col = "grey",
            cluster = FALSE,
            cell_width = rep(1,nrow(heatmap)),
            cell.border = NA,
            track.height = track_height[2]
        )
        # add color blocks for factor heatmap (group labels)
        circos.track(
            track.index = get.current.track.index(),
            panel.fun = function(x, y) {
                if(CELL_META$sector.numeric.index == nrow(clusters)) {
                    for(i in 1:ncol(heatmap_factor)) {
                        circos.rect(
                            CELL_META$cell.xlim[2] + convert_x(4,"mm"),i - 1, # might have to make 4 and 11 dynamic
                            CELL_META$cell.xlim[2] + convert_x(11,"mm"),i,
                            col = group_colors[i], border = NA
                        )
                    }
                }
            },bg.border = NA,track.height = track_height[2]
        )
    }
    # add dendrogram track to circos plot
    circos.track(
        ylim = c(0,max_height),
        panel.fun = function(x,y) {
            sector.index = get.cell.meta.data("sector.index")
            dend = dend_list[[as.numeric(sector.index) + 1]] #  + 1
            circos.dendrogram(
                dend,
                max_height = max_height,
                facing = "outside"
            )
        },bg.border = NA,track.height = track_height[3]
    )
    # check if the index for highlighting is not null
    if(!is.null(highlight_index)) {
        if(!is.null(heatmap) & !is.null(heatmap_factor)) {
            total_tracks <- 4
        } else {
            total_tracks <- 3
        }
        # draw the highlight for the user input sector
        draw.sector(
            get.cell.meta.data("cell.start.degree", sector.index = highlight_index) + 1,
            get.cell.meta.data("cell.end.degree", sector.index = highlight_index) - 1,
            rou1 = get.cell.meta.data("cell.top.radius", track.index = 1) + 0.1,
            rou2 = get.cell.meta.data("cell.bottom.radius", track.index = total_tracks) - 0.05,
            col = adjustcolor(highlight_color,0.20),
            border = NA
        )
    }
    # if(!is.null(zoom)) {
    #     placeholder for adding a zoom option to visualize a single section
    # }
    circos.clear()
  
}

# function to prepare data for polar dendrogram
polar_dendrogram <- function(dend_labs,dend2,heatmap,squish_bounds,cluster_label_size,labels_size,group_colors,track_height,highlight_index,highlight_color,legend_title,plot_file,height = 20,width = 20) {
  
    # # get columns want for heatmap
    # heatmap <- heatmap[,c("label","value","source")]
    # create data frame with dendrogram labels and heatmap
    heatmap <- merge(dend_labs,heatmap,all.y = TRUE,by = "label")
    # order by the original row order
    heatmap <- heatmap[order(heatmap$x,heatmap$y,decreasing = FALSE),]
    clusters <- unique(heatmap$col.x)
    clusters <- data.frame(
      "color" = clusters
    )
    clusters$cluster <- 0:(nrow(clusters) - 1)
    # change column names
    split <- unique(heatmap[,c("label","cluster")])$cluster
    heatmap <- heatmap[,c("label","value","source","cluster","row_order.x")]
    # change label to factor
    heatmap$label <- factor(heatmap$label,levels = unique(heatmap$label))
    # change source to factor
    heatmap$source <- factor(heatmap$source,levels = unique(heatmap$source))
    # transform data before running acast
    heatmap_list <- heatmap %>% group_split(source)
    # loop through list and give TRUE/FALSE
    factor_cols <- do.call("c",lapply(heatmap_list,function(x) {
        if(all(is.na(x$value)|x$value == 1)) {
            return(TRUE)
        } else {
            return(FALSE)
        }
    }))
    if(any(factor_cols)) {
        heatmap_factor <- data.frame(do.call("rbind",heatmap_list[factor_cols]))
        heatmap_factor$value <- as.character(heatmap_factor$value)
        heatmap_factor <- acast(heatmap_factor,label ~ source,var = "value")
        heatmap <- data.frame(do.call("rbind",heatmap_list[!factor_cols]))
    }
    # grab user defined group colors
    if(!is.null(group_colors)) {
        group_colors <- do.call("c",strsplit(group_colors,","))  
    }
    # transform data for polar dendrogram
    if(!is.null(heatmap)) {
        heatmap <- acast(heatmap,label ~ source,var = "value") 
        if(ncol(heatmap_factor) >= 1) {
            if(is.null(group_colors)) {
                group_colors <- primary.colors(n = ncol(heatmap) + ncol(heatmap_factor) + 1,no.white = TRUE)
            }
        } else {
            if(is.null(group_colors)) {
                group_colors <- primary.colors(n = ncol(heatmap) + 1,no.white = TRUE)
            }
        }
    } else {
        if(is.null(group_colors)) {
            group_colors <- primary.colors(n = ncol(heatmap_factor) + 1,no.white = TRUE)  
        }
    }
    # remove black from the default colors
    group_colors <- group_colors[which(group_colors != "#000000")]
    # color scale for continuous columns
    if(!is.null(squish_bounds)) {
        squish_bounds <- as.numeric(do.call("c",strsplit(squish_bounds,",")))
        col_fun1 <- colorRamp2(c(squish_bounds[1],0,squish_bounds[2]),c("blue","white","red"))  
    } else {
        col_fun1 <- colorRamp2(c(-2,0,2),c("blue","white","red"))
    }
    # change clusters to match sub-clustering!
    dend_list <- get_subdendrograms(dend2,nrow(clusters))
    max_height <- max(sapply(dend_list,function(x) attr(x, "height")))
    pdf(plot_file,height = height,width = width)
    circlize(
        dend_list = dend_list,
        heatmap = heatmap,
        heatmap_factor = heatmap_factor,
        clusters = clusters,
        cluster_label_size = cluster_label_size, # cluster_label_size = 3
        split = split,
        col_fun1 = col_fun1,
        labels_size = labels_size, # labels_size = 0.75
        group_colors = group_colors, # rev(brewer.pal(8,"Set1")[c(1:5,8)])
        max_height = max_height,
        track_height = track_height, # track_height = "0.2,0.2,0.2"
        highlight_index = highlight_index,
        highlight_color = highlight_color # highlight_color = "#34EBDC"
    )
    if(!is.null(heatmap)) {
        cont_legend <- Legend(
            at = c(squish_bounds[1],0,squish_bounds[2]),
            col_fun = col_fun1,
            title_position = "topleft",
            title = legend_title[1], # callt his somthing else
            title_gp = gpar(fontface = 1,cex = 1.5),
            labels_gp = gpar(fontface = 1,cex = 1.25),
            title_gap = unit(10,"mm"),
            grid_height = unit(2,"inches"),
            grid_width = unit(0.25,"inches")
        )
    }
    legend(
        # x = 1.25,y = 0,
        x = 1.15,y = 0.25,
        legend = c(colnames(heatmap_factor),colnames(heatmap)),
        fill = c(group_colors[(ncol(heatmap)+1):length(group_colors)],group_colors[1:ncol(heatmap)]),
        title = legend_title[2], # make dynamic
        xjust = 0,
        bty = "n",
        cex = 1.5,
        border = c(group_colors[(ncol(heatmap)+1):length(group_colors)],group_colors[1:ncol(heatmap)]),
        text.width = strwidth(colnames(heatmap))[1]*2,
        title.adj = 0.05,
        inset = -0.01,
        xpd = TRUE
    )
    if(!is.null(heatmap_factor)) {
        legend(
            # x = 1.25,y = -0.35,
            x = 1.15,y = -0.4,
            legend = c("absent","present"),
            fill = c("grey","black"),
            title = "",
            xjust = 0,
            bty = "n",
            cex = 1.5,
            border = c("grey","black"),
            text.width = strwidth(colnames(heatmap))[1]*2,
            title.adj = 0.05,
            inset = -0.01,
            xpd = TRUE
        )  
    }
    if(!is.null(heatmap)) {
        draw(
            cont_legend,
            x = unit(0.865,"npc"), # x = unit(0.9,"npc")
            y = unit(0.7,"npc") # y = unit(0.6,"npc")
        )
    }
    dev.off()
  
}
