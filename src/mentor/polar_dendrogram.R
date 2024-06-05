
# function to create polar dendrogram
circlize <- function(dend_list,heatmap,heatmap_factor,clusters,cluster_label_size,split,col_fun1,labels_size,group_colors,max_height,track_height,highlight_index,highlight_color,height,width) {
  
    # clear circos
    circos.clear()
    # setting up global parameters for circos plot
    circos.par(
        #clock.wise = FALSE, # uncomment to make sectors appear in anti-clockwise direction
        #start.degree = 219.5, # uncomment to change the start degree of the first cluster
        "canvas.xlim" = c(-1,1),
        "canvas.ylim" = c(-1,1),
        gap.degree = c(rep(5,nrow(clusters) - 1),8),
        points.overflow.warning = FALSE
    )
    # set up track heights & initialize
    if(is.null(heatmap) & is.null(heatmap_factor)) {
        circos.initialize(split,xlim = c(0,1))
        dend_th <- track_height
    } else if(!is.null(heatmap) & is.null(heatmap_factor)) {
        circos.heatmap.initialize(heatmap,split = split,cluster = FALSE)
        heatmap_th <- track_height[1]
        dend_th <- track_height[2]
    } else if(is.null(heatmap) & !is.null(heatmap_factor)) {
        circos.heatmap.initialize(heatmap_factor,split = split,cluster = FALSE)
        heatmap_factor_th <- track_height[1]
        dend_th <- track_height[2]
    } else {
        circos.heatmap.initialize(heatmap,split = split,cluster = FALSE)
        heatmap_th <- track_height[1]
        heatmap_factor_th <- track_height[2]
        dend_th <- track_height[3]
    }
    # adding cluster labels
    circos.track(
        ylim = c(0,1),
        panel.fun = function(x,y) {
            circos.text(
                CELL_META$xcenter,
                CELL_META$cell.ylim[2],
                paste0("C",CELL_META$sector.index),
                facing = "inside",
                cex = cluster_label_size,
                col = clusters[as.numeric(CELL_META$sector.index) + 1,]$color,
                adj = c(0.5,0),
                niceFacing = TRUE
            )
        },bg.border = NA,track.height = 0.005
    )
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
            track.height = heatmap_th
        )
        # add color blocks for continuous heatmap (group labels)
        circos.track(
            track.index = get.current.track.index(),
            panel.fun = function(x,y) {
                if(CELL_META$sector.numeric.index == nrow(clusters)) {
                    for(i in 1:ncol(heatmap)) {
                        circos.rect(
                            CELL_META$cell.xlim[2] + convert_x(5.5,"mm"),i - 1, # might have to make 5.5 and 15.5 dynamic
                            CELL_META$cell.xlim[2] + convert_x(15.5,"mm"),i,
                            col = group_colors[i],
                            border = NA
                        )
                    }
                }
            },bg.border = NA,track.height = heatmap_th
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
            track.height = heatmap_factor_th
        )
        # add color blocks for factor heatmap (group labels)
        circos.track(
            track.index = get.current.track.index(),
            panel.fun = function(x,y) {
                if(CELL_META$sector.numeric.index == nrow(clusters)) {
                    for(i in 1:ncol(heatmap_factor)) {
                        circos.rect(
                            CELL_META$cell.xlim[2] + convert_x(4,"mm"),i - 1, # might have to make 4 and 11 dynamic
                            CELL_META$cell.xlim[2] + convert_x(11,"mm"),i,
                            col = group_colors[i],
                            border = NA
                        )
                    }
                }
            },bg.border = NA,track.height = heatmap_factor_th
        )
    }
    # add dendrogram track to circos plot
    circos.track(
        ylim = c(0,max_height),
        panel.fun = function(x,y) {
            sector.index = get.cell.meta.data("sector.index")
            dend = dend_list[[as.numeric(sector.index) + 1]]
            circos.dendrogram(
                dend,
                max_height = max_height,
                facing = "outside"
            )
        },bg.border = NA,track.height = dend_th
    )
    # check if the index for highlighting is not null
    if(!is.null(highlight_index)) {
        if(!is.null(heatmap) & !is.null(heatmap_factor)) {
            total_tracks <- 6
        } else if(is.null(heatmap)|is.null(heatmap_factor)) {
            total_tracks <- 4
        } else {
            total_tracks <- 2
        }
        highlight_index <- do.call("c",strsplit(highlight_index,","))
        # draw the highlight for the user input sector
        lapply(unique(highlight_index),function(x){
            draw.sector(
                get.cell.meta.data("cell.start.degree",sector.index = x) + 2.5,
                get.cell.meta.data("cell.end.degree",sector.index = x) - 2.5,
                rou1 = get.cell.meta.data("cell.top.radius",track.index = 1) + 0.05,
                rou2 = get.cell.meta.data("cell.bottom.radius",track.index = total_tracks),
                col = adjustcolor(highlight_color,0.20),
                border = NA
            )
        })
    }
    # if(!is.null(zoom)) {
    #     placeholder for adding a zoom option to visualize a single section
    # }
    circos.clear()
  
}

# function to prepare data for polar dendrogram
polar_dendrogram <- function(dend_labs,dend2,heatmap,squish_bounds,cluster_label_size,labels_size,group_colors,track_height,highlight_index,highlight_color,legend_title,plot_file,height,width) {
  
    # create data frame with dendrogram labels and heatmap
    heatmap <- merge(dend_labs,heatmap,all.y = TRUE,by = "label")
    # order by the original row order
    heatmap <- heatmap[order(heatmap$x,heatmap$y,decreasing = FALSE),]
    # get nuique clusters
    clusters <- unique(heatmap[,c("col.x","cluster")])
    # change column names
    names(clusters) <- c("color","cluster")
    # create split for circos plot
    split <- unique(heatmap[,c("label","cluster")])$cluster
    # get heatmap columns used for circos plot
    heatmap <- heatmap[,c("label","value","source","cluster","row_order")]
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
        if(any(factor_cols)) {
            if(is.null(group_colors)) {
                group_colors <- primary.colors(n = ncol(heatmap) + ncol(heatmap_factor) + 1,no.white = TRUE)
                group_colors <- group_colors[which(group_colors != "#000000")]
                legend_cols <- c(rev(colnames(heatmap_factor)),rev(colnames(heatmap)))
                fill_cols <- c(group_colors[(ncol(heatmap)+1):length(group_colors)],group_colors[1:ncol(heatmap)])
            }
        } else {
            if(is.null(group_colors)) {
                group_colors <- primary.colors(n = ncol(heatmap) + 1,no.white = TRUE)
                group_colors <- group_colors[which(group_colors != "#000000")]
                legend_cols <- rev(colnames(heatmap))
                fill_cols <- group_colors
            }
        }
    } else {
        if(is.null(group_colors)) {
            group_colors <- primary.colors(n = ncol(heatmap_factor) + 1,no.white = TRUE) 
            group_colors <- group_colors[which(group_colors != "#000000")]
            legend_cols <- rev(colnames(heatmap_factor))
            fill_cols <- group_colors
        }
    }
    # color scale for continuous columns
    if(!is.null(squish_bounds)) {
        squish_bounds <- as.numeric(do.call("c",strsplit(squish_bounds,",")))
        col_fun1 <- colorRamp2(c(squish_bounds[1],0,squish_bounds[2]),c("blue","white","red"))  
    } else {
        col_fun1 <- colorRamp2(c(-2,0,2),c("blue","white","red"))
    }
    # create list of dendrograms based on our sub-/clustering
    dend_list <- lapply(unique(dend_labs$cluster),function(x) {
      dend_labs[dend_labs$cluster == x,]$ensembl
      #dend_labs[dend_labs$cluster == x,]$label
    })
    dend_list <- lapply(dend_list,function(x) {
      find_dendrogram(dend2,selected_labels = x)
    })
    class(dend_list) <- "dendlist"
    if(!any(factor_cols)) heatmap_factor <- NULL
    if(all(factor_cols)) heatmap <- NULL
    max_height <- max(sapply(dend_list,function(x) attr(x,"height")))
    track_height <- as.numeric(do.call("c",strsplit(track_height,",")))
    par(mai=c(0.005,0.005,0.005,0.005))
    pdf(plot_file,height = height,width = width)
    circlize(
        dend_list = dend_list,
        heatmap = heatmap,
        heatmap_factor = heatmap_factor,
        clusters = clusters,
        cluster_label_size = cluster_label_size,
        split = split,
        col_fun1 = col_fun1,
        labels_size = labels_size,
        group_colors = group_colors,
        max_height = max_height,
        track_height = track_height,
        highlight_index = highlight_index,
        highlight_color = highlight_color,
        height = height,
        width = width
    )
    if(!is.null(heatmap)) {
        #cont_legend <- Legend(
        #    at = c(squish_bounds[1],0,squish_bounds[2]),
        #    col_fun = col_fun1,
        #    title_position = "topleft",
        #    title = TeX("$\\log_{2}(FC)$"), # fix this legend_title[1]
        #    title_gp = gpar(fontface = 1,cex = 1.5),
        #    labels_gp = gpar(fontface = 1,cex = 1.25),
        #    title_gap = unit(10,"mm"),
        #    grid_height = unit(2,"inches"),
        #    grid_width = unit(0.25,"inches")
        #)
        # get max label size to add to first cont legend number (empty string)
        max_nchar <- max(nchar(legend_cols))
        lgd <- rep(NA, 11)
	# need to put in a check here for squish bounds
	if(!is.null(squish_bounds)) {
	    if(squish_bounds[1] < 1) {
	        lgd[c(1,6,11)] <- c(squish_bounds[2],0,squish_bounds[1])
	    } else {
	        lgd[c(1,6,11)] <- c(squish_bounds[1],0,squish_bounds[2])
	    }
	} else {
	    lgd[c(1,6,11)] <- c(2,0,-2)
	}
        lgd[1] <- paste0(lgd[1],paste0(rep(" ",max_nchar + 20),collapse = ""))
        cont_lgd <- legend(
            "right",
            legend = lgd,
            fill = colorRampPalette(colors = c('red','white','blue'))(11),
            title = "",
            xjust = 0,
            bty = "n",
            border = NA,
            cex = 1.25,
            xpd = TRUE,
            y.intersp = 0.5,
            plot = FALSE
        )
        legend(
            "right",
            legend = lgd,
            fill = colorRampPalette(colors = c('red','white','blue'))(11),
            title = "",
            xjust = 0,
            bty = "n",
            border = NA,
            cex = 1.25,
            xpd = TRUE,
            y.intersp = 0.5
        )
        legend(
            x = cont_lgd$rect$left,y = -0.15,
            legend = legend_cols,
            fill = fill_cols,
            title = legend_title[2],
            xjust = 0,
            bty = "n",
            cex = 1.25,
            border = fill_cols,
            title.adj = 0.05,
            xpd = TRUE
        )
    }
    if(!is.null(heatmap_factor)) {
        if(!is.null(heatmap)) {
            legend(
                x = cont_lgd$rect$left,y = 0.25,
                legend = c("absent","present"),
                fill = c("grey","black"),
                title = "",
                xjust = 0,
                bty = "n",
                cex = 1.5,
                border = c("grey","black"),
                title.adj = 0.05,
                xpd = TRUE
            ) 
        } else {
            legend(
                "right",
                legend = c("absent","present"),
                fill = c("grey","black"),
                title = "",
                xjust = 0,
                bty = "n",
                cex = 1.5,
                border = c("grey","black"),
                title.adj = 0.05,
                xpd = TRUE
            )
        }
    }
    #if(!is.null(heatmap)) {
       #draw(
       #    cont_legend,
       #    x = unit(0.875,"npc"), # x = unit(0.9,"npc") 0.868
       #    y = unit(0.68,"npc") # y = unit(0.6,"npc")
       #)
    #}
    whatever <- dev.off()
  
}
