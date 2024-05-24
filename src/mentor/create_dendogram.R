# Create dendograms from a dissimilarity matrix

# Authors:
  # Mirko Pavicic
  # Alice Townsend
  # Kyle Sullivan
  # Manesh Shah
  # Erica T. Prates

# Created: 2023-08-21
# Last Modified: 2023-10-08

# Version 0.2 (2023-08-25)
# - Added log2fc heatmap
# Version 0.3 (2023-08-29)
# - Added subclustering function

#################### Load libraries ##############################

suppressWarnings(suppressPackageStartupMessages(require(tidyverse)))
suppressWarnings(suppressPackageStartupMessages(require(RColorBrewer)))
suppressWarnings(suppressPackageStartupMessages(require(dendextend)))
suppressWarnings(suppressPackageStartupMessages(require(cowplot)))
suppressWarnings(suppressPackageStartupMessages(require(optparse)))
suppressWarnings(suppressPackageStartupMessages(require(ggnewscale)))
suppressWarnings(suppressPackageStartupMessages(require(reshape2)))
suppressWarnings(suppressPackageStartupMessages(require(circlize)))
suppressWarnings(suppressPackageStartupMessages(require(ComplexHeatmap)))
suppressWarnings(suppressPackageStartupMessages(require(colorRamps)))
suppressWarnings(suppressPackageStartupMessages(require(gridBase)))
suppressWarnings(suppressPackageStartupMessages(require(grid)))
suppressWarnings(suppressPackageStartupMessages(require(latex2exp)))
suppressWarnings(suppressPackageStartupMessages(require(readr)))
# library(colorRamps)

#################### Argument parser ##############################

option_list <- list(
  make_option(
    # c("-d","--distances"), 
    "--distances",
    type = "character",
    default = NULL, 
    help = "Full path to the distance matrix file", 
    metavar = "character"
  ),
  make_option(
    # c("-k","--clusters"),
    "--clusters",
    type = "integer",
    default = 3,
    help = "Number of initial clusters desired in dendrogram (initial # if using --subcluster)", 
    metavar = "character"
  ),
  make_option(
    # c("-x","--map"), 
    "--map",
    type = "character",
    default = NULL, 
    help = "Full path to the gene symbol mapping file", 
    metavar = "character"
  ),
  make_option(
    # c("-o","--outdir"), 
    "--outdir",
    type = "character",
    default = NULL, 
    help = "Location to save dendrogram and clusters", 
    metavar = "character"
  ),
  make_option(
    # c("-f","--outfile"),
    "--outfile",
    type = "character",
    default = NULL,
    help = "Description to append to filenames",
    metavar = "character"
  ),
  make_option(
    # c("-s","--subcluster"), 
    "--subcluster",
    action = "store_true",
    default = FALSE, 
    help = "Specify if want to subcluster dendrogram", 
    metavar = "character"
  ),
  make_option(
    # c("-i","--increment"), 
    "--increment",
    type = "integer",
    default = 5, 
    help = "Increment to use in subclustering", 
    metavar = "character"
  ),
  make_option(
    # c("-m","--maxsize"), 
    "--maxsize",
    type = "integer",
    default = 40, 
    help = "Maximum size for clades if subclustering", 
    metavar = "character"
  ),
  make_option(
    # c("-z","--heatmaps"), 
    "--heatmaps",
    type = "character",
    default = NULL, 
    help = "Full path to the heatmap file if you want to include a heatmap", 
    metavar = "character"
  ),
  make_option(
    # c("-z","--heatmaps"), 
    "--plottype",
    type = "character",
    default = "rectangular", 
    help = "Type of dendrogram (either rectangular or polar)", 
    metavar = "character"
  ),
  make_option(
    # c("-l","--reordercols"),
    "--reordercols",
    action = "store_true",
    default = FALSE,
    help = "Specify if you want to perform hierarchical clustering on the columns of the heatmap",
    metavar = "character"
  ),
  make_option(
    # c("-p","--legendtitle"),
    "--legendtitle",
    type = "character",
    default = "Value,Group",
    help = "Titles to give to legends; for rectangular dendrogram title to give to continuous legend; for polar dendrogram titles to give to continuous legend and factor legend",
    metavar = "character"
  ),
  make_option(
    # c("-q","--squish"),
    "--squish",
    type = "character",
    default = NULL,
    help = "If continuous columns in heatmap present squish the upper and lower bounds of color scheme to lower_bound,upper_bound",
    metavar = "character"
  ),
  make_option(
    # c("-r","--relwidths"),
    "--relwidths",
    type = "character",
    default = "1,1",
    help = "If including a heatmap then adjust the relative widths of the dendrogram to the heatmap with dend_width,heatmap_width; only used for rectangular dendrogram",
    metavar = "character"
  ),
  make_option(
    # c("-e","--plotwidth"),
    "--clusterlabelsize",
    type = "double",
    default = 2.5,
    help = "Size of cluster labels associated with polar dendrogram",
    metavar = "character"
  ),
  make_option(
    # c("-e","--plotwidth"),
    "--heatmaplabelsize",
    type = "double",
    default = 0.75,
    help = "Size of the labels associated with the polar dendrogram heatmap (gene names)",
    metavar = "character"
  ),
  make_option(
    # c("-e","--plotwidth"),
    "--groupcolors",
    type = "character",
    default = NULL,
    help = "Colors for the group blocks in polar dendrogram",
    metavar = "character"
  ),
  make_option(
    # c("-r","--relwidths"),
    "--trackheight",
    type = "character",
    default = "0.2,0.2,0.2",
    help = "Width of the tracks in the polar dendrogram (heatmap1,heatmap2,dendrogram)",
    metavar = "character"
  ),
  make_option(
    # c("-e","--plotwidth"),
    "--highlightindex",
    type = "character",
    default = NULL,
    help = "Sector(s) to highlight on polar dendrogram (sector1,sector2,...,sectorN)",
    metavar = "character"
  ),
  make_option(
    # c("-e","--plotwidth"),
    "--highlightcolor",
    type = "character",
    default = "#34EBDC",
    help = "Color to use for highlighted sectors on polar dendrogram",
    metavar = "character"
  ),
  make_option(
    # c("-e","--plotwidth"),
    "--plotwidth",
    type = "integer",
    default = 30,
    help = "Width of the plot",
    metavar = "character"
  ),
  make_option(
    # c("d","--plotheight"),
    "--plotheight",
    type = "integer",
    default = NULL,
    help = "Height of the plot",
    metavar = "character"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser,positional_arguments = FALSE)

#################### Create dendogram ##############################


create_dendogram <- function(
    
  dis_mat,
  k,
  map,
  out_dir,
  out_file,
  subcluster,
  k_increment,
  max_size,
  heatmaps,
  plot_type,
  reordercols,
  legendtitle,
  squish_bounds,
  relative_widths,
  cluster_label_size,
  labels_size,
  group_colors,
  track_height,
  highlight_index,
  highlight_color,
  plot_width,
  plot_height
  
) {
 
  cat("\n\nreading dissimilarity matrix")     	
  # import dissimilarity_matrix as dm
  dm <- suppressMessages(read_tsv(dis_mat, col_names = TRUE, show_col_types = FALSE))
  # make first column row names
  dm_matrix <- dm  %>%
    column_to_rownames(var = "...1")  %>%
    as.matrix()
  # check if k > number of genes
  if(k > ncol(dm_matrix)) {
    cat("\n\nerror: number of clusters specified (k) is greater than the number of genes in geneset")
    quit(save = "no")
  }
  # convert to distance matrix
  dm_dist <- dm_matrix  %>% as.dist()
  
  # get current directory of path to source other files
  initial.options <- commandArgs(trailingOnly = FALSE)
  file.arg.name <- "--file="
  script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
  script.basename <- dirname(script.name)
  cat("\n\nplotting dendrogram")
  
  if(!subcluster) {
   
    # create cluster file
    cluster_file <- paste0(plot_type,"_clusters_k",k,".tsv") 
    # create plot filename
    plot_file <- paste0(plot_type,"_dendrogram_k",k,".pdf")
    # source the dendrogram.R file
    dendrogram.path <- file.path(script.basename, "/dendrogram.R")
    source(dendrogram.path)
    # create dendrogram
    dendo <- dendrogram(dis_mat = dm_dist,k,map,plot_type)
    
  } else {

    cat("\n\nsubclustering")
    # create cluster file
    cluster_file <- paste0(plot_type,"_clusters_sk",k,"_ki",k_increment,"_ms",max_size,".tsv")
    # create plot filename
    plot_file <- paste0(plot_type,"_subcluster_dendrogram_sk",k,"_ki",k_increment,"_ms",max_size,".pdf")
    # source the subclustered_dendrogram.R file
    subclustered.path <- file.path(script.basename, "/subclustered_dendrogram.R")
    source(subclustered.path)
    # create sub-clustered dendrogram
    dendo <- subclustered_dendrogram(dis_mat = dm_dist,k,k_increment,max_size,map,plot_type)
    
  }
  
  # grab the dendrogram and labels from result
  dendrogram <- dendo$plot
  dend2 <- dendo$dendrogram
  dend_labs <- dendo$dendrogram_labels
  
  # get legend title collapse into vector
  legendtitle <- do.call("c",strsplit(legendtitle,","))
 
  if(!is.null(heatmaps)) {
    
    cat("\n\nadding heatmap")
    # adjust plot filename
    plot_file <- gsub(".pdf","_heatmap.pdf",plot_file)
    # source the heatmaps.R file
    heatmaps.path <- file.path(script.basename, "/heatmaps.R")
    source(heatmaps.path)
    # create the heatmap to add to dendrogram
    heat <- heatmap_(plot_type = plot_type,heatmap = heatmaps,dend_labs,reordercols,legendtitle,squish_bounds)
    if(plot_type == "rectangular") {
      # grab relative widths for final plot
      relative_widths <- do.call("c",strsplit(relative_widths,","))
      # add to dendrogram
      dendrogram <- plot_grid(
        dendrogram,
        heat$heat,
        nrow = 1,
        align = "h",
        axis = "l",
        rel_widths = as.numeric(relative_widths)
      )  
    }
    
  }
  
  # export the clusters and dendrogram
  if(!is.null(out_file)) {
    cluster_file <- paste0(out_file,"_",cluster_file)
    plot_file <- paste0(out_file,"_",plot_file)    
  }
  dend_labs$row_order <- 1:nrow(dend_labs)
  groups <- data.frame(
    "col" = unique(dend_labs$col),
    "cluster" = 0:(length(unique(dend_labs$col)) - 1)
  )
  dend_labs <- merge(dend_labs,groups,by = "col",all.x = TRUE)
  dend_labs <- dend_labs[order(dend_labs$row_order,decreasing = FALSE),]
  dend_labs <- dend_labs[,c("ensembl","label","cluster","col")] # ,"x","y","row_order"
  cat("\n\nexporting clusters")
  write.table(
    dend_labs,
    paste0(out_dir,cluster_file),
    sep = "\t",
    col.names = TRUE,
    row.names = FALSE,
    quote = FALSE
  )
  
  if(plot_type == "polar") {
    
    polar_dendrogram.path <- file.path(script.basename, "/polar_dendrogram.R")
    source(polar_dendrogram.path)
    if(is.null(plot_height)) {
      plot_height = 20
    }
    circos.clear()
    polar_dendrogram(
      dend_labs,
      dend2,
      heatmap = heat$heat_labs,
      squish_bounds,
      cluster_label_size,
      labels_size,
      group_colors,
      track_height,
      highlight_index,
      highlight_color,
      legend_title = legendtitle,
      plot_file = paste0(out_dir,plot_file),
      height = plot_height,
      width = plot_width
    )
    
  } else {
    
    # adjust height of plot based on user input or by the number of genes
    if(is.null(plot_height)) {
      plot_height = nrow(dend_labs) * 0.6
    }
    # export the ggplot dendrogram
    cat("\n\nsaving visualization")
    ggsave(
      paste0(out_dir,plot_file),
      plot = dendrogram,
      width = plot_width,
      height = plot_height,
      units = "cm",
      limitsize = FALSE
    )
    
  }
  
}

#################### Parse & apply ##############################


# need to add
# cluster_label_size,rownames_size,track_height,highlight_index,highlight_color
create_dendogram(
    
  dis_mat = opt$distances,
  k = opt$clusters,
  map = opt$map,
  out_dir = opt$outdir,
  out_file = opt$outfile,
  subcluster = opt$subcluster,
  k_increment = opt$increment,
  max_size = opt$maxsize,
  heatmaps = opt$heatmaps,
  plot_type = opt$plottype,
  reordercols = opt$reordercols,
  legendtitle = opt$legendtitle, # change to legend titles (3 elements for polar 1 for rectangular)
  squish_bounds = opt$squish,
  relative_widths = opt$relwidths, # only for rectangular dendrogram 
  cluster_label_size = opt$clusterlabelsize,
  labels_size = opt$heatmaplabelsize,
  group_colors = opt$groupcolors,
  track_height = opt$trackheight,
  highlight_index = opt$highlightindex,
  highlight_color = opt$highlightcolor,
  plot_width = opt$plotwidth,
  plot_height = opt$plotheight
  
)

