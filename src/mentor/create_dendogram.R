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
suppressWarnings(suppressPackageStartupMessages(require(latex2exp)))

#################### Argument parser ##############################

option_list <- list(
  make_option(
    c("-d","--distances"), 
    type = "character",
    default = NULL, 
    help = "full path to the distance matrix file", 
    metavar = "character"
  ),
  make_option(
    c("-k","--clusters"),
    type = "integer",
    default = 3,
    help = "the number of clusters desired in the dendrogram", 
    metavar = "character"
  ),
  make_option(
    c("-x","--map"), 
    type = "character",
    default = NULL, 
    help = "full path to the gene symbol mapping file", 
    metavar = "character"
  ),
  make_option(
    c("-o","--outdir"), 
    type = "character",
    default = NULL, 
    help = "full path to output the dendrogram", 
    metavar = "character"
  ),
  make_option(
    c("-f","--outfile"),
    type = "character",
    default = NULL,
    help = "filename to append",
    metavar = "character"
  ),
  make_option(
    c("-s","--subcluster"), 
    action = "store_true",
    default = FALSE, 
    help = "specify if want to subcluster", 
    metavar = "character"
  ),
  make_option(
    c("-i","--increment"), 
    type = "integer",
    default = 5, 
    help = "increment to use in subclustering", 
    metavar = "character"
  ),
  make_option(
    c("-m","--maxsize"), 
    type = "integer",
    default = 40, 
    help = "maximum size for clades if subclustering", 
    metavar = "character"
  ),
  make_option(
    c("-z","--heatmaps"), 
    type = "character",
    default = NULL, 
    help = "full path to the heatmap file if you want to include a heatmap", 
    metavar = "character"
  ),
  make_option(
    c("-l","--reordercols"),
    action = "store_true",
    default = FALSE,
    help = "specify if want to subcluster",
    metavar = "character"
  ),
  make_option(
    c("-a","--legend"),
    type = "character",
    default = NULL,
    help = "title to give to continuous legend",
    metavar = "character"
  ),
  make_option(
    c("-q","--squish"),
    type = "character",
    default = NULL,
    help = "if you have FC values squish the upper and lower bounds to two numbers formatted as 'lower_bound,upper_bound'",
    metavar = "character"
  ),
  make_option(
    c("-r","--relwidths"),
    type = "character",
    default = "1,1",
    help = "if you have heatmap then provide two numbers 'x,y' for the relative widths in the combined plot",
    metavar = "character"
  ),
  make_option(
    c("-e","--plotwidth"),
    type = "integer",
    default = 30,
    help = "width of the final plot",
    metavar = "character"
  ),
  make_option(
    c("d","--plotheight"),
    type = "integer",
    default = NULL,
    help = "height of the final plot",
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
  reordercols,
  legend,
  squish_bounds,
  relative_widths,
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
    cluster_file <- paste0("clusters_k",k,".tsv") 
    # create plot filename
    plot_file <- paste0("dendrogram_k",k,".pdf")
    # source the dendrogram.R file
    dendrogram.path <- file.path(script.basename, "/dendrogram.R")
    source(dendrogram.path)
    # create dendrogram
    dendo <- dendrogram(dis_mat = dm_dist,k,map)
    
  } else {

    cat("\n\nsubclustering")
    # create cluster file
    cluster_file <- paste0("clusters_sk",k,"_ki",k_increment,"_ms",max_size,".tsv")
    # create plot filename
    plot_file <- paste0("subcluster_dendrogram_sk",k,"_ki",k_increment,"_ms",max_size,".pdf")
    # source the subclustered_dendrogram.R file
    subclustered.path <- file.path(script.basename, "/subclustered_dendrogram.R")
    source(subclustered.path)
    # create sub-clustered dendrogram
    dendo <- subclustered_dendrogram(dis_mat = dm_dist,k,k_increment,max_size,map)
    
  }
  
  # grab the dendrogram and labels from result
  dendrogram <- dendo$dendrogram
  dend_labs <- dendo$dendrogram_labels
 
  if(!is.null(heatmaps)) {
    
    cat("\n\nadding heatmap")
    # adjust plot filename
    plot_file <- gsub(".pdf","_logfc_heatmap.pdf",plot_file)
    # source the heatmaps.R file
    heatmaps.path <- file.path(script.basename, "/heatmaps.R")
    source(heatmaps.path)
    # create the heatmap to add to dendrogram
    heat <- heatmap(heatmap = heatmaps,dend_labs,reordercols,legend,squish_bounds)
    # grab relative widths for final plot
    relative_widths <- do.call("c",strsplit(relative_widths,","))
    # add to dendrogram
    dendrogram <- plot_grid(
      dendrogram,
      heat,
      nrow = 1,
      align = "h",
      axis = "l",
      rel_widths = as.numeric(relative_widths)
    )
    
  }
  
  # export the clusters and dendrogram
  if(!is.null(out_file)) {
    cluster_file <- paste0(out_file,"_",cluster_file)
    plot_file <- paste0(out_file,"_",plot_file)
  }
  dend_labs$row_order <- 1:nrow(dend_labs)
  groups <- data.frame(
    "col" = unique(dend_labs$col),
    "cluster" = 0:(length(unique(dend_labs$col))-1)
  )
  dend_labs <- merge(dend_labs,groups,by = "col",all.x = TRUE)
  dend_labs <- dend_labs[order(dend_labs$row_order,decreasing = FALSE),]
  dend_labs <- dend_labs[,c("label","cluster")]
  cat("\n\nexporting clusters")
  write.table(dend_labs,paste0(out_dir,cluster_file),sep = "\t",col.names = TRUE,row.names = FALSE,quote = FALSE)
  # adjust heigh of plot based on user input or by the number of genes
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

#################### Parse & apply ##############################

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
  reordercols = opt$reordercols,
  legend = opt$legend,
  squish_bounds = opt$squish,
  relative_widths = opt$relwidths,
  plot_width = opt$plotwidth,
  plot_height = opt$plotheight
  
)

