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
# suppressPackageStartupMessages(require(grid))
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
#  make_option(
#    c("-e","--export"), 
#    action = "store_true",
#    default = FALSE, 
#    help = "export the plot", 
#    metavar = "character"
#  ),
  make_option(
    c("-z","--heatmaps"), 
    type = "character",
    default = NULL, 
    help = "full path to the heatmap file if you want to include a heatmap", 
    metavar = "character"
  ),
  make_option(
    c("-p","--pcutoff"),
    type = "double",
    default = NULL,
    help = "adjusted p-value cutoff to use if p-values present in the heatmap table",
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
    c("-w","--plotwidth"),
    type = "integer",
    default = 30,
    help = "width of the final plot",
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
  subcluster,
  k_increment,
  max_size,
  #export,
  heatmaps,
  p_cutoff,
  squish_bounds,
  relative_widths,
  plot_width
  
) {
  
  # import dissimilarity_matrix as dm
  dm <- suppressMessages(read_tsv(dis_mat, col_names = TRUE, show_col_types = FALSE))
  # make first column row names
  dm_matrix <- dm  %>%
    column_to_rownames(var = "...1")  %>%
    as.matrix()
  # convert to distance matrix
  dm_dist <- dm_matrix  %>% as.dist()
  
  # get current directory of path to source other files
  initial.options <- commandArgs(trailingOnly = FALSE)
  file.arg.name <- "--file="
  script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
  script.basename <- dirname(script.name)
  
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
    
    # adjust plot filename
    plot_file <- gsub(".pdf","_logfc_heatmap.pdf",plot_file)
    # source the heatmaps.R file
    heatmaps.path <- file.path(script.basename, "/heatmaps.R")
    source(heatmaps.path)
    # create the heatmap to add to dendrogram
    heat <- heatmap(heatmap = heatmaps,dend_labs,p_cutoff,squish_bounds)
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
  
  #if(export) {
    
    # export the clusters
    dend_labs$row_order <- 1:nrow(dend_labs)
    groups <- data.frame(
      "col" = unique(dend_labs$col),
      "cluster" = 0:(length(unique(dend_labs$col))-1)
    )
    dend_labs <- merge(dend_labs,groups,by = "col",all.x = TRUE)
    dend_labs <- dend_labs[order(dend_labs$row_order,decreasing = FALSE),]
    dend_labs <- dend_labs[,c("label","cluster")]
    write.table(dend_labs,paste0(out_dir,cluster_file),sep = "\t",col.names = TRUE,row.names = FALSE,quote = FALSE)
    # export the ggplot dendrogram
    ggsave(
      paste0(out_dir,plot_file),
      plot = dendrogram,
      width = plot_width,
      height = nrow(dend_labs) * 0.6,
      units = "cm",
      limitsize = FALSE
    )
    
  #} else {
    
  #  print("Warning: to save dendrogram be sure to use --export in your command")  
  #  return(dendrogram)
    
  #}
  
}

#################### Parse & apply ##############################

create_dendogram(
    
  dis_mat = opt$distances,
  k = opt$clusters,
  map = opt$map,
  out_dir = opt$outdir,
  subcluster = opt$subcluster,
  k_increment = opt$increment,
  max_size = opt$maxsize,
  #export = opt$export,
  heatmaps = opt$heatmaps,
  p_cutoff = opt$pcutoff,
  squish_bounds = opt$squish,
  relative_widths = opt$relwidths,
  plot_width = opt$plotwidth
  
)

