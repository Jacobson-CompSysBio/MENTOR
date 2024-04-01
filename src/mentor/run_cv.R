########################################################################
# Perform K-fold Cross Validation on a gene set using RWR to find the RWR rank of the left-out genes
# - Input: Pre-computed multiplex network and a geneset
# - Output: Table with the ranking of each gene in the gene set when left out, along with AUPRC and AUROC curves
# Author: David Kainer
########################################################################

suppressWarnings(suppressPackageStartupMessages(require(igraph)))
suppressWarnings(suppressPackageStartupMessages(require(tidyverse)))
suppressWarnings(suppressPackageStartupMessages(require(data.table)))
suppressWarnings(suppressPackageStartupMessages(require(BiocManager)))
suppressWarnings(suppressPackageStartupMessages(require(doParallel)))
suppressWarnings(suppressPackageStartupMessages(require(foreach)))
suppressWarnings(suppressPackageStartupMessages(require(optparse)))
suppressWarnings(suppressPackageStartupMessages(require(patchwork)))
suppressWarnings(suppressPackageStartupMessages(require(iterators)))
suppressWarnings(suppressPackageStartupMessages(require(RandomWalkRestartMH)))
suppressWarnings(suppressPackageStartupMessages(require(IRkernel)))
suppressWarnings(suppressPackageStartupMessages(require(stringr)))
suppressWarnings(suppressPackageStartupMessages(require(Matrix)))
suppressWarnings(suppressPackageStartupMessages(require(dplyr)))

parse_arguments <- function() {
  
  suppressPackageStartupMessages(require(optparse))
  option_list <- list(
    make_option(c("-d","--data"),
      action = "store",
      default = NULL,
      type = "character",
      help = "The path to the .Rdata file for your combo of underlying
              functional networks. This file is produced by RWR_make_multiplex."
    ),
    make_option(c("-g","--geneset"),
      action = "store",
      default = NULL,
      type = "character",
      help = "The path to the gene set file. It must have the following
                    first two columns with no headers tab-delimited:
                    <setid> <gene> <weight>."
    ),
    make_option(c("-o","--outdir"),
      action = "store",
      default = NULL,
      type = "character",
      help = "Path to the output directory. Both 'fullranks' and
              'meanranks' will be saved here with auto-generated filenames.
              (--out-fullranks and --out-meanranks override this.)"
    ),
    make_option(c("-t","--threads"),
      action = "store",
      default = parallel::detectCores() - 1,
      type = "numeric",
      help = "Number of threads to use. Default for your system is
                    all cores - 1. [default %default]"
    ),
    make_option(c("-v","--verbose"),
      action = "store_true",
      default = FALSE,
      help = "Verbose mode. [default %default]"
    )
  )
  opt <- parse_args(OptionParser(option_list = option_list),convert_hyphens_to_underscores = TRUE)
  errors <- 0
  if (is.null(opt$data)) {
    message("ERROR:: --data is required.")
    errors <- errors + 1
  }
  if (is.null(opt$geneset)) {
    message("ERROR:: --geneset is required.")
    errors <- errors + 1
  }
  if (opt$verbose) {
    #print(opt)
  }
  if (errors > 0) {
    quit()
  }
  return(opt)
  
}


main <- function(opt) {

  initial.options <- commandArgs(trailingOnly = FALSE)
  file.arg.name <- "--file="
  script.name <- sub(file.arg.name,"",initial.options[grep(file.arg.name,initial.options)])
  script.basename <- dirname(script.name)
  utils_path = file.path(script.basename,'rwr_utils.R')
  rwr_path = file.path(script.basename,'rwr_cv.R')
  source(rwr_path)
  if (file.exists(utils_path)) {
    source(utils_path)
  } else {
    message(sprintf('ERROR: Cannot find utils.R at path: %s',utils_path))
    return(1)
  }
  opt <- parse_arguments()
  RWR_CV(
    data = opt$data,
    geneset_path = opt$geneset,
    method = "singletons",
    folds = 5,
    restart = 0.7,
    tau = 1.0,
    numranked = 1.0,
    outdir = opt$outdir,
    threads = opt$threads,
    verbose = opt$verbose,
    write_to_file = TRUE
  )
  return(0)
  
}

status <- main()
quit(save = "no",status = status)
