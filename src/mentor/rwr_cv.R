########################################################################
# Perform K-fold Cross Validation on a gene set using RWR to find the RWR rank
# of the left-out genes:
# - Input: Pre-computed multiplex network and a geneset
# - Output: Table with the ranking of each gene in the gene set when left out,
#           along with AUPRC and AUROC curves
########################################################################
# Internal Functions
########################################################################

update_folds_by_method <- function(geneset,method,num_folds,verbose = FALSE) {
  
  chunks <- NULL
  folds <- nrow(geneset)
  cat("\nCross-validation method is ",
    method,
    ", i.e, folds=",
    folds,
    "\n",
    sep = "",
    file = stderr()
  )
  return(list(folds=folds,geneset=geneset,chunks=chunks,method=method))
  
}

extract_lo_and_seed_genes_cv <- function(geneset,method,r,chunks = NULL) {
  
  seed_genes <- geneset %>%
    dplyr::slice(r) %>%
    dplyr::pull(gene) # nolint       Get one seed gene.  
  leftout <- geneset %>%
    dplyr::filter(gene != seed_genes) %>%     #nolint
    dplyr::pull(gene)
  list(leftout, seed_genes)
  
}

create_rankings_cv <- function(rwr,networks,r,name,geneset,method,seed_genes,leftout,num_nodes_in_multiplex) {
  
  # turn scores into rankings: highest score gets lowest rank
  rwr$RWRM_Results <- rwr$RWRM_Results %>%
    dplyr::mutate(
      rank = dplyr::row_number(-Score),
      InValset = as.numeric(rwr$RWRM_Results$NodeNames %in% leftout)
    ) %>%
    dplyr::mutate(rank = dplyr::if_else(Score == 0,
      true = as.integer(num_nodes_in_multiplex),
      false = rank
    ))
  num_in_network <- nrow(geneset)
  num_seeds <- length(seed_genes)
  num_leftout <- length(leftout)
  networks <- networks
  fold <- r
  modname <- name
  geneset <- geneset$setid[1]
  seed <- seed_genes
  leftout <- "AllBut1"
  rwr$RWRM_Results <- rwr$RWRM_Results %>%
    dplyr::mutate(
      num_in_network = num_in_network,
      num_seeds = num_seeds,
      num_leftout = num_leftout,
      networks = networks,
      fold = fold,
      modname = modname,
      geneset = geneset,
      seed = seed,
      leftout = leftout,
      method = method
    )
  rwr$RWRM_Results
  
}

# For 'loo' and 'kfold', each gene is leftout once, so gets ranked once.
RWR <- function(geneset,
                adjnorm,
                mpo,
                method,
                num_folds,
                chunks,
                restart = 0.7,
                tau = c(1, 1),
                name = "default",
                threads = 1,
                verbose = FALSE) {
  
  networks <- paste(names(mpo)[1:mpo$Number_of_Layers], collapse = "_")
  doParallel::registerDoParallel(cores = threads)
  res <- foreach::foreach(r = 1:num_folds) %dopar% {
    lo_seed_genelist <- extract_lo_and_seed_genes_cv(
      geneset,
      method,
      r,
      chunks
    )
    leftout <- lo_seed_genelist[[1]]
    seed_genes <- lo_seed_genelist[[2]]
    ### run RWR on a fold
    rwr <- RandomWalkRestartMH::Random.Walk.Restart.Multiplex(
      x = adjnorm,
      MultiplexObject = mpo,
      Seeds = seed_genes,
      r = restart,
      tau = tau
    )
    ranking_results <- create_rankings_cv(
      rwr,
      networks,
      r,
      name,
      geneset,
      method,
      seed_genes,
      leftout,
      mpo$Number_of_Nodes_Multiplex
    )
    ranking_results
  }
  doParallel::stopImplicitCluster()
  res
  
}

calc_metrics_cv <- function(res_combined,res_avg) {
  
  ### calculate ROC/PRC per-fold
  res_combined <- res_combined %>%
    dplyr::mutate(fold = as.factor(fold)) %>%
    dplyr::group_by(fold) %>%
    dplyr::group_modify(~ calc_ROCPRC(df=.x, scorescol = "rank", labelscol = "InValset"))
  # 1. precision @ numleftout (aka R-PREC)
  output <- res_combined %>%
    dplyr::group_by(fold) %>%
    dplyr::filter(rank==num_leftout) %>%
    dplyr::summarise(value = PREC,measure = "P@NumLeftOut")
  # 2. Avg PRC (uses Average Precision, not interpolated precision)
  output <- rbind(output, res_combined %>%
    dplyr::group_by(fold) %>%
    dplyr::filter(TP==1) %>%
    dplyr::summarise(value = sum(PREC)/sum(InValset),measure = "AvgPrec"))
  # 3. AUPRC (provide both the AUPRC per fold, and the AUPRC expected by an unskilled model)
  output <- rbind(output, res_combined %>%
    dplyr::group_by(fold) %>%
    dplyr::summarise(value = area_under_curve(REC, PREC, method="trapezoid", ties="max"), measure="AUPRC"))
  output <- rbind(output, res_combined %>%
    dplyr::group_by(fold) %>% 
    dplyr::summarise(value = dplyr::first(num_leftout)/dplyr::n(), measure="ExpectedAUPRC"))
  # 4. AUROC
  output <- rbind(output, res_combined %>%
    dplyr::group_by(fold) %>%
    dplyr::summarise(value = sum(REC)/dplyr::n(), measure="AUROC"))
  ### get metrics based on reranking of mean ranks
  res_avg <- calc_ROCPRC(res_avg, scorescol = "rerank",labelscol = "InValset")
  # 5. AvgPrec of mean ranks (uses Average Precision, not interpolated avg precision)
  output <- rbind(output, res_avg %>%
    dplyr::filter(TP==1) %>%
    dplyr::summarise(fold="meanrank",value = sum(PREC)/sum(InValset),measure="AvgPrec"))
  # 6. AUPRC of mean ranks
  output <- rbind(output, res_avg %>%
    dplyr::summarise(fold="meanrank",value = area_under_curve(REC,PREC,method="trapezoid",ties="max"),measure="AUPRC"))
  # 7. AUROC of mean ranks
  output <- rbind(output, res_avg %>%
    dplyr::summarise(fold="meanrank",value = sum(REC)/dplyr::n(),measure="AUROC"))
  output$geneset <- res_combined$geneset[1]
  return(list(summary = output,res_combined = res_combined,res_avg = res_avg))

}

calculate_max_precision <- function(pr, metrics) {
  
  maxprec <- foreach::foreach(
    f = pr$fold@.Data,
    r = pr$recall,
    .combine = c
  ) %do% {
    prec <- metrics$res_combined %>%
      dplyr::filter(fold == f, REC >= r) %>%
      dplyr::pull(PREC)
    if (length(prec) > 0) {
      max(prec)
    } else {
      return(0)
    }
  }
  maxprec
  
}

post_process_rwr_output_cv <- function(res,extras,folds,nw.mpo) {
  
  res_method <- res[[1]]$method[1]
  extras_exist <- !is.null(extras)
  if (extras_exist) {
    for (i in 0:(nrow(extras) - 1))
    {
      res.tmp <- res[[(i %% folds) + 1]] # Get res[ i mod extras ]'th res
      extra_ith_row <- data.frame(
        NodeNames = extras$gene[i + 1],
        Score = 0,
        rank = max(res.tmp$rank) + 1,
        InValset = 1,
        num_in_network = dplyr::first(res.tmp$num_in_network),
        num_seeds = dplyr::first(res.tmp$num_seeds),
        num_leftout = dplyr::first(res.tmp$num_leftout),
        networks = dplyr::first(res.tmp$networks),
        fold = dplyr::first(res.tmp$fold),
        modname = dplyr::first(res.tmp$modname),
        geneset = dplyr::first(res.tmp$geneset),
        seed = "missing",
        leftout = "missing",
        method = dplyr::first(res.tmp$method)
      )
      res[[(i %% folds) + 1]] <- rbind(res.tmp, extra_ith_row)
    }
  }
  res_combined <- dplyr::bind_rows(res) %>% dplyr::arrange(rank)
  res_combined
  
}

calculate_average_rank_across_folds_cv <- function(res_combined){

  res_avg <- res_combined %>%
    dplyr::group_by(NodeNames) %>%
    dplyr::summarise(meanrank = mean(rank), 
      InValset = dplyr::first(InValset), 
      geneset=dplyr::first(geneset),
      num_in_network=dplyr::first(num_in_network)) %>%
    dplyr::arrange(meanrank) %>%
    dplyr::mutate(rerank = rank(meanrank, ties.method = "min"),.after=meanrank) # here we rerank the mean ranks of the folds to go from 1:nrow again
  res_avg
  
}

########################################################################
# Main Function
########################################################################

RWR_CV <- function(data = NULL,
                   geneset_path = NULL,
                   method = "kfold",
                   folds = 5,
                   restart = 0.7,
                   tau = 1.0,
                   numranked = 1.0,
                   outdir = NULL,
                   modname = "default",
                   out_full_ranks = NULL,
                   out_mean_ranks = NULL,
                   threads = 1,
                   verbose = FALSE,
                   write_to_file = FALSE) {

  data_list <- load_multiplex_data(data)
  nw_mpo <- data_list$nw.mpo
  nw_adjnorm <- data_list$nw.adjnorm
  if ((!is.null(outdir)|!is.null(out_full_ranks)|!is.null(out_mean_ranks)) & write_to_file == FALSE) {
    warning(sprintf("write_to_file was set to false, however, an output file path was set. write_to_file has been updated to TRUE.\n")) #nolint warning
    write_to_file <- TRUE
  }
  if (is.null(outdir)){
    outdir <- './'
  }
  tau <- get_or_set_tau(nw_mpo, tau)
  geneset_list <- load_geneset(geneset_path,nw_mpo,verbose = verbose)
  geneset <- geneset_list$geneset
  extras <- geneset_list$extras
  updated_data_list <- update_folds_by_method(geneset,method,folds)
  folds <- updated_data_list$folds
  geneset <- updated_data_list$geneset
  chunks <- updated_data_list$chunks
  method <- updated_data_list$method
  res <- RWR(geneset,nw_adjnorm,nw_mpo,method,folds,chunks,restart,tau,modname,threads,verbose)
  res_combined <- post_process_rwr_output_cv(res,extras,folds,nw_mpo)
  res_avg <- calculate_average_rank_across_folds_cv(res_combined)
  metrics <- calc_metrics_cv(res_combined, res_avg)
  if (!is.null(out_full_ranks)) {
    out_path <- out_full_ranks
  } else {
    out_path <- get_file_path("RWR-CV_",
      res_combined$geneset[1],
      get_base_name(data),
      modname,
      outdir = outdir,
      ext = ".fullranks.tsv"
    )
  }
  if (write_to_file) {
    if (!file.exists(outdir)) {
      print("Creating directory")
      dir.create(outdir,recursive = TRUE)
    }
    combined <- res_combined %>%
      dplyr::group_by(fold) %>%
      dplyr::slice_head(prop = numranked)
    write_table(combined,out_path)
  }
  if (!is.null(out_mean_ranks)) {
    out_path <- out_mean_ranks
  } else {
    out_path <- get_file_path("RWR-CV_",
      res_avg$geneset[1],
      get_base_name(data),
      modname,
      outdir = outdir,
      ext = ".meanranks.tsv"
    )
  }
  if (write_to_file) {
    write_table(
      res_avg %>%
        dplyr::slice_head(prop = numranked),
      out_path
    )
  }
  out_path <- get_file_path("RWR-CV_",
    metrics$res_combined$geneset[1],
    get_base_name(data),
    modname,
    outdir = outdir,
    ext = ".metrics.tsv"
  )
  if (write_to_file) {
    write_table(metrics$res_avg, out_path)
  }
  out_path <- get_file_path("RWR-CV_",
    metrics$res_combined$geneset[1],
    get_base_name(data),
    modname,
    outdir = outdir,
    ext = ".summary.tsv"
  )
  if (write_to_file) {
    write_table(metrics$summary, out_path)
  }
  return(
    list(
      "fullranks" = res_combined,
       "meanranks" = res_avg,
       "metrics" = metrics$res_avg,
       "summary" = metrics$summary
    )
  )
  
}
