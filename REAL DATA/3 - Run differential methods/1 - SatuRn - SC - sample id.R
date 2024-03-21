alias R='R_LIBS=/home/stiber/R/x86_64-pc-linux-gnu-library/4.3 /home/Shared/simone/DifferentialRegulation/software/R-4.3.0/bin/R'

cd /home/Shared_sherborne/simone/Diff_Velo/discovery

R
rm(list = ls())

TIMES = system.time({
  library(satuRn)
  library(SingleCellExperiment)
  library(BiocParallel)
  
  DF_SatuRn = list()
  
  # load USA counts with pre-computed cell-type information:
  name = file.path("sce_cell_type.rds")
  sce = readRDS(name)
  
  # define list of all possible group separations:
  sample_ids <-  paste0("organoid",c(1:3, 16:18))
  group = ifelse(sample_ids %in% paste0("organoid",1:3), "3mon", "6mon")
  
  clusters = unique(sce$cell_types)
  clusters
  
  sort(table(sce$cell_types))
  
  min_counts_per_group = 10
  
  for(cl in seq_along(clusters)){
    cluster = clusters[cl]
    
    # filter cell types:
    sce_one_cluster = sce[,sce$cell_types == cluster ]
    
    n_samples = length(sample_ids);
    n_samples
    
    SUA = list()
    for(i in 1:n_samples){
      sel_cells = sce_one_cluster$sample_id == sample_ids[i]
      sce_one_sample = sce_one_cluster[, sel_cells]
      S = rowSums(assays(sce_one_sample)$spliced)
      U = rowSums(assays(sce_one_sample)$unspliced)
      A = rowSums(assays(sce_one_sample)$ambiguous)
      
      SUA[[i]] = cbind(S, U, A)
    }
    x = do.call(cbind, SUA)
    
    # start PI_counts from relative counts estimates
    # start PI_SU from 50:50 or relative abundance from SU counts estimates
    genes = rownames(sce_one_cluster)
    n_genes = length(genes)
    
    sel_A = 1:9 # select columns of group A
    sel = (rowSums(x[,sel_A]) >= min_counts_per_group) & 
      (rowSums(x[,-sel_A]) >= min_counts_per_group)
    sel_genes = genes[sel]
    n_genes = length(sel_genes)
    sum(sel)
    
    # filter transcripts based on sel
    sce_one_cluster = sce_one_cluster[sel,]
    
    SUA = rbind(as.matrix(assays(sce_one_cluster)$spliced),
                as.matrix(assays(sce_one_cluster)$unspliced),
                as.matrix(assays(sce_one_cluster)$ambiguous))
    
    # prepare design matrix
    spliced_unspliced_ambiguous <- factor(c(rep("S", n_genes), rep("U", n_genes),  rep("A", n_genes)))
    
    rownames(SUA) = paste( rep(sel_genes, 3),
                           spliced_unspliced_ambiguous,
                           sep  = "-")
    ############################################################
    # satuRn:
    ############################################################
    # Load truth table:
    design <- data.frame(row.names = colnames(SUA),
                         sample_id = sce_one_cluster$sample_id,
                         group = sce_one_cluster$group)
    
    # the columns with transcript identifiers is names isoform_id, 
    # while the column containing gene identifiers should be named gene_id
    txInfo = data.frame(isoform_id = rownames(SUA),
                        gene_id = rep(sel_genes, 3))
    
    # Create a SummarizedExperiment object
    sumExp <- SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = SUA),
      colData = design,
      rowData = txInfo
    )
    
    # for sake of completeness: specify design formula from colData
    metadata(sumExp)$formula <- ~ 0 + 
      as.factor(colData(sumExp)$group) +
      as.factor(colData(sumExp)$sample_id)
    sumExp
    
    BPPARAM = MulticoreParam(6)
    
    # We fit the model
    # ?satuRn::fitDTU
    sumExp <- satuRn::fitDTU(
      object = sumExp,
      formula = ~ 0 + group + sample_id,
      parallel = FALSE,
      BPPARAM = BPPARAM,
      verbose = TRUE
    )
    
    # construct design matrix
    group <- as.factor(design$group)
    sample_id = as.factor(design$sample_id)
    design_full <- model.matrix(~ 0 + group + sample_id)
    colnames(design_full)[1:2] <- c("A", "B")
    
    # initialize contrast matrix
    L <- limma::makeContrasts(
      Contrast1 = A - B,
      levels = design_full
    )
    
    # We test the genes
    # ?satuRn::testDTU
    sumExp <- satuRn::testDTU(
      object = sumExp,
      contrasts = L,
      diagplot1 = TRUE,
      diagplot2 = TRUE,
      sort = FALSE
    )
    res <- rowData(sumExp)[["fitDTUResult_Contrast1"]]
    DF_SatuRn[[cl]] = data.frame(empirical_pval = res$empirical_pval, 
                           transcript_id = rownames(res),
                           cluster_id = clusters[cl])
  }
})
############################################################
# save results:
############################################################

full_name = file.path("02_results/SatuRnSCsampleDesign.RData")

save(DF_SatuRn, TIMES, file = full_name)
