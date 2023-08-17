alias R='R_LIBS=/home/stiber/R/x86_64-pc-linux-gnu-library/4.3 /home/Shared/simone/DifferentialRegulation/software/R-4.3.0/bin/R'

cd /home/Shared_sherborne/simone/Diff_Velo/discovery

R
rm(list = ls())

TIMES = system.time({
  library(DRIMSeq)
  library(SingleCellExperiment)
  library(BiocParallel)
  
  DF_DRIMSeq = list()
  
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
    
    x = x[sel,]
    dim(x)
    S = x[, seq(1, length.out = n_samples, by = 3) ]
    U = x[, seq(2, length.out = n_samples, by = 3) ]
    A = x[, seq(3, length.out = n_samples, by = 3) ]
    head(S, 20); head(U, 20); head(A, 20)
    colSums(S);
    colSums(U);
    colSums(A);
    
    SUA = rbind(S, U, A)
    dim(SUA); colSums(SUA)
    
    # prepare design matrix
    spliced_unspliced_ambiguous <- factor(c(rep("S", n_genes), rep("U", n_genes),  rep("A", n_genes)))
    
    rownames(SUA) = paste( rep(sel_genes, 3),
                           spliced_unspliced_ambiguous,
                           sep  = "-")
    ############################################################
    # DRIMSeq:
    ############################################################
    # Load truth table:
    colnames(SUA) = sample_ids
    
    design <- data.frame(sample_id = sample_ids,
                         group = group)
    
    # save transcripts as "gene_id"
    counts_df = data.frame(SUA, 
                           gene_id = rep(sel_genes, 3),
                           feature_id = rownames(SUA))
    
    # Create a dmDSdata object
    d <- dmDSdata(counts = counts_df, samples = design)
    # with 46647 genes and 6 samples
    
    design_full <- model.matrix(~ group, data = DRIMSeq::samples(d))
    design_full
    
    # set parallel cores:
    BPPARAM = MulticoreParam(6)
    
    # infer the precision parameters:
    d <- dmPrecision(d, genewise_precision = TRUE, 
                     design = design_full,
                     BPPARAM = BPPARAM)
    
    # We fit the model
    # ?dmFit
    d <- dmFit(d, design = design_full, 
               verbose = 1,
               BPPARAM = BPPARAM)
    
    # We test the genes
    # ?dmTest
    d <- dmTest(d, coef = "group6mon", verbose = 1,
                BPPARAM = BPPARAM)
    
    results_gene = results(d, level = "gene")
    
    # initialize results data frame
    DF_DRIMSeq[[cl]] <- data.frame(Gene_id = results_gene$gene_id,
                                   Cluster_id = cluster,
                                   p_DRIMSeq = results_gene$pvalue,
                                   p_DRIMSeq_adj = results_gene$adj_pvalue)
  }
  DF_DRIMSeq = do.call(rbind, DF_DRIMSeq)
})
############################################################
# save results:
############################################################

full_name = file.path("02_results/DRIMSeq.RData")

save(DF_DRIMSeq, TIMES, file = full_name)
