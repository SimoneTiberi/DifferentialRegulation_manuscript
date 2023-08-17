alias R='R_LIBS=/home/stiber/R/x86_64-pc-linux-gnu-library/4.3 /home/Shared/simone/DifferentialRegulation/software/R-4.3.0/bin/R'

cd /home/Shared_sherborne/simone/Diff_Velo/NEW_simulation/06_alevin_fry_USA

R
rm(list = ls())

library(DifferentialRegulation)
library(DRIMSeq)
library(SingleCellExperiment)
library(BiocParallel)

# specify 4 samples ids:
sample_ids = paste0("normal", seq_len(4))

min_counts_per_group = 10

clusters = c("Adipocytes", "Epithelial_cells", "Hepatocytes")

#for(DGE in c(FALSE,TRUE)){
for(DGE in c(FALSE)){
  
  TIMES = list()  
  
  if(DGE){
    data_dir = "simulation_DGE"
  }else{
    data_dir = "simulation"
  }
  
  DF_DRIMSeq = list()
  
  for(cl in seq_along(clusters)){
    
    TIMES[[cl]] = system.time({
      cluster = clusters[cl]
      
      # set directories of each sample input data (obtained via alevin-fry):
      base_dir = file.path(data_dir, paste0(sample_ids, "_", cluster) )
      file.exists(base_dir)
      
      # set paths to objects:
      path_to_counts = file.path(base_dir,"/alevin/quants_mat.mtx")
      path_to_cell_id = file.path(base_dir,"/alevin/quants_mat_rows.txt")
      path_to_gene_id = file.path(base_dir,"/alevin/quants_mat_cols.txt")
      
      file.exists(path_to_counts)
      file.exists(path_to_cell_id)
      file.exists(path_to_gene_id)
      
      # load USA counts:
      sce = load_USA(path_to_counts,
                     path_to_cell_id,
                     path_to_gene_id,
                     sample_ids)
      
      n_samples = nlevels(sce$sample_id); n_samples
      
      SUA = list()
      for(i in 1:n_samples){
        sel_cells = sce$sample_id == levels(sce$sample_id)[i]
        sce_one_sample = sce[, sel_cells]
        S = rowSums(assays(sce_one_sample)$spliced)
        U = rowSums(assays(sce_one_sample)$unspliced)
        A = rowSums(assays(sce_one_sample)$ambiguous)
        
        SUA[[i]] = cbind(S, U, A)
      }
      x = do.call(cbind, SUA)
      
      # start PI_counts from relative counts estimates
      # start PI_SU from 50:50 or relative abundance from SU counts estimates
      genes = rownames(sce)
      n_genes = length(genes)
      
      group = factor(c("A", "B", "A", "B"))
      sel_A = c(1:3, 7:9) # select columns of group A
      sel = (rowSums(x[,sel_A]) >= min_counts_per_group) & (rowSums(x[,-sel_A]) >= min_counts_per_group)
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
      design <- data.frame(condition = factor(group))
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
      
      # infer the precision parameters:
      d <- dmPrecision(d, genewise_precision = TRUE, 
                       design = design_full)
      
      # We fit the model
      # ?dmFit
      d <- dmFit(d, design = design_full, 
                 verbose = 1)
      
      # We test the genes
      # ?dmTest
      d <- dmTest(d, coef = "groupB", verbose = 1)
      
      results_gene = results(d, level = "gene")
      
      # initialize results data frame
      DF_DRIMSeq[[cl]] <- data.frame(Gene_id = results_gene$gene_id,
                                     Cluster_id = cluster,
                                     p_DRIMSeq = results_gene$pvalue,
                                     p_DRIMSeq_adj = results_gene$adj_pvalue)
    })
  }
  DF_DRIMSeq = do.call(rbind, DF_DRIMSeq)
  
  if(DGE){
    name = paste0("results_DGE_DRIMSeq.RData")
  }else{
    name = paste0("results_NO_DGE_DRIMSeq.RData")
  }
  
  full_name = file.path("../07_results", name)
  
  save(DF_DRIMSeq, TIMES, file = full_name)
}

