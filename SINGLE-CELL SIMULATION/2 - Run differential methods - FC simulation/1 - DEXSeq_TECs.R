alias R='R_LIBS=/home/stiber/R/x86_64-pc-linux-gnu-library/4.3 /home/Shared/simone/DifferentialRegulation/software/R-4.3.0/bin/R'

cd /home/Shared_sherborne/simone/Diff_Velo/NEW_simulation/06_alevin_fry_USA

R
rm(list = ls())

library(DifferentialRegulation)
library(DEXSeq)
library(SingleCellExperiment)
library(BiocParallel)

# specify 4 samples ids:
sample_ids = paste0("normal", seq_len(4))

min_counts_per_group = 10

clusters = c("Adipocytes", "Epithelial_cells", "Hepatocytes")

for(FC in c(6, 9)){
  TIMES = list()  
  
  if(FC == 6){
    data_dir = "simulation_FC6"
  }else{
    data_dir = "simulation_FC9"
  }
  
  DF_DEXSeq = list()
  
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
      path_to_EC_counts = file.path(base_dir,"/alevin/geqc_counts.mtx")
      path_to_EC = file.path(base_dir,"/alevin/gene_eqclass.txt.gz")
      
      file.exists(path_to_counts)
      file.exists(path_to_cell_id)
      file.exists(path_to_gene_id)
      file.exists(path_to_EC_counts)
      file.exists(path_to_EC)
      
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
      
      # set parallel cores:
      #library(BiocParallel)
      #BPPARAM = MulticoreParam(6)
      
      # analyze each cluster separately:
      dxd <- DEXSeqDataSet(countData = round(SUA),
                           sampleData = design,
                           design = ~sample + exon + condition:exon,
                           featureID = spliced_unspliced_ambiguous,
                           groupID = rep(sel_genes, 3))
      
      dxd <- estimateSizeFactors(dxd)
      dxd <- estimateDispersions(dxd) #, BPPARAM = BPPARAM)
      dxd <- testForDEU(dxd, reducedModel = ~sample + exon) #, BPPARAM = BPPARAM)
      
      dxr <- DEXSeqResults(dxd, independentFiltering = FALSE)
      
      gene_q_val = perGeneQValue(dxr)
      
      print(head(gene_q_val))
      
      DF_DEXSeq[[cl]] = data.frame(qval = gene_q_val, gene_id = names(gene_q_val), cluster_id = clusters[cl])
    })
  }
  
  DF_DEXSeq = do.call(rbind, DF_DEXSeq)
  
  if(FC == 6){
    name = paste0("results_FC6_DEXSeq_USA.RData")
  }else{
    name = paste0("results_FC9_DEXSeq_USA.RData")
  }
  
  full_name = file.path("../07_results", name)
  
  save(DF_DEXSeq, TIMES, file = full_name)
}

