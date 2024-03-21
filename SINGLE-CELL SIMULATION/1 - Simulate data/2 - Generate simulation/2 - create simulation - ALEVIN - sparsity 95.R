# cd /home/Shared_sherborne/simone/Diff_Velo/NEW_simulation

setwd("~/Desktop/Differential Velocity/FINAL SCRIPTS/SINGLE-CELL SIMULATION")

library(SingleCellExperiment)

# clean-up environment
rm(list = ls())

# load original simulation ()
sce = readRDS("03_data/mouse_simulation_data-ALEVIN.rds")

# original rate of zeros (about 90%)
mean(assays(sce)$TOT_counts == 0)

##########################################
# CHANGE SPARSITY LEVELS:
##########################################
CLUSTERS <- c("Epithelial cells", "Adipocytes", "Hepatocytes")

SAMPLES = levels(sce$sample_id)

assays(sce)$logcounts = NULL
assays(sce)$counts = NULL

non_zero_cells = 0.05

sce_final = sce
for(cluster in CLUSTERS){
  set.seed(169612)
  sce_one_cluster = sce[, sce$cell_type == cluster]
  
  for(sample in SAMPLES){
    sce_one_cluster_one_sample = sce_one_cluster[, sce_one_cluster$sample_id == sample]
    # here loop inside sample id, to ensure each sample-cluster combination has the total overall S and U counts
    # we only change how they are distributed across cells, to have 1, 5 and 10 % of non-zero cells
    
    spliced = assays(sce_one_cluster_one_sample)$spliced
    unspliced = assays(sce_one_cluster_one_sample)$unspliced
    
    S = round(rowSums(spliced))
    U = round(rowSums(unspliced))
    
    n_genes = nrow(sce_one_cluster_one_sample)
    n_cells = ncol(sce_one_cluster_one_sample)
    
    # set all to zero:
    spliced[,] = 0
    unspliced[,] = 0
    
    for(gene in 1:n_genes){
      # for each gene, select XX % of non-zero cells
      n_non_zero_cells = round(n_cells * non_zero_cells)
      sel_cells = sample.int(n_cells, n_non_zero_cells)
      
      if(S[gene] > 0){
        spliced[gene,sel_cells] = rmultinom(n = 1, size = S[gene], prob = rep(1, n_non_zero_cells))
      }
      if(U[gene] > 0){
        unspliced[gene,sel_cells] = rmultinom(n = 1, size = U[gene], prob = rep(1, n_non_zero_cells))
      }
    }
    
    assays(sce_one_cluster_one_sample)$spliced = spliced
    assays(sce_one_cluster_one_sample)$unspliced = unspliced
    assays(sce_one_cluster_one_sample)$TOT_counts = spliced + unspliced
    
    sce_one_cluster[, sce_one_cluster$sample_id == sample] = sce_one_cluster_one_sample
    
    print(sample)
  }
  sce_final[, sce_final$cell_type == cluster] = sce_one_cluster
  
  print(cluster)
}

sce_final; sce

# original rate of zeros (about 90%)
mean(assays(sce)$TOT_counts == 0)

# NEW rate of zeros (defined above):
mean(assays(sce_final)$TOT_counts == 0)

# Total number of counts (should be similar)
sum(assays(sce)$TOT_counts);sum(assays(sce_final)$TOT_counts)

# Total number of counts (should be similar in each sample):
sum(assays(sce)$TOT_counts[sce$sample_id==1]); sum(assays(sce_final)$TOT_counts[sce_final$sample_id==1])
sum(assays(sce)$TOT_counts[sce$sample_id==2]); sum(assays(sce_final)$TOT_counts[sce_final$sample_id==2])
sum(assays(sce)$TOT_counts[sce$sample_id==3]); sum(assays(sce_final)$TOT_counts[sce_final$sample_id==3])
sum(assays(sce)$TOT_counts[sce$sample_id==4]); sum(assays(sce_final)$TOT_counts[sce_final$sample_id==4])

rm(sce); rm(sce_one_cluster_one_sample); rm(sce_one_cluster)
rm(spliced); rm(unspliced); rm(S); rm(U)

non_zero_cells
saveRDS(sce_final, file = "03_data/mouse_simulation_data-ALEVIN-dropout95.rds")

##########################################
# PLOTS
##########################################
# make and plot UMAP of the simulated data:
library(scater)
assays(sce_final)$counts = assays(sce_final)$spliced
sce_final <- computeLibraryFactors(sce_final)
sce_final <- logNormCounts(sce_final)
sce_final <- runUMAP(sce_final)

sce_final$sample_id = factor(sce_final$sample_id)
plotUMAP(sce_final, colour_by = "cell_type")
plotUMAP(sce_final, colour_by = "sample_id")
plotUMAP(sce_final, colour_by = "group")
# a group difference is not clear: it was introduced by the simulation
