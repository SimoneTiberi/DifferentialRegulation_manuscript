alias R='R_LIBS=/home/stiber/R/x86_64-pc-linux-gnu-library/4.3 /home/Shared/simone/DifferentialRegulation/software/R-4.3.0/bin/R'

cd /home/Shared_sherborne/simone/Diff_Velo/discovery

R
rm(list = ls())

library(SingleCellExperiment)
library(muscat)
library(eisaR)

TIMES = system.time({
  # load USA counts with pre-computed cell-type information:
  name = file.path("sce_cell_type.rds")
  sce = readRDS(name)
  
  # assign ambiguous counts to S or U
  assays(sce)$spliced   = assays(sce)$spliced + 0.5 * assays(sce)$ambiguous
  assays(sce)$unspliced = assays(sce)$unspliced + 0.5 * assays(sce)$ambiguous
  
  # define list of all possible group separations:
  sample_ids <- levels(sce$sample_id)
  group = ifelse(sample_ids %in% paste0("organoid",1:3), "3mon", "6mon")
  
  clusters = unique(sce$cell_types)
  clusters
  
  sort(table(sce$cell_type))
  
  min_counts_per_group = 10
  
  DF_eisaR =  list()
  
  for(cl in seq_along(clusters)){
    cluster = clusters[cl]
    
    # filter cell types:
    sce_one_cluster = sce[,sce$cell_types == cluster ]
    
    # prepare sce for pseudobulk
    sce_one_cluster <- prepSCE(sce_one_cluster, 
                               kid = "cell_types", # subpopulation assignments
                               sid = "sample_id", # sample IDs (ctrl/stim.1234)
                               gid = "group",
                               drop = TRUE) # drop all other colData columns
    
    # compute pseudo-bulk (PB) spliced and unspliced counts
    pb_S <- aggregateData(sce_one_cluster,
                          assay = "spliced", 
                          fun = "sum",
                          by = c("cluster_id", "sample_id"))
    
    pb_U <- aggregateData(sce_one_cluster,
                          assay = "unspliced", 
                          fun = "sum",
                          by = c("cluster_id", "sample_id"))
    
    # run eisaR
    R_ex <- assays(pb_S)[[1]]
    R_in <- assays(pb_U)[[1]]
    
    # min "min_counts_per_group" counts in each group:
    x = R_ex + R_in
    sel_A = group == "3mon" # select columns of group A
    sel = (rowSums(x[,sel_A]) >= min_counts_per_group) & (rowSums(x[,-sel_A]) >= min_counts_per_group)
    R_ex = R_ex[sel,]
    R_in = R_in[sel,]
    
    set.seed(169612)
    tab <- runEISA(cntEx  = R_ex, cntIn = R_in, cond = group, geneSelection = "none")$tab.ExIn
    
    # initialize results data frame
    DF_eisaR[[cl]] <- data.frame(Gene_id = rownames(tab),
                                 Cluster_id = cluster,
                                 p_eisaR = tab$PValue,
                                 p_eisaR_adj = p.adjust(tab$PValue, method = "BH"))
  }
  DF_eisaR = do.call(rbind, DF_eisaR)
})

full_name = file.path("02_results/eisaR.RData")

save(DF_eisaR, TIMES, file = full_name)
