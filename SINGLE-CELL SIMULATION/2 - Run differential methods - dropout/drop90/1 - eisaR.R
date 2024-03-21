alias R='R_LIBS=/home/stiber/R/x86_64-pc-linux-gnu-library/4.3 /home/Shared/simone/DifferentialRegulation/software/R-4.3.0/bin/R'

cd /home/Shared_sherborne/simone/Diff_Velo/NEW_simulation/06_alevin_fry_USA

R
rm(list = ls())

library(DifferentialRegulation)
library(SingleCellExperiment)
library(muscat)
library(eisaR)

# specify 4 samples ids:
sample_ids = paste0("normal", seq_len(4))

min_counts_per_group = 10

clusters = c("Adipocytes", "Epithelial_cells", "Hepatocytes")

group = factor(c("A", "B", "A", "B"))

TIMES = list()
data_dir = "simulation_drop90"

DF_eisaR = list()

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
    
    assays(sce)$spliced = assays(sce)$spliced + 0.5 * assays(sce)$ambiguous
    assays(sce)$unspliced = assays(sce)$unspliced + 0.5 * assays(sce)$ambiguous
    
    assays(sce)$ambiguous = NULL
    assays(sce)$counts = NULL
    
    sce$group = 1
    sce$cluster = 1
    
    # prepare sce for pseudobulk
    sce <- prepSCE(sce, 
                   kid = "cluster", # subpopulation assignments
                   sid = "sample_id", # sample IDs (ctrl/stim.1234)
                   gid = "group",
                   drop = TRUE) # drop all other colData columns
    
    # compute pseudo-bulk (PB) spliced and unspliced counts
    pb_S <- aggregateData(sce,
                          assay = "spliced", 
                          fun = "sum",
                          by = c("cluster_id", "sample_id"))
    
    pb_U <- aggregateData(sce,
                          assay = "unspliced", 
                          fun = "sum",
                          by = c("cluster_id", "sample_id"))
    
    # run eisaR
    R_ex <- assays(pb_S)[[1]]; R_in <- assays(pb_U)[[1]]
    
    # min "min_counts_per_group" counts in each group:
    x = R_ex + R_in
    sel_A = group == "A" # select columns of group A
    sel = (rowSums(x[,sel_A]) >= min_counts_per_group) & (rowSums(x[,-sel_A]) >= min_counts_per_group)
    R_ex = R_ex[sel,]
    R_in = R_in[sel,]
    
    tab <- runEISA(cntEx  = R_ex, cntIn = R_in, cond = group, geneSelection = "none")$tab.ExIn
    
    # initialize results data frame
    DF_eisaR[[cl]] <- data.frame(Gene_id = rownames(tab),
                                 Cluster_id = cluster,
                                 p_eisaR = tab$PValue,
                                 p_eisaR_adj = p.adjust(tab$PValue, method = "BH"))
  })
}
DF_eisaR = do.call(rbind, DF_eisaR)

name = paste0("results_drop90_eisaR.RData")

full_name = file.path("../07_results", name)

save(DF_eisaR, TIMES, file = full_name)
