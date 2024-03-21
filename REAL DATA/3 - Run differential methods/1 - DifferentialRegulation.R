alias R='R_LIBS=/home/stiber/R/x86_64-pc-linux-gnu-library/4.3 /home/Shared/simone/DifferentialRegulation/software/R-4.3.0/bin/R'

cd /home/Shared_sherborne/simone/Diff_Velo/discovery/05_alevin_fry

R

rm(list = ls())

seed_123 = FALSE

TIMES = system.time({
  library(DifferentialRegulation)
  library(SingleCellExperiment)
  
  sample_ids <- paste0("organoid", c(1:3, 16:18) )
  
  # set directories of each sample input data (obtained via alevin-fry):
  file.exists(sample_ids)
  
  # set paths to objects:
  path_to_cell_id = file.path(sample_ids,"/alevin/quants_mat_rows.txt")
  path_to_gene_id = file.path(sample_ids,"/alevin/quants_mat_cols.txt")
  path_to_EC_counts = file.path(sample_ids,"/alevin/geqc_counts.mtx")
  path_to_EC = file.path(sample_ids,"/alevin/gene_eqclass.txt.gz")
  
  file.exists(path_to_cell_id)
  file.exists(path_to_gene_id)
  file.exists(path_to_EC_counts)
  file.exists(path_to_EC)
  
  # load USA counts with pre-computed cell-type information:
  name = file.path("/home/Shared_sherborne/simone/Diff_Velo/results/discovery/sce_cell_type.rds")
  sce = readRDS(name)
  
  sum(assays(sce)$ambiguous)/ ( sum(assays(sce)$spliced) + sum(assays(sce)$unspliced) + sum(assays(sce)$ambiguous) )
  # 0.1563477
  
  #GB = 1073741824
  #object.size(sce)/GB
  # 1.1 GB
  
  # load EC counts:
  EC_list = load_EC(path_to_EC_counts,
                    path_to_EC,
                    path_to_cell_id,
                    path_to_gene_id,
                    sample_ids)
  # The percentage of multi-gene mapping reads in sample 'organoid1' is: 21.46
  # The percentage of multi-gene mapping reads in sample 'organoid2' is: 21.22
  # The percentage of multi-gene mapping reads in sample 'organoid3' is: 21.18
  # The percentage of multi-gene mapping reads in sample 'organoid16' is: 23.43
  # The percentage of multi-gene mapping reads in sample 'organoid17' is: 19.63
  # The percentage of multi-gene mapping reads in sample 'organoid18' is: 18.53
  
  # define the design of the study:
  design = data.frame(sample = sample_ids,
                      group = ifelse(sample_ids %in% paste0("organoid",1:3), "3mon", "6mon") )
  design
  
  # compute PB counts
  PB_counts = compute_PB_counts(sce = sce,
                                EC_list = EC_list,
                                design =  design,
                                sample_col_name = "sample",
                                group_col_name = "group",
                                sce_cluster_name = "cell_types",
                                min_cells_per_cluster = 0, 
                                min_counts_per_gene_per_group = 10)
  rm(sce); rm(EC_list)
  
  # EC-based test:
  if(seed_123){
    set.seed(123)
  }else{
    set.seed(169612)
  }
  results_EC = DifferentialRegulation(PB_counts,
                                      n_cores = 6,
                                      trace = TRUE)
})

if(seed_123){
  full_name = file.path("../02_results/DifferentialRegulation_WITH_trace_seed123.RData")
}else{
  full_name = file.path("../02_results/DifferentialRegulation_WITH_trace.RData")
}

save(results_EC, TIMES, file = full_name)
