alias R='R_LIBS=/home/stiber/R/x86_64-pc-linux-gnu-library/4.3 /home/Shared/simone/DifferentialRegulation/software/R-4.3.0/bin/R'

cd /home/Shared_sherborne/simone/Diff_Velo/NEW_simulation/06_alevin_fry_USA

R

rm(list = ls())

undersampling_int = 10

N_MCMC = 2000
burn_in = N_MCMC/4

library(DifferentialRegulation)

# specify 4 samples ids:
sample_ids = paste0("normal", seq_len(4))

results_EC = list()

data_dir = "simulation_NULL"

TIMES = list()  

clusters = c("Adipocytes", "Epithelial_cells", "Hepatocytes")


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
    
    # load EC counts:
    EC_list = load_EC(path_to_EC_counts,
                      path_to_EC,
                      path_to_cell_id,
                      path_to_gene_id,
                      sample_ids)
    
    # define the design of the study:
    design = data.frame(sample = sample_ids,
                        group = c("A", "B", "A", "B") )
    design
    
    # add cell type information (1 cell type only):
    sce$cell_type = cluster
    
    # compute PB counts
    PB_counts = compute_PB_counts(sce = sce,
                                  EC_list = EC_list,
                                  design =  design,
                                  sample_col_name = "sample",
                                  group_col_name = "group",
                                  sce_cluster_name = "cell_type",
                                  min_cells_per_cluster = 0, 
                                  min_counts_per_gene_per_group = 10)
    rm(sce); rm(EC_list)
    
    # EC-based test:
    set.seed(16961)
    results_EC[[cl]] = DifferentialRegulation(PB_counts,
                                              n_cores = 1,
                                              undersampling_int = undersampling_int,
                                              N_MCMC = N_MCMC,
                                              burn_in = burn_in)
  })
  print(cl)
}

name = "DifferentialRegulation_NULL.RData"

full_name = file.path("../07_results", name)

save(results_EC, TIMES, file = full_name)
