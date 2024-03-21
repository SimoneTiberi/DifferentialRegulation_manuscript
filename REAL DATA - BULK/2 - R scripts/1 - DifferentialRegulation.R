alias R='R_LIBS=/home/stiber/R/x86_64-pc-linux-gnu-library/4.3 /home/Shared/simone/DifferentialRegulation/software/R-4.3.0/bin/R'

cd /home/Shared/simone/Diff_Velo/discovery_bulk

R

rm(list =ls())

seed_123 = FALSE

library(DifferentialRegulation)

############################################################
# load data:
############################################################
library(tximport)
TIME = system.time({
  if(seed_123){
    set.seed(123)
  }else{
    set.seed(169612)
  }
  data_dir = "03_Salmon"
  
  sample_names = paste0("ERR45792", 67:72)
  quant_files = file.path(data_dir, sample_names, "quant.sf")
  file.exists(quant_files)
  
  equiv_classes_files = file.path(data_dir, sample_names, "aux_info/eq_classes.txt.gz")
  file.exists(equiv_classes_files)
  
  ############################################################
  # DR pipeline:
  ############################################################
  # load EC:
  EC_list = load_bulk_EC(path_to_eq_classes = equiv_classes_files)
  
  # load US estimated counts:
  sce = load_bulk_US(quant_files,
                     sample_names)
  
  group_names = rep(c("A", "B"), each = 3)
  design = data.frame(sample = sample_names,
                      group = group_names)
  design
  
  res_EC = DifferentialRegulation_bulk(sce, 
                                       EC_list = EC_list,
                                       design = design, 
                                       n_cores = 1,
                                       trace = TRUE)
})
############################################################
# save results:
############################################################
data_dir = "04_results"

if(seed_123){
  name = paste("DifferentialRegulation_WITH_trace_seed123.RData")
}else{
  name = paste("DifferentialRegulation_WITH_trace.RData")
}

save(res_EC, file = file.path(data_dir,name))

#OUTSIDE THE LOOP, save time:
name = paste0("TIME_", name)
save(TIME, file = file.path(data_dir, name))
