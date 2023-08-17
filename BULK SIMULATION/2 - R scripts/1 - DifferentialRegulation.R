alias R='R_LIBS=/home/stiber/R/x86_64-pc-linux-gnu-library/4.3 /home/Shared/simone/DifferentialRegulation/software/R-4.3.0/bin/R'

cd /home/Shared/simone/DifferentialRegulation

R

rm(list =ls())
library(DifferentialRegulation)

undersampling_int = 10

N_MCMC = 2000
burn_in = N_MCMC/4

############################################################
# choose directories:
############################################################
# TRANSCRIPT ESTIMATED COUNTS:
directories = c("de0_ds0_dr0",
                "de0_ds0_dr2000",
                "de2000_ds0_dr2000",
                "de0_ds2000_dr2000")

directories_results = c("NULL",
                        "DR",
                        "DR + DGE",
                        "DR + DAS")

############################################################
# load data:
############################################################
library(tximport)
TIME = list()
for(type in 2:4){
  TIME[[type]] = system.time({
    set.seed(169612)
    
    data_dir = file.path("1 - data", directories[type],
                         "2_quants/salmon" )
    
    sample_names = paste0("sample", seq_len(6))
    
    quant_files = file.path(data_dir, sample_names, "quant.sf")
    equiv_classes_files = file.path(data_dir, sample_names, "aux_info/eq_classes.txt.gz")
    
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
                                         undersampling_int = undersampling_int,
                                         N_MCMC = N_MCMC,
                                         burn_in = burn_in)
  })
  ############################################################
  # save results:
  ############################################################
  data_dir = file.path("3 - results",
                       directories_results[type])
  
  name = paste("DifferentialRegulation")
  
  name = paste0(name, "_", N_MCMC)
  
  name = paste0(name, "_", undersampling_int)
  
  name = paste0(name, ".RData")
  
  save(res_EC, file = file.path(data_dir,name))
}
#OUTSIDE THE LOOP, save time:
name = paste0("TIME_", name)
save(TIME, file = file.path("3 - results", name))

do.call(rbind, TIME[2:4])

mean(do.call(rbind, TIME[2:4])[,3])/60
# 63.56871 # undersampling_int [1] 10
# 245.4225 # undersampling_int [1] 1
245.4225/63.56871
3.860744 faster

1-63.56871/245.4225
# 75/% speed up