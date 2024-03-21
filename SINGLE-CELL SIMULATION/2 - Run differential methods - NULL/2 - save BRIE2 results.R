alias R='R_LIBS=/home/stiber/R/x86_64-pc-linux-gnu-library/4.3 /home/Shared/simone/DifferentialRegulation/software/R-4.3.0/bin/R'

cd /home/Shared_sherborne/simone/Diff_Velo/NEW_simulation/06_alevin_fry_USA/BRIE2_USA

R
rm(list = ls())

library(SingleCellExperiment)

CLUSTERS <- c("Adipocytes", "Epithelial_cells", "Hepatocytes")

cols <- c("GeneID", "is_A_pval", "is_A_FDR")

RES <- lapply(CLUSTERS, FUN = function(cell_type) {
  name = paste0("isA_", cell_type, ".brie_ident.tsv")
  temp <- read.table(file = paste0("simulation_NULL/", name),
                     sep = "\t",
                     header = TRUE)[, cols]
  temp$Cell_type <- cell_type
  colnames(temp) <- c("Gene_id", "p_BRIE2", "p_BRIE2_adj", "Cell_type")
  
  return(temp)
})
RES <- do.call("rbind", RES)

real  =  1498 * 60 + 34.205
user  =  8099 * 60 + 2.286
sys   =  142 * 60 + 57.664

TIME = c(user, sys, real)

save(RES, TIME, file = "../../07_results/BRIE2_USA_NULL.RData")

rm(RES); rm(TIME)
