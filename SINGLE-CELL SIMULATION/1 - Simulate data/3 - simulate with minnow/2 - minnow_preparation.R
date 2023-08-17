alias R='R_LIBS=/home/stiber/R/x86_64-pc-linux-gnu-library/4.3 /home/Shared/simone/DifferentialRegulation/software/R-4.3.0/bin/R'

cd /home/Shared_sherborne/simone/Diff_Velo/NEW_simulation

R

# clean environment
rm(list = ls())

# load libraries
suppressPackageStartupMessages({
  library(SingleCellExperiment)
})

# prepare for minnow simulation
sim_dir <- "/home/Shared_sherborne/simone/Diff_Velo/NEW_simulation/04_simulation_minnow/simulation/"
sim_dir_DGE <- "/home/Shared_sherborne/simone/Diff_Velo/NEW_simulation/04_simulation_minnow/simulation_DGE/"

sce <- readRDS("/home/Shared_sherborne/simone/Diff_Velo/NEW_simulation/03_data/mouse_simulation_data-ALEVIN.rds")
sce_DGE <- readRDS("/home/Shared_sherborne/simone/Diff_Velo/NEW_simulation/03_data/mouse_simulation_DGE_data-ALEVIN.rds")

cell_types <- c("Adipocytes", "Epithelial cells", "Hepatocytes")

prepare_minnow <- function (sce, sim_dir) {
  for (sample_id in 1:4) {
    for (cell_type in cell_types) {
      temp <- sce[, sce$sample_id == sample_id & sce$cell_type == cell_type]
      
      gnames <- rownames(temp)
      bcs <- sub(paste0(sample_id, "."), "", colnames(temp))
      
      cell_type <- ifelse(cell_type == "Epithelial cells", "Epithelial_cells", cell_type)
      
      s <- assay(temp, "spliced")
      u <- assay(temp, "unspliced")
      
      pat <- paste0(sample_id, ".")
      out_dir <- file.path(sim_dir, paste0("normal", sample_id), cell_type)
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      
      colnames(s) <- sub(pat, "", colnames(s))
      colnames(u) <- sub(pat, "", colnames(u))
      M <- as.matrix(rbind(s, u))
      
      write.table(M, file = file.path(out_dir, "quants_mat.csv"),
                  sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)
      write.table(c(gnames, paste0(gnames, "-U")), file = file.path(out_dir, "quants_mat_rows.txt"), 
                  quote = FALSE, row.names = FALSE, col.names = FALSE)
      write.table(bcs, file = file.path(out_dir, "quants_mat_cols.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
    }
  }
}

prepare_minnow(sce, sim_dir)
prepare_minnow(sce_DGE, sim_dir_DGE)