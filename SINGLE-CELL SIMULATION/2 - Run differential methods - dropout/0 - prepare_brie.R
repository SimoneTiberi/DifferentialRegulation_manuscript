alias R='R_LIBS=/home/stiber/R/x86_64-pc-linux-gnu-library/4.3 /home/Shared/simone/DifferentialRegulation/software/R-4.3.0/bin/R'

cd /home/Shared_sherborne/simone/Diff_Velo/NEW_simulation/06_alevin_fry_USA

R

# clean environment
rm(list = ls())

# load packages
suppressPackageStartupMessages({
  # load libraries
  library(anndata)
  library(fishpond)
  library(SingleCellExperiment)
})

# load abundances
sample_ids <- 1:4
cell_types <- c("Adipocytes", "Epithelial_cells", "Hepatocytes")

load_abundances <- function (sample_id, cell_type, drop) {
  path <- ifelse(drop == 95, paste0("simulation_drop95/", paste0("normal", sample_id),
                             "_", cell_type), 
                 paste0("simulation_drop99/", paste0("normal", sample_id),
                        "_", cell_type))
  temp <- loadFry(path, outputFormat = "raw", nonzero = FALSE, quiet = TRUE)
  temp$sample_id <- sample_id
  colnames(temp) <- paste0(sample_id, ".", colnames(temp))
  return (temp)
}

clusters = c("Adipocytes", "Epithelial_cells", "Hepatocytes")

for(drop in c(95,99)){
  for(cl in clusters){
    sces <- lapply(sample_ids, load_abundances, cell_type = cl, drop = drop)
    
    sces <- unlist(sces)
    sce <- do.call("cbind", sces)
    assays(sce)
    
    assays(sce)$TOT_counts <- assays(sce)$spliced + assays(sce)$unspliced + assays(sce)$ambiguous
    
    # remove undetected genes
    sce <- sce[rowSums(assays(sce)$TOT_counts > 0) > 0, ]
    
    # set group parameter
    GROUP <- c("A", "B", "A", "B")
    
    # FILTER SCE BEFORE RUNNING BRIE2:
    min_count = 10
    
    # set group attribute
    sce$group <- ifelse(sce$sample_id %in% which(GROUP == "A"), "A", "B")
    
    # remove lowly expressed genes: at least 10 non-zero cells:
    filter <- rowSums(assays(sce)$TOT_counts[, sce$group == "A"]) >= min_count & 
      rowSums(assays(sce)$TOT_counts[, sce$group == "B"]) >= min_count
    
    sce <- sce[filter, ]
    
    assays(sce)$TOT_counts = NULL
    assays(sce)$counts = assays(sce)$spliced
    
    library(zellkonverter)
    ad = SCE2AnnData(sce,
                     X_name = "counts")
    
    if (drop == 95) {
      ad$write_h5ad(paste0("BRIE2_USA/simulation_drop95/cell_info_",
                           cl, ".h5ad"))
      write.table(data.frame(index = colnames(sce), 
                             is_A = ifelse(sce$group == "A", 1, 0)),
                  file = paste0("BRIE2_USA/simulation_drop95/cell_info_", 
                                cl, ".tsv"),
                  quote = FALSE,
                  sep = "\t", row.names = FALSE)
    } else {
      ad$write_h5ad(paste0("BRIE2_USA/simulation_drop99/cell_info_",
                           cl, ".h5ad"))
      write.table(data.frame(index = colnames(sce), 
                             is_A = ifelse(sce$group == "A", 1, 0)),
                  file = paste0("BRIE2_USA/simulation_drop99/cell_info_", 
                                cl, ".tsv"),
                  quote = FALSE,
                  sep = "\t", row.names = FALSE)
    }
  }
}
