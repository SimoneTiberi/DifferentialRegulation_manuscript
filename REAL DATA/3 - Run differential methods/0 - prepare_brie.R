alias R='R_LIBS=/home/stiber/R/x86_64-pc-linux-gnu-library/4.3 /home/Shared/simone/DifferentialRegulation/software/R-4.3.0/bin/R'

cd /home/Shared_sherborne/simone/Diff_Velo/discovery

R

# clean environment
rm(list = ls())

# load packages
suppressPackageStartupMessages({
  # load libraries
  library(SingleCellExperiment)
})

# load USA counts with pre-computed cell-type information:
name = file.path("sce_cell_type.rds")
sce = readRDS(name)

# define list of all possible group separations:
sample_ids <- levels(sce$sample_id)
group = ifelse(sample_ids %in% paste0("organoid",1:3), "3mon", "6mon")

clusters = unique(sce$cell_types)
clusters

reducedDim(sce) = NULL
colData(sce)$sizeFactor = NULL

min_counts_per_group = 10

for(cl in seq_along(clusters)){
  cluster = clusters[cl]
  
  # filter cell types:
  sce_one_cluster = sce[,sce$cell_types == cluster ]
  
  assays(sce_one_cluster)$TOT_counts <- assays(sce_one_cluster)$spliced + assays(sce_one_cluster)$unspliced + assays(sce_one_cluster)$ambiguous
  
  # remove undetected genes
  sce_one_cluster <- sce_one_cluster[rowSums(assays(sce_one_cluster)$TOT_counts > 0) > 0, ]
  
  # FILTER SCE BEFORE RUNNING BRIE2:
  min_count = 10
  
  # remove lowly expressed genes: at least 10 non-zero cells:
  filter <- rowSums(assays(sce_one_cluster)$TOT_counts[, sce_one_cluster$group == "3mon"]) >= min_count & 
    rowSums(assays(sce_one_cluster)$TOT_counts[, sce_one_cluster$group == "6mon"]) >= min_count
  
  sce_one_cluster <- sce_one_cluster[filter, ]
  
  assays(sce_one_cluster)$TOT_counts = NULL
  assays(sce_one_cluster)$counts = assays(sce_one_cluster)$spliced
  
  assays(sce_one_cluster)$spliced = assays(sce_one_cluster)$spliced
  assays(sce_one_cluster)$unspliced = assays(sce_one_cluster)$unspliced
  assays(sce_one_cluster)$ambiguous = assays(sce_one_cluster)$ambiguous

  # ERROR -> rowData was missing!!
  # it needs to be added to the sce!
  rowData(sce_one_cluster) = data.frame(rownames = rownames(sce_one_cluster))
  
  library(zellkonverter)
  ad = SCE2AnnData(sce_one_cluster,
                   X_name = "counts")

  ad$write_h5ad(paste0("BRIE2_USA/", cluster, ".h5ad"))
  
  write.table(data.frame(index = colnames(sce_one_cluster), 
                         is_A = ifelse(sce_one_cluster$group == "3mon", 1, 0)),
              file = paste0("BRIE2_USA/", cluster, ".tsv"),
              quote = FALSE,
              sep = "\t", row.names = FALSE)
}
