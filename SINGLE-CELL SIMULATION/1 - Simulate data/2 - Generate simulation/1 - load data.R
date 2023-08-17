alias R='R_LIBS=/home/stiber/R/x86_64-pc-linux-gnu-library/4.3 /home/Shared/simone/DifferentialRegulation/software/R-4.3.0/bin/R'

cd /home/Shared_sherborne/simone/Diff_Velo/NEW_simulation

R

# clean environment
rm(list = ls())

# load packages
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scater)
  library(celldex)
  library(SingleR)
  library(tximeta)
})

# load abundances
loadLinkedTxome(paste0("01_annotation/", 
                       "gencode.vM24.annotation.expanded.json"))

import_abundances <- function (sample_id){
  txi <- tximeta(coldata = data.frame(
    names = sample_id,
    files = paste0("02_alevin/", sample_id, "/alevin/quants_mat.gz"), 
    stringsAsFactors = FALSE), type = "alevin")
  cg <- read.delim(paste0("01_annotation/", "gencode.vM24.annotation.expanded.features.tsv"),
                   header = TRUE, as.is = TRUE)
  colnames(cg)[colnames(cg) == "intron"] <- "unspliced"
  txis <- splitSE(txi, cg, assayName = "counts")
  txis <- as(txis, "SingleCellExperiment")
  assays(txis) <- list(
    counts = assay(txis, "spliced"),
    spliced = assay(txis, "spliced"),
    unspliced = assay(txis, "unspliced"))
  return (txis)
}

sample_ids <- paste0("normal", 1:4)
txis <- lapply(sample_ids, import_abundances)
names(txis) <- sample_ids

for (i in 1:4) {
  colnames(txis[[i]]) <- paste0(i, ".", colnames(txis[[i]]))
  colData(txis[[i]])$sample_id <- i
}

sce <- do.call("cbind", txis)
colData(sce)$sample_id <- factor(colData(sce)$sample_id)

assays(sce)$TOT_counts <- assays(sce)$spliced + assays(sce)$unspliced

# remove undetected genes
sce <- sce[rowSums(assays(sce)$TOT_counts > 0) > 0, ]

# calculate per-cell quality control (QC) metrics
qc <- perCellQCMetrics(sce)

# remove cells with few or many detected genes
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
sce <- sce[, !ol]

# remove lowly expressed genes: at least 10 non-zero cells:
sce <- sce[rowSums(assays(sce)$TOT_counts > 0) >= 10, ]

# compute sum-factors & normalize
sce <- computeLibraryFactors(sce)
sce <- logNormCounts(sce)

# assign gene names and cell types
ref <- celldex::MouseRNAseqData(ensembl = TRUE, cell.ont = "all")

# check how many gene names match btw the 2 datasets:
gene_ids <- rownames(sce)
rownames(sce) <- substr(rownames(sce), 1, 18)

pred <- SingleR(test = sce, ref = ref, labels = ref$label.main)
sce$cell_type <- pred$labels
rownames(sce) <- gene_ids

metadata(sce) <- list()

# run dimensionality reduction algorithm
sce <- runPCA(sce)
sce <- runUMAP(sce)

# save data
# saveRDS(sce, file = "03_data/mouse_data_alevin.rds")

#sce = readRDS("03_data/mouse_data_alevin.rds")

library(scater)
plotUMAP(sce, colour_by = "cell_type")


sce$sample_id = factor(sce$sample_id)
plotUMAP(sce, colour_by = "sample_id")

#sce_tmp = sce[,sce$sample_id %in% c(1,3)]
#plotUMAP(sce_tmp, colour_by = "sample_id")
