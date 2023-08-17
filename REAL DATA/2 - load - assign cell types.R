alias R='R_LIBS=/home/stiber/R/x86_64-pc-linux-gnu-library/4.2 /usr/local/R/R-4.2.0/bin/R'

# velo folder:
cd /home/Shared_sherborne/simone/Diff_Velo/RNAVelo/brain_human/05_alevin_fry

R
rm(list =ls())

library(SingleCellExperiment)
library(DifferentialRegulation)
library(scater)


sample_ids <- paste0("organoid", c(1:3, 16:18) )

# set directories of each sample input data (obtained via alevin-fry):
file.exists(sample_ids)

# set paths to objects:
path_to_counts = file.path(sample_ids,"/alevin/quants_mat.mtx")
path_to_cell_id = file.path(sample_ids,"/alevin/quants_mat_rows.txt")
path_to_gene_id = file.path(sample_ids,"/alevin/quants_mat_cols.txt")
path_to_EC_counts = file.path(sample_ids,"/alevin/geqc_counts.mtx")
path_to_EC = file.path(sample_ids,"/alevin/gene_eqclass.txt.gz")

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

########################################################
# QC and cell filtering
########################################################
# calculate total counts
assays(sce)$TOT_counts = assays(sce)$spliced + assays(sce)$unspliced + assays(sce)$ambiguous

########################################################
# remove undetected genes
########################################################
sce <- sce[rowSums(assays(sce)$TOT_counts > 0) > 0, ]

########################################################
# remove low quality cells
########################################################
# calculate per-cell quality control (QC) metrics
qc <- perCellQCMetrics(sce)

# remove cells with few or many detected genes
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
sce <- sce[, !ol]

########################################################
# remove lowly expressed genes: at least 10 non-zero cells
########################################################
sce <- sce[rowSums(assays(sce)$TOT_counts > 0) >= 10, ]

########################################################
# compute sum-factors & normalize
########################################################
sce <- computeLibraryFactors(sce)
sce <- logNormCounts(sce)

########################################################
# assign cell types
########################################################
# assign cell types
md <- read.csv("/home/Shared_sherborne/simone/Diff_Velo/results/discovery/meta_combined.txt", sep = "\t")
md <- md[grepl("PGP1", md$Batch), ]
md$SEQ <- sapply(md$NAME, FUN = function (x) strsplit(x, split = "_")[[1]][3])

table(md$Organoid)
# keep organoids 1:3 and 16:18 only
md = md[ md$Organoid %in% c("1", "2", "3", "16", "17", "18"), ]
table(md$Organoid)

# rename cell id (to include sample name):
md$ID <- paste0("organoid",md$Organoid, ".", md$SEQ)
md <- md[, c("ID", "CellType")]
head(md$ID); tail(md$ID)

matches = match(colnames(sce), md$ID)
head(colnames(sce)); head(md$ID[matches]) # OK!
tail(colnames(sce)); tail(md$ID[matches]) # OK!

sce$cell_types <- md$CellType[matches]
rm(md); rm(matches)

# count how many there are before filtering them out (if too many smt wrong):
sel_NAs = is.na(sce$cell_types)
sum( sel_NAs ); mean( sel_NAs)
# [1] 0.1640522

# filter out cells without cell types
sce <- sce[, !sel_NAs]
rm(sel_NAs)

sel_Unknown = sce$cell_types == "Unknown"
sum( sel_Unknown ); mean( sel_Unknown)
# 0.01539293

# filter out cells with unknown cell type
sce <- sce[, !sel_Unknown]
rm(sel_Unknown)

########################################################
# filter sce with cell types present in both groups:
########################################################
# define group
sce$group = ifelse(sce$sample_id %in% paste0("organoid",1:3), "3mon", "6mon")
table(sce$group, sce$sample_id)

tab = table(sce$cell_types, sce$group)
tab
# at least 100 cells in both groups
sel = rowSums(tab >= 100) == 2

sel_cell_types = rownames(tab)[sel]
sel_cell_types

sel = sce$cell_types %in% sel_cell_types

sce = sce[,sel]

table(sce$cell_types, sce$sample_id)

rm(sel)

########################################################
# UMAP:
########################################################
sce <- runUMAP(sce)

AA = plotUMAP(sce, colour_by = "sample_id")

BB = plotUMAP(sce, colour_by = "cell_types")

CC = plotUMAP(sce, colour_by = "group")

ggsave(filename = "plots/UMAP_sample.pdf",
       plot = AA,
       device = "pdf",
       width = 8,
       height = 8,
       units = "in",
       dpi = 300,
       limitsize = TRUE)

ggsave(filename = "plots/UMAP_cell.pdf",
       plot = BB,
       device = "pdf",
       width = 8,
       height = 8,
       units = "in",
       dpi = 300,
       limitsize = TRUE)

ggsave(filename = "plots/UMAP_group.pdf",
       plot = CC,
       device = "pdf",
       width = 8,
       height = 8,
       units = "in",
       dpi = 300,
       limitsize = TRUE)

########################################################
# remove unnecessary objects from sce
########################################################
assays(sce)$logcounts = NULL
assays(sce)$counts = NULL
assays(sce)$TOT_counts = NULL

########################################################
# save:
########################################################
name = file.path("/home/Shared_sherborne/simone/Diff_Velo/results/discovery/sce_cell_type.rds")

saveRDS(sce, file = name)

