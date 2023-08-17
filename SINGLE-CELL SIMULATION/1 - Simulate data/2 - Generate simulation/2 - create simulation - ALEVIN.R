cd /home/Shared_sherborne/simone/Diff_Velo/NEW_simulation

# clean-up environment
rm(list = ls())

# load mouse data
sce <- readRDS("03_data/mouse_data_alevin.rds")

# filter cell types
CLUSTERS <- c("Epithelial cells", "Adipocytes", "Hepatocytes")
sce <- sce[, sce$cell_type %in% CLUSTERS]

# set arbitrary condition
GROUP <- c("A", "B", "A", "B")
sce$group <- ifelse(sce$sample_id %in% c("1", "3"), "A", "B")

# set parameters for simulation
p_genes <- 0.2
p_genes_DGE <- 0.2
mu_FC <- 3

assays(sce)

# load scripts to simulate data:
source("functions.R")

# simulate mouse data based on preset parameters
set.seed(123)
sim <- mouse_simulation(sce = sce, CLUSTERS = CLUSTERS, 
                        p_genes = p_genes)

# save simulated data sets
saveRDS(sim, file = "03_data/mouse_simulation_data-ALEVIN.rds")

set.seed(1234)
sim_2 <- mouse_simulation(sce = sce, CLUSTERS = CLUSTERS, 
                        p_genes = p_genes)
set.seed(12345)
sim_DGE <- add_DGE(sce = sim_2, CLUSTERS = CLUSTERS, 
                   p_genes_DGE = p_genes_DGE)

# save simulated data sets
saveRDS(sim_DGE, file = "03_data/mouse_simulation_DGE_data-ALEVIN.rds")

# make and plot UMAP of the simulated data:
library(scater)
assays(sim)$counts = assays(sim)$spliced
sim <- computeLibraryFactors(sim)
sim <- logNormCounts(sim)
sim <- runUMAP(sim)

sim$sample_id = factor(sim$sample_id)
plotUMAP(sim, colour_by = "cell_type")
plotUMAP(sim, colour_by = "sample_id")
plotUMAP(sim, colour_by = "group")
# a group difference is not clear: it was introduced by the simulation

library(scater)
assays(sim_DGE)$counts = assays(sim_DGE)$spliced
sim_DGE <- computeLibraryFactors(sim_DGE)
sim_DGE <- logNormCounts(sim_DGE)
sim_DGE <- runUMAP(sim_DGE)

sim_DGE$sample_id = factor(sim_DGE$sample_id)

plotUMAP(sim_DGE, colour_by = "cell_type")
plotUMAP(sim_DGE, colour_by = "sample_id")
plotUMAP(sim_DGE, colour_by = "group")


