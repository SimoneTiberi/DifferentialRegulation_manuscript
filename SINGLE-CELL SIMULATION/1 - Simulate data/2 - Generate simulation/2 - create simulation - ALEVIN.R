# cd /home/Shared_sherborne/simone/Diff_Velo/NEW_simulation

setwd("~/Desktop/Differential Velocity/FINAL SCRIPTS/SINGLE-CELL SIMULATION")
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
source("1 - Simulate data/2 - Generate simulation/functions.R")

##########################################
# simulate mouse data based on preset parameters
##########################################
set.seed(123)
sim <- mouse_simulation(sce = sce, CLUSTERS = CLUSTERS, 
                        p_genes = p_genes)

# save simulated data sets
saveRDS(sim, file = "03_data/mouse_simulation_data-ALEVIN.rds")

##########################################
# batch simulation
##########################################
batches = c("1", "1", "2", "2")
table(GROUP, batches)

sim$batch <- ifelse(sim$sample_id %in% c("1", "2"), "1", "2")
sim$batch
table(sim$batch, sim$group)

table(sim$sample_id)

set.seed(12)
sim_batch <- add_batch(sce = sim, 
                       CLUSTERS = CLUSTERS)

table(metadata(sim_batch)$truth$truth); table(metadata(sim_batch)$truth$truth_batch)
sum(metadata(sim_batch)$truth$truth); sum(metadata(sim_batch)$truth$truth_batch)

# save simulated data sets
saveRDS(sim_batch, file = "03_data/mouse_simulation_batch_data-ALEVIN.rds")

##########################################
# NULL SIMULATION
##########################################
set.seed(1)
sim_null <- mouse_simulation(sce = sce, CLUSTERS = CLUSTERS, 
                             p_genes = 0)

sim_null

# save simulated data sets
saveRDS(sim_null, file = "03_data/mouse_simulation_data-NULL-ALEVIN.rds")

##########################################
# DGE simulation
##########################################
set.seed(1234)
sim_2 <- mouse_simulation(sce = sce, CLUSTERS = CLUSTERS, 
                          p_genes = p_genes)
set.seed(12345)
sim_DGE <- add_DGE(sce = sim_2, CLUSTERS = CLUSTERS, 
                   p_genes_DGE = p_genes_DGE,
                   mu_FC = 3)

# save simulated data sets
saveRDS(sim_DGE, file = "03_data/mouse_simulation_DGE_data-ALEVIN.rds")

##########################################
# DGE simulation - FC 6
##########################################
set.seed(123456)
sim_3 <- mouse_simulation(sce = sce, CLUSTERS = CLUSTERS, 
                          p_genes = p_genes)
set.seed(1234567)
sim_DGE_FC6 <- add_DGE(sce = sim_3, CLUSTERS = CLUSTERS, 
                       p_genes_DGE = p_genes_DGE,
                       mu_FC = 6)

# save simulated data sets
saveRDS(sim_DGE_FC6, file = "03_data/mouse_simulation_DGE_FC6_data-ALEVIN.rds")

##########################################
# DGE simulation - FC 9
##########################################
set.seed(12345678)
sim_4 <- mouse_simulation(sce = sce, CLUSTERS = CLUSTERS, 
                          p_genes = p_genes)
set.seed(123456789)
sim_DGE_FC9 <- add_DGE(sce = sim_4, CLUSTERS = CLUSTERS, 
                       p_genes_DGE = p_genes_DGE,
                       mu_FC = 9)

# save simulated data sets
saveRDS(sim_DGE_FC9, file = "03_data/mouse_simulation_DGE_FC9_data-ALEVIN.rds")

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
plotUMAP(sim, colour_by = "batch")
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


