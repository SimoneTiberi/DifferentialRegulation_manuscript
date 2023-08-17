alias R='R_LIBS=/home/stiber/R/x86_64-pc-linux-gnu-library/4.3 /home/Shared/simone/DifferentialRegulation/software/R-4.3.0/bin/R'

cd /home/Shared/simone/Diff_Velo/NEW_simulation

R
rm(list = ls())
setwd("~/Desktop/Differential Velocity/FINAL SCRIPTS/SINGLE-CELL SIMULATION")

ROC_ALL = list()
library(iCOBRA)
library(ggplot2)
library(SingleCellExperiment)

source("08_plots/3 - plot theme.R")
for(DGE in c(FALSE, TRUE)){
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  # load ground truth:
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  if(DGE){
    sim = readRDS("03_data/mouse_simulation_DGE_data-ALEVIN.rds")
  }else{
    sim = readRDS("03_data/mouse_simulation_data-ALEVIN.rds")
  }
  DF_merged = metadata(sim)$truth
  
  DF_merged$gene_id =  DF_merged$Gene_id
  DF_merged$cluster_id =  DF_merged$Cell_type
  
  DF_merged = DF_merged[,-c(1:2)]
  
  DF = list()
  
  clusters = unique(sim$cell_type)
  for(cl in 1:length(clusters)){
    sel_cl = sim$cell_type == clusters[cl]
    sim_one_cl = sim[,sel_cl]
    
    sel_A = sim_one_cl$group == "A"
    tot_counts_A = rowSums(assays(sim_one_cl)$TOT_counts[,sel_A])
    tot_counts_B = rowSums(assays(sim_one_cl)$TOT_counts[,!sel_A])
    
    tot_S_A = rowSums(assays(sim_one_cl)$spliced[,sel_A])
    tot_S_B = rowSums(assays(sim_one_cl)$spliced[,!sel_A])
    
    DF[[cl]] = data.frame(gene_id = rownames(sim_one_cl),
                          cluster_id = clusters[cl],
                          tot_counts_A = tot_counts_A,
                          tot_counts_B = tot_counts_B,
                          pi_S_A = tot_S_A/tot_counts_A,
                          pi_S_B = tot_S_B/tot_counts_B)
  }
  DF = do.call(rbind, DF)
  
  DF_merged = merge(DF_merged, DF, by = c("gene_id", "cluster_id"))
  head(DF_merged)
  
  # keep results with at least 20 counts per group:
  sel_counts = (DF_merged$tot_counts_A >= 10) &  (DF_merged$tot_counts_B >= 10)
  
  DF_merged = DF_merged[sel_counts,]
  
  rm(sim_one_cl); rm(DF); 
  
  DF_merged$cluster_id[DF_merged$cluster_id == "Epithelial cells"]= "Epithelial_cells"
  table(DF_merged$cluster_id)
  
  table(DF_merged$truth)
  
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  # DifferentialRegulation
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  if(DGE){
    save_name = paste0("07_results/DifferentialRegulation_DGE_2000_10.RData")
  }else{
    save_name = paste0("07_results/DifferentialRegulation_NO_DGE_2000_10.RData")
  }
  load(save_name)
  
  results_EC = lapply(results_EC, function(X) X[[1]])
  results_EC = do.call(rbind, results_EC)
  
  colnames(results_EC)[1:2] = c("gene_id", "cluster_id")
  colnames(results_EC)[3:6] = paste0("DifferentialRegulation_", colnames(results_EC)[3:6])
  
  results_EC[,6] = sapply(results_EC[,6], function(x){
    2*min(x, 1-x)
  })
  
  DF_merged = merge(DF_merged, results_EC, by = c("gene_id", "cluster_id"))
  head(DF_merged)
  
  rm(results_EC);
  
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  # BRIE2 USA without sample_id
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  if(DGE){
    load("07_results/BRIE2_USA_DGE.RData")
    
    res = data.frame(gene_id = RES$Gene_id,
                     cluster_id =  RES$Cell_type,
                     BRIE2 = RES$p_BRIE2,
                     BRIE2_FDR = RES$p_BRIE2_adj)
  }else{
    load("07_results/BRIE2_USA.RData")
    
    res = data.frame(gene_id = RES$Gene_id,
                     cluster_id =  RES$Cell_type,
                     BRIE2 = RES$p_BRIE2,
                     BRIE2_FDR = RES$p_BRIE2_adj)
  }
  
  DF_merged = merge(DF_merged, res, by = c("gene_id", "cluster_id"), all.x = TRUE)
  head(DF_merged)
  
  rm(RES)
  
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  # BRIE2 USA WITH sample_id
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  if(DGE){
    load("07_results/BRIE2_USA_sample_DGE.RData")
    
    res = data.frame(gene_id = RES$Gene_id,
                     cluster_id =  RES$Cell_type,
                     BRIE2_sample = RES$p_BRIE2,
                     BRIE2_sample_FDR = RES$p_BRIE2_adj)
  }else{
    load("07_results/BRIE2_USA_sample.RData")
    
    res = data.frame(gene_id = RES$Gene_id,
                     cluster_id =  RES$Cell_type,
                     BRIE2_sample = RES$p_BRIE2,
                     BRIE2_sample_FDR = RES$p_BRIE2_adj)
  }
  
  DF_merged = merge(DF_merged, res, by = c("gene_id", "cluster_id"), all.x = TRUE)
  head(DF_merged)
  
  rm(RES)
  
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  # iCOBRA:
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  colSums(is.na(DF_merged))
  
  DF_merged[is.na(DF_merged)] = 1
  
  DF_merged$TOT_counts = DF_merged$tot_counts_A + DF_merged$tot_counts_B
  
  summary(DF_merged$TOT_counts)

  DF_COBRA <- COBRAData(
    pval = data.frame(
      BRIE2 = DF_merged$BRIE2,
      BRIE2_sample = DF_merged$BRIE2_sample
    ),
    truth = data.frame(status = DF_merged$truth))
  
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  # ROC, FDR plots:
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  perf <- calculate_performance(DF_COBRA, binary_truth = "status",
                                aspects = c("roc", "fpc"),
                                thrs = c(0.01, 0.05, 0.1, 0.2))
  
  cobra_plot <- prepare_data_for_plot(perf, 
                                      colorscheme = c(all_colours[1], "#E31A1C"), 
                                      incloverall = FALSE,
                                      facetted = TRUE,
                                      conditionalfill = FALSE)
  
  # plot ROC curve
  ROC = plot_roc(cobra_plot,linewidth=2) + 
    my_theme +
    geom_abline(slope = 1, intercept = 0 )
  
  print(DGE)
  if(DGE){
    ROC_ALL[[2]] = ROC
  }else{
    ROC_ALL[[1]] = ROC
  }
}

# make an overall plot of the 6 images in 
# ROC_ALL[[1]][2:4] # NON-DGE
# ROC_ALL[[2]][2:4] # DGE
# 1 unique legend at the bottom
# 3 rows = low, mid and high abundance (sel)
# 2 cols = type (DR, DR + DGE)
my_theme2 = theme(text = element_text(size = 9),
                  plot.title = element_text(hjust = 0, face = "bold", size=35),
                  plot.subtitle = element_text(hjust = 0.5, face = "bold", size=30))

AA = ggpubr::ggarrange(ROC_ALL[[1]] + xlab("") + 
                          labs(subtitle = "DR") + my_theme2,
                       ROC_ALL[[2]] + xlab("") + 
                          ylab("") + 
                          labs(subtitle = "DR + DGE") + my_theme2,
                        legend = "bottom",
                       common.legend = TRUE,
                        ncol = 2, nrow = 1)

save_path = file.path("08_plots")
ggsave(filename = paste0('ROC_BRIE2.pdf'),
       plot = AA,
       device = "pdf",
       path = save_path,
       width = 30,
       height = 16,
       units = "in",
       dpi = 300,
       limitsize = TRUE)
