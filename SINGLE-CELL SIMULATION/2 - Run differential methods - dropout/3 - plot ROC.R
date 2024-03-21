#alias R='R_LIBS=/home/stiber/R/x86_64-pc-linux-gnu-library/4.3 /home/Shared/simone/DifferentialRegulation/software/R-4.3.0/bin/R'

#cd /home/Shared/simone/Diff_Velo/NEW_simulation

#R
rm(list = ls())

setwd("~/Desktop/Differential Velocity/FINAL SCRIPTS/SINGLE-CELL SIMULATION")

topN_ALL = ROC_ALL = list()
library(iCOBRA)
library(ggplot2)
library(SingleCellExperiment)

source("08_plots/3 - plot theme.R")
all_colours = c(all_colours, "white")

for(drop in c(90, 95, 99)){
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  # load ground truth:
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  #if(drop == 90){
  sim = readRDS("03_data/mouse_simulation_data-ALEVIN-dropout90.rds")
  #}
  #if(drop == 95){
  #  sim = readRDS("03_data/mouse_simulation_data-ALEVIN-dropout95.rds")
  #}
  #if(drop == 99){
  #  sim = readRDS("03_data/mouse_simulation_data-ALEVIN-dropout99.rds")
  #}
  
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
    tot_counts_A = rowSums(as.matrix(assays(sim_one_cl)$TOT_counts[,sel_A]))
    tot_counts_B = rowSums(as.matrix(assays(sim_one_cl)$TOT_counts[,!sel_A]))
    
    tot_S_A = rowSums(as.matrix(assays(sim_one_cl)$spliced[,sel_A]))
    tot_S_B = rowSums(as.matrix(assays(sim_one_cl)$spliced[,!sel_A]))
    
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
  if(drop == 90){
    save_name = paste0("07_results/DifferentialRegulation_drop90.RData")
  }
  if(drop == 95){
    save_name = paste0("07_results/DifferentialRegulation_drop95.RData")
  }
  if(drop == 99){
    save_name = paste0("07_results/DifferentialRegulation_drop99.RData")
  }
  
  load(save_name)
  
  results_EC = lapply(results_EC, function(X) X[[1]])
  results_EC = do.call(rbind, results_EC)
  
  colnames(results_EC)[1:2] = c("gene_id", "cluster_id")
  colnames(results_EC)[3:6] = paste0("DifferentialRegulation_", colnames(results_EC)[3:6])
  
  results_EC[,6] = sapply(results_EC[,6], function(x){
    2*min(x, 1-x)
  })
  
  DF_merged = merge(DF_merged, results_EC, by = c("gene_id", "cluster_id"), all.x = TRUE)
  head(DF_merged)
  
  rm(results_EC);
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  # eisaR:
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  if(drop == 90){
    save_name = paste0("07_results/results_drop90_eisaR.RData")
  }
  if(drop == 95){
    save_name = paste0("07_results/results_drop95_eisaR.RData")
  }
  if(drop == 99){
    save_name = paste0("07_results/results_drop99_eisaR.RData")
  }
  load(save_name)
  
  colnames(DF_eisaR) = c("gene_id", "cluster_id",
                         "eisaR", "eisaR_FDR")
  
  DF_merged = merge(DF_merged, DF_eisaR, by = c("gene_id", "cluster_id"), all.x = TRUE)
  head(DF_merged)
  
  rm(DF_eisaR);
  
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  # DEXSeq - USA:
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  if(drop == 90){
    save_name = paste0("07_results/results_drop90_DEXSeq_USA.RData")
  }
  if(drop == 95){
    save_name = paste0("07_results/results_drop95_DEXSeq_USA.RData")
  }
  if(drop == 99){
    save_name = paste0("07_results/results_drop99_DEXSeq_USA.RData")
  }
  load(save_name)
  
  colnames(DF_DEXSeq) = c("DEXSeq_USA", "gene_id", "cluster_id")
  
  DF_merged = merge(DF_merged, DF_DEXSeq, 
                    by = c("gene_id", "cluster_id"),
                    all.x = TRUE)
  head(DF_merged)
  
  rm(DF_DEXSeq);
  
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  # DRIMSeq - USA:
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  if(drop == 90){
    save_name = paste0("07_results/results_drop90_DRIMSeq.RData")
  }
  if(drop == 95){
    save_name = paste0("07_results/results_drop95_DRIMSeq.RData")
  }
  if(drop == 99){
    save_name = paste0("07_results/results_drop99_DRIMSeq.RData")
  }
  load(save_name)
  
  colnames(DF_DRIMSeq) = c("gene_id", "cluster_id",
                           "DRIMSeq", "DRIMSeq_FDR")
  
  DF_merged = merge(DF_merged, DF_DRIMSeq, 
                    by = c("gene_id", "cluster_id"),
                    all.x = TRUE)
  head(DF_merged)
  
  rm(DF_DRIMSeq);
  
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  # BRIE2 50:50 without sample_id
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  if(drop == 90){
    load("07_results/BRIE2_USA_drop90.RData")
  }
  if(drop == 95){
    load("07_results/BRIE2_USA_drop95.RData")
  }
  if(drop == 99){
    load("07_results/BRIE2_USA_drop99.RData")
  }
  
  res = data.frame(gene_id = RES$Gene_id,
                   cluster_id =  RES$Cell_type,
                   BRIE2 = RES$p_BRIE2,
                   BRIE2_FDR = RES$p_BRIE2_adj)
  
  DF_merged = merge(DF_merged, res, by = c("gene_id", "cluster_id"), all.x = TRUE)
  head(DF_merged)
  
  rm(RES)
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  # satuRn:
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  if(drop == 90){
    load("07_results/results_drop90_SatuRn.RData")
  }
  if(drop == 95){
    load("07_results/results_drop95_SatuRn.RData")
  }
  if(drop == 99){
    load("07_results/results_drop99_SatuRn.RData")
  }
  
  clusters = c("Adipocytes", "Epithelial_cells", "Hepatocytes")
  
  satuRn = list()
  for(cl in 1:3){
    res = DF_SatuRn[[cl]]
    tr_names = res$transcript_id
    tr_names = strsplit(tr_names, "-")
    tr_names = sapply(tr_names, function(x){ x[[1]]})
    
    res$empirical_pval[is.na(res$empirical_pval)] = 1
    res_by_tr = split(res$empirical_pval, tr_names)
    min_p_val = sapply(res_by_tr, min, na.rm = TRUE)
    
    satuRn[[cl]] = data.frame(gene_id = names(min_p_val),
                              cluster_id =  clusters[cl],
                              satuRn = min_p_val) 
  }
  satuRn = do.call(rbind, satuRn)
  
  DF_merged = merge(DF_merged, satuRn, by = c("gene_id", "cluster_id"), all.x = TRUE)
  head(DF_merged)
  tail(DF_merged)
  
  rm(satuRn); rm(DF_SatuRn);
  rm(min_p_val); rm(tr_names); rm(res)
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  # satuRn SC:
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  if(drop == 90){
    load("07_results/results_drop90_SatuRnSCsampleDesign.RData")
  }
  if(drop == 95){
    load("07_results/results_drop95_SatuRnSCsampleDesign.RData")
  }
  if(drop == 99){
    load("07_results/results_drop99_SatuRnSCsampleDesign.RData")
  }
  
  clusters = c("Adipocytes", "Epithelial_cells", "Hepatocytes")
  
  satuRn = list()
  for(cl in 1:3){
    res = DF_SatuRn[[cl]]
    tr_names = res$transcript_id
    tr_names = strsplit(tr_names, "-")
    tr_names = sapply(tr_names, function(x){ x[[1]]})
    
    res$empirical_pval[is.na(res$empirical_pval)] = 1
    res_by_tr = split(res$empirical_pval, tr_names)
    min_p_val = sapply(res_by_tr, min, na.rm = TRUE)
    
    satuRn[[cl]] = data.frame(gene_id = names(min_p_val),
                              cluster_id =  clusters[cl],
                              satuRn_SC = min_p_val) 
  }
  satuRn = do.call(rbind, satuRn)
  
  DF_merged = merge(DF_merged, satuRn, by = c("gene_id", "cluster_id"), all.x = TRUE)
  head(DF_merged)
  tail(DF_merged)
  
  rm(satuRn); rm(DF_SatuRn);
  rm(min_p_val); rm(tr_names); rm(res)
  
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  # iCOBRA:
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  DF_merged[is.na(DF_merged)] = 1
  
  pval = data.frame(
    BRIE2 = DF_merged$BRIE2,
    DEXSeq_TECs = DF_merged$DEXSeq_USA,
    DifferentialRegulation = DF_merged$`DifferentialRegulation_Prob-B-UP`,
    DifferentialRegulation_Wald = DF_merged$DifferentialRegulation_p_val,
    DRIMSeq = DF_merged$DRIMSeq,
    eisaR = DF_merged$eisaR,
    satuRn = DF_merged$satuRn,
    satuRn_SC = DF_merged$satuRn_SC
  )
  
  DF_COBRA <- COBRAData(
    pval = pval,
    truth = data.frame(status = DF_merged$truth))
  
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  # ROC, FDR plots:
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  perf <- calculate_performance(DF_COBRA, binary_truth = "status", 
                                aspects = c("roc", "fpc"),
                                thrs = c(0.01, 0.05, 0.1, 0.2))
  
  cobra_plot <- prepare_data_for_plot(perf, colorscheme = all_colours, incloverall = FALSE,
                                      facetted = TRUE,conditionalfill = FALSE,
                                      keepmethods = methods_all)
  
  # plot ROC curve
  ROC = plot_roc(cobra_plot,linewidth=2) + 
    scale_color_manual(values = all_colours,
                       name = "",
                       breaks=methods_all,
                       labels=methods_all)+
    my_theme +
    geom_abline(slope = 1, intercept = 0 )
  
  # top N results:
  TOP = plot_fpc(cobra_plot, maxnfdc = 3000, linewidth = 2 ) +
    theme(legend.position="bottom") + my_theme +
    scale_color_manual(values = all_colours,
                       name = "",
                       breaks=methods_all,
                       labels=methods_all) +
    coord_cartesian(
      ylim = c(0,200),
      expand = TRUE,
      default = FALSE,
      clip = "on")
  
  if(drop == 90){
    ROC_ALL[[1]] = ROC
    topN_ALL[[1]] = TOP
  }
  if(drop == 95){
    ROC_ALL[[2]] = ROC
    topN_ALL[[2]] = TOP
  }
  if(drop == 99){
    ROC_ALL[[3]] = ROC
    topN_ALL[[3]] = TOP
  }
}

my_theme2 = theme(legend.position = "none", 
                  text = element_text(size = 9),
                  plot.title = element_text(hjust = 0, face = "bold", size=35),
                  plot.subtitle = element_text(hjust = 0.5, face = "bold", size=30))

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# 3 ROC and 2 topN together:
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
GG_ROC_topN = ggpubr::ggarrange(ROC_ALL[[1]] + labs(subtitle = "zeros 90") + my_theme2,
                                ROC_ALL[[2]] + ylab("") + labs(subtitle = "zeros 95") + my_theme2,
                                ROC_ALL[[3]] + ylab("") + labs(subtitle = "zeros 99") + my_theme2,
                                topN_ALL[[1]] + my_theme2,
                                topN_ALL[[2]] + ylab("") + my_theme2,
                                topN_ALL[[3]] + ylab("") + my_theme2,
                                legend = "bottom", 
                                common.legend = TRUE,
                                ncol = 3, nrow = 2)

save_path = file.path("08_plots")
ggsave(filename = paste0('DGE_drop_sc.pdf'),
       plot = GG_ROC_topN,
       device = "pdf",
       path = save_path,
       width = 30,
       height = 22,
       units = "in",
       dpi = 300,
       limitsize = TRUE)
