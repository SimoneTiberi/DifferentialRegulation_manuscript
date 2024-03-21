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

AUC = list()

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
  # eisaR:
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  if(DGE){
    save_name = paste0("07_results/results_DGE_eisaR.RData")
  }else{
    save_name = paste0("07_results/results_NO_DGE_eisaR.RData")
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
  if(DGE){
    save_name = paste0("07_results/results_DGE_DEXSeq_USA.RData")
  }else{
    save_name = paste0("07_results/results_NO_DGE_DEXSeq_USA.RData")
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
  if(DGE){
    save_name = paste0("07_results/results_DGE_DRIMSeq.RData")
  }else{
    save_name = paste0("07_results/results_NO_DGE_DRIMSeq.RData")
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
  # satuRn:
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  if(DGE){
    load("07_results/results_DGE_SatuRn.RData")
  }else{
    load("07_results/results_NO_DGE_SatuRn.RData")
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
    
    qval = p.adjust(min_p_val, method = "BH")
    
    satuRn[[cl]] = data.frame(gene_id = names(min_p_val),
                              cluster_id =  clusters[cl],
                              satuRn = min_p_val,
                              satuRn_FDR = qval) 
  }
  satuRn = do.call(rbind, satuRn)
  
  DF_merged = merge(DF_merged, satuRn, by = c("gene_id", "cluster_id"), all.x = TRUE)
  head(DF_merged)
  tail(DF_merged)
  
  rm(satuRn); rm(DF_SatuRn); rm(qval);
  rm(min_p_val); rm(tr_names); rm(res)
  
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  # satuRn SINGLE-CELL DATA
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  if(DGE){
    load("07_results/results_DGE_SatuRnSCsampleDesign.RData")
  }else{
    load("07_results/results_NO_DGE_SatuRnSCsampleDesign.RData")
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
    
    qval = p.adjust(min_p_val, method = "BH")
    
    satuRn[[cl]] = data.frame(gene_id = names(min_p_val),
                              cluster_id =  clusters[cl],
                              satuRn_SC = min_p_val,
                              satuRn_SC_FDR = qval) 
  }
  satuRn = do.call(rbind, satuRn)
  
  DF_merged = merge(DF_merged, satuRn, by = c("gene_id", "cluster_id"), all.x = TRUE)
  head(DF_merged)
  tail(DF_merged)
  
  rm(satuRn); rm(DF_SatuRn); rm(qval);
  rm(min_p_val); rm(tr_names); rm(res)
  
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  # iCOBRA:
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  DF_merged[is.na(DF_merged)] = 1
  
  DF_merged$TOT_counts = DF_merged$tot_counts_A + DF_merged$tot_counts_B
  
  summary(DF_merged$TOT_counts)
  
  sel_rows = list()
  
  # all results:
  sel_rows[[1]] = rep(TRUE, length(DF_merged) )
  
  # split in low, medium and high counts:
  cuts = quantile(DF_merged$TOT_counts, c(1/3, 2/3))
  
  sel_rows[[2]] = DF_merged$TOT_counts <= cuts[1]
  sel_rows[[3]] = (DF_merged$TOT_counts > cuts[1]) & (DF_merged$TOT_counts <= cuts[2])
  sel_rows[[4]] = DF_merged$TOT_counts > cuts[2]
  
  ROC = TOP = list()
  
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
  
  # AUC calculation:
  library(pROC)
  if(DGE){
    AUC[[2]] = apply(pval, 2, function(x){
      auc(DF_merged$truth, x)[1]
    })
  }else{
    AUC[[1]] = apply(pval, 2, function(x){
      auc(DF_merged$truth, x)[1]
    })
  }
  #}
  
  for(sel in 1:4){
    DF_sel = DF_merged[sel_rows[[sel]], ]
    
    DF_COBRA <- COBRAData(
      pval = data.frame(
        BRIE2 = DF_sel$BRIE2,
        DEXSeq_TECs = DF_sel$DEXSeq_USA,
        DifferentialRegulation = DF_sel$`DifferentialRegulation_Prob-B-UP`,
        DifferentialRegulation_Wald = DF_sel$DifferentialRegulation_p_val,
        DRIMSeq = DF_sel$DRIMSeq,
        eisaR = DF_sel$eisaR,
        satuRn = DF_sel$satuRn,
        satuRn_SC = DF_sel$satuRn_SC
      ),
      truth = data.frame(status = DF_sel$truth))
    
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
    ROC[[sel]] = plot_roc(cobra_plot,linewidth=2) + 
      scale_color_manual(values = all_colours,
                         name = "",
                         breaks=methods_all,
                         labels=methods_all)+
      my_theme +
      geom_abline(slope = 1, intercept = 0 )
    
    # top N results:
    TOP[[sel]] = plot_fpc(cobra_plot, maxnfdc = ifelse(sel == 1, 3000, 1000), linewidth = 2 ) +
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
    
    # limit y axis (BRIE2 has too many FPs)
  }
  
  print(DGE)
  if(DGE){
    ROC_ALL[[2]] = ROC
    topN_ALL[[2]] = TOP
  }else{
    ROC_ALL[[1]] = ROC
    topN_ALL[[1]] = TOP
  }
  
  if(DGE){
    save_path = file.path("08_plots/DGE")
  }else{
    save_path = file.path("08_plots/non-DGE")
  }
  
  ggsave(filename = paste0('ROC_sc.pdf'),
         plot = ROC[[1]],
         device = "pdf",
         path = save_path,
         width = 25,
         height = 20,
         units = "in",
         dpi = 300,
         limitsize = TRUE)
  
  
  ggsave(filename = paste0('TOP_sc.pdf'),
         plot = TOP[[1]],
         device = "pdf",
         path = save_path,
         width = 20,
         height = 15,
         units = "in",
         dpi = 300,
         limitsize = TRUE)
  
  # make 1 plot (1 row and 2 columns) with ROC[[1]] + topN[[1]]
  # 1 unique legend at the bottom
  ROC_topN <- ggpubr::ggarrange(ROC[[1]], TOP[[1]], ncol=2, nrow=1, common.legend = TRUE, legend="bottom")
  ggsave(filename = paste0('ROC_topN_sc.pdf'),
         plot = ROC_topN,
         device = "pdf",
         path = save_path,
         width = 25,
         height = 15,
         units = "in",
         dpi = 300,
         limitsize = TRUE)
  
  ggsave(filename = paste0('ROC_topN_sc2.pdf'),
         plot = ROC_topN,
         device = "pdf",
         path = save_path,
         width = 35,
         height = 20,
         units = "in",
         dpi = 300,
         limitsize = TRUE)
}

my_theme2 = theme(legend.position = "none", 
                  text = element_text(size = 9),
                  plot.title = element_text(hjust = 0, face = "bold", size=35),
                  plot.subtitle = element_text(hjust = 0.5, face = "bold", size=30))

#### #### #### #### #### #### #### 
# 2 ROC curves
#### #### #### #### #### #### #### 
GG_ROC = ggpubr::ggarrange(ROC_ALL[[1]][[1]] + labs(subtitle = "DR") + my_theme2,
                           ROC_ALL[[2]][[1]] + ylab("") + labs(subtitle = "DR + DGE") + my_theme2,
                           legend = "bottom", 
                           common.legend = TRUE,
                           ncol = 2, nrow = 1)

save_path = file.path("08_plots")
ggsave(filename = paste0('2_ROC_sc.pdf'),
       plot = GG_ROC,
       device = "pdf",
       path = save_path,
       width = 20,
       height = 11,
       units = "in",
       dpi = 300,
       limitsize = TRUE)

#### #### #### #### #### #### #### 
# 2 topN curves
#### #### #### #### #### #### #### 
GG_topN = ggpubr::ggarrange(topN_ALL[[1]][[1]] + labs(subtitle = "DR") + my_theme2,
                            topN_ALL[[2]][[1]] + ylab("") + labs(subtitle = "DR + DGE") + my_theme2,
                            legend = "bottom", 
                            common.legend = TRUE,
                            ncol = 2, nrow = 1)

save_path = file.path("08_plots")
ggsave(filename = paste0('2_topN_sc.pdf'),
       plot = GG_topN,
       device = "pdf",
       path = save_path,
       width = 20,
       height = 11,
       units = "in",
       dpi = 300,
       limitsize = TRUE)

#### #### #### #### #### #### #### 
# 6 ROC curves (DGE and NON-DGE for 3 abundance levels)
#### #### #### #### #### #### #### 
AA1 = ggpubr::ggarrange(ROC_ALL[[1]][[2]] + xlab("") + labs(title = "Low abundance", subtitle = "DR") + my_theme2,
                        ROC_ALL[[2]][[2]] + xlab("") + ylab("") + labs(subtitle = "DR + DGE") + my_theme2,
                        ncol = 2, nrow = 1)

AA2 = ggpubr::ggarrange(ROC_ALL[[1]][[3]] + xlab("") + labs(title = "Mid abundance") + my_theme2,
                        ROC_ALL[[2]][[3]] + xlab("") + ylab("")  + my_theme2,
                        ncol = 2, nrow = 1)

AA3 = ggpubr::ggarrange(ROC_ALL[[1]][[4]] + labs(title = "High abundance") + my_theme2,
                        ROC_ALL[[2]][[4]] + ylab("") + my_theme2,
                        ncol = 2, nrow = 1)

legend <- ggpubr::get_legend(ROC_ALL[[1]][[2]] +
                               #guides(colour = guide_legend(nrow = 1, byrow = FALSE,
                               #                             override.aes = list(shape = shape_fill,
                               #                                                 fill = all_colours) ) ) +
                               theme(legend.key.height=unit(3,"line"),
                                     legend.text=element_text(size=35),
                                     legend.key.width=unit(6,"line")))

ggpubr::as_ggplot(legend)

AA = ggpubr::ggarrange( AA1, AA2, AA3,legend,
                        heights=c(1.6,1.6,1.6,0.5),
                        widths = c(1,1,1,1),align = "v",
                        ncol = 1, nrow = 4 )
save_path = file.path("08_plots")
ggsave(filename = paste0('Sup_ROC_sc_2.pdf'),
       plot = AA,
       device = "pdf",
       path = save_path,
       width = 30,
       height = 40,
       units = "in",
       dpi = 300,
       limitsize = TRUE)

#### #### #### #### #### #### #### 
# 6 topN curves (DGE and NON-DGE for 3 abundance levels)
#### #### #### #### #### #### #### 
AA1 = ggpubr::ggarrange(topN_ALL[[1]][[2]] + xlab("") + labs(title = "Low abundance", subtitle = "DR") + my_theme2,
                        topN_ALL[[2]][[2]] + xlab("") + ylab("") + labs(subtitle = "DR + DGE") + my_theme2,
                        ncol = 2, nrow = 1)

AA2 = ggpubr::ggarrange(topN_ALL[[1]][[3]] + xlab("") + labs(title = "Mid abundance") + my_theme2,
                        topN_ALL[[2]][[3]] + xlab("") + ylab("")  + my_theme2,
                        ncol = 2, nrow = 1)

AA3 = ggpubr::ggarrange(topN_ALL[[1]][[4]] + labs(title = "High abundance") + my_theme2,
                        topN_ALL[[2]][[4]] + ylab("") + my_theme2,
                        ncol = 2, nrow = 1)

legend <- ggpubr::get_legend(topN_ALL[[1]][[2]] +
                               #guides(colour = guide_legend(nrow = 1, byrow = FALSE,
                               #                             override.aes = list(shape = shape_fill,
                               #                                                 fill = all_colours) ) ) +
                               theme(legend.key.height=unit(3,"line"),
                                     legend.text=element_text(size=35),
                                     legend.key.width=unit(6,"line")))

ggpubr::as_ggplot(legend)

AA = ggpubr::ggarrange( AA1, AA2, AA3,legend,
                        heights=c(1.6,1.6,1.6,0.5),
                        widths = c(1,1,1,1),align = "v",
                        ncol = 1, nrow = 4 )

save_path = file.path("08_plots")
ggsave(filename = paste0('Sup_topN_sc_2.pdf'),
       plot = AA,
       device = "pdf",
       path = save_path,
       width = 30,
       height = 40,
       units = "in",
       dpi = 300,
       limitsize = TRUE)


#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# 2 ROC and 2 topN together:
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
GG_ROC_topN = ggpubr::ggarrange(ROC_ALL[[1]][[1]] + labs(subtitle = "DR") + my_theme2,
                                ROC_ALL[[2]][[1]] + ylab("") + labs(subtitle = "DR + DGE") + my_theme2,
                                topN_ALL[[1]][[1]] + labs(subtitle = "DR") + my_theme2,
                                topN_ALL[[2]][[1]] + ylab("") + labs(subtitle = "DR + DGE") + my_theme2,
                                legend = "bottom", 
                                common.legend = TRUE,
                                ncol = 2, nrow = 2)

ggsave(filename = paste0('2_ROC_2_topN_sc.pdf'),
       plot = GG_ROC_topN,
       device = "pdf",
       path = save_path,
       width = 20,
       height = 22,
       units = "in",
       dpi = 300,
       limitsize = TRUE)


# AUC
A = do.call(rbind, AUC)
overall = colMeans(A)

A = A[,order(overall)]
# library(xtable)
# xtable(A, digits =3)
# \begin{table}[ht]
# \centering
# \begin{tabular}{rrrrrrrrr}
# \hline
# & BRIE2 & DRIMSeq & satuRn\_SC & DifferentialRegulation & DifferentialRegulation\_Wald & satuRn & eisaR & DEXSeq\_TECs \\ 
# \hline
# 1 & 0.568 & 0.593 & 0.789 & 0.770 & 0.770 & 0.777 & 0.774 & 0.791 \\ 
# 2 & 0.563 & 0.597 & 0.710 & 0.763 & 0.764 & 0.764 & 0.768 & 0.786 \\ 
# \hline
# \end{tabular}
# \end{table}
# 
# library(xtable)
# xtable(A, digits =2)
# \begin{table}[ht]
# \centering
# \begin{tabular}{rrrrrrrrr}
# \hline
# & BRIE2 & DRIMSeq & satuRn\_SC & DifferentialRegulation & DifferentialRegulation\_Wald & satuRn & eisaR & DEXSeq\_TECs \\ 
# \hline
# 1 & 0.57 & 0.59 & 0.74 & 0.77 & 0.77 & 0.78 & 0.77 & 0.79 \\ 
# 2 & 0.56 & 0.60 & 0.68 & 0.76 & 0.76 & 0.76 & 0.77 & 0.79 \\ 
# \hline
# \end{tabular}
# \end{table}

