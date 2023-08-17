alias R='R_LIBS=/home/stiber/R/4.3.1 /usr/local/R/R-4.3.1/bin/R'
cd /home/Shared/simone/DifferentialRegulation

R

rm(list = ls())
setwd("~/Desktop/Differential Velocity/FINAL SCRIPTS/BULK SIMULATION")
#setwd("/home/Shared/simone/DifferentialRegulation")

library(data.table)
library(iCOBRA)
library(ggplot2)

source("2 - R scripts/2 - plot theme.R")
all_colours = c(all_colours, "white")

############################################################
# directories:
############################################################
directories = c("de0_ds0_dr0",
                "de0_ds0_dr2000",
                "de2000_ds0_dr2000",
                "de0_ds2000_dr2000")

directories_results = c("NULL",
                        "DR",
                        "DR + DGE",
                        "DR + DAS")
############################################################
# specify type of data:
############################################################
# save the 3 plots below (ROC, FDR and top_N)
plot_ROC_list <- list()
plot_topN_list <- list()

ROC_type = topN_type = list()

AUC = list()

for(type in 2:4){
  # null
  # DR
  # DR + DGE
  # DR + DAS
  
  ############################################################
  # load ground truth:
  ############################################################
  data_dir = file.path("1 - data", directories[type],
                       "3_truth/simulation_details.txt")
  
  sim = fread(data_dir)
  
  DF_merged = data.frame(Transcript_id = sim$transcript_id,
                         truth = sim$transcript_dr_status,
                         counts = sim$expected_count)
  
  data_input = file.path("3 - results", directories_results[type])
  
  ############################################################
  # load results:
  ############################################################
  # DifferentialRegulation:
  name = paste("DifferentialRegulation.RData")
  file = file.path(data_input,name)
  if(file.exists(file)){
    load(file)
    DR = res_EC$Differential_results
    names(DR)[2:4]  = c("p_val_DifferentialRegulation", "p_adj_DifferentialRegulation",
                        "Prob_DifferentialRegulation")
    DR = DR[,1:4]
    
    DR[,4] = sapply(DR[,4], function(x){
      2*min(x, 1-x)
    })
    DF_merged = merge(DR, DF_merged, by = "Transcript_id")
    head(DF_merged)
    rm(DR); rm(res_EC); rm(name)
  }
  
  # eisaR:
  name = paste("eisaR.RData")
  file = file.path(data_input,name)
  if(file.exists(file)){
    load(file)
    DF_eisaR = data.frame(Transcript_id = DF_eisaR$Transcript_id,
                          p_val_eisaR = DF_eisaR$p_eisaR,
                          p_adj_eisaR = DF_eisaR$p_eisaR_adj)
    DF_merged = merge(DF_eisaR, DF_merged, by = "Transcript_id", all.y = TRUE)
    head(DF_merged)
    rm(DF_eisaR); rm(name)
  }
  
  # DEXSeq - TECs:
  name = paste("DEXSeq_TECs.RData")
  file = file.path(data_input,name)
  if(file.exists(file)){
    load(file)
    DEXSeq_TECs = data.frame(Transcript_id = res_gene$gene,
                             p_val_DEXSeq_TECs = res_gene$qval,
                             p_adj_DEXSeq_TECs = res_gene$qval)
    DF_merged = merge(DEXSeq_TECs, DF_merged, by = "Transcript_id", all.y = TRUE)
    head(DF_merged)
    rm(DEXSeq_TECs); rm(res); rm(name)
  }
  
  # DEXSeq - ECs:
  name = paste("DEXSeq_ECCs.RData")
  file = file.path(data_input,name)
  if(file.exists(file)){
    load(file)
    DEXSeq_ECs = data.frame(Transcript_id = res_gene$gene,
                            p_val_DEXSeq_ECs = res_gene$qval,
                            p_adj_DEXSeq_ECs = res_gene$qval)
    DF_merged = merge(DEXSeq_ECs, DF_merged, by = "Transcript_id", all.y = TRUE)
    head(DF_merged)
    rm(DEXSeq_ECs); rm(res_gene); rm(name)
  }
  
  # DRIMSeq:
  name = paste("DRIMSeq.RData")
  file = file.path(data_input,name)
  if(file.exists(file)){
    load(file)
    DRIMSeq = data.frame(Transcript_id = results_gene$gene_id,
                         p_val_DRIMSeq = results_gene$pvalue,
                         p_adj_DRIMSeq = results_gene$adj_pvalue)
    DF_merged = merge(DRIMSeq, DF_merged, by = "Transcript_id", all.y = TRUE)
    head(DF_merged)
    rm(DRIMSeq); rm(results_gene); rm(name)
  }
  
  # satuRn:
  name = paste("satuRn.RData")
  file = file.path(data_input,name)
  if(file.exists(file)){
    load(file)
    
    tr_names = rownames(res)
    tr_names = strsplit(tr_names, "-")
    tr_names = sapply(tr_names, function(x){ x[[1]]})
    
    res$empirical_pval[is.na(res$empirical_pval)] = 1
    res_by_tr = split(res$empirical_pval, tr_names)
    min_p_val = sapply(res_by_tr, min, na.rm = TRUE)
    
    qval = p.adjust(min_p_val, method = "BH")
    
    #theta = sort(unique(min_p_val))
    #qval = DEXSeq:::perGeneQValueExact(min_p_val, theta, res_by_tr)
    
    #res$empirical_FDR[is.na(res$empirical_FDR)] = 1
    #res$empirical_FDR[is.na(res$empirical_FDR)] = 1
    #res_by_tr = split(res$empirical_FDR, tr_names)
    #min_FDR = sapply(res_by_tr, min, na.rm = TRUE)
    
    satuRn = data.frame(Transcript_id = names(min_p_val),  #rownames(res),
                        p_val_satuRn = min_p_val, #res$empirical_pval,
                        p_adj_satuRn = qval)   #res$empirical_FDR)
    DF_merged = merge(satuRn, DF_merged, by = "Transcript_id", all.y = TRUE)
    head(DF_merged)
    rm(satuRn); rm(res); rm(name); rm(qval); #rm(theta); 
    rm(min_p_val); rm(tr_names)
  }
  
  # SUPPA2:
  name = paste("SUPPA2.RData")
  file = file.path(data_input,name)
  if(file.exists(file)){
    load(file)
    
    res_SUPPA[is.na(res_SUPPA)] = 1
    qval = p.adjust(res_SUPPA, method = "BH")
    #theta = sort(unique(res_SUPPA))
    #res_by_tr = lapply(res_SUPPA, function(i){
    #  1:2
    #})
    # qval = DEXSeq:::perGeneQValueExact(res_SUPPA, theta, res_by_tr)
    
    SUPPA2 = data.frame(Transcript_id = names(res_SUPPA),
                        p_val_SUPPA2 = res_SUPPA,
                        p_adj_SUPPA2 = qval)
    colnames(SUPPA2)[1] <- "Transcript_id"
    DF_merged = merge(SUPPA2, DF_merged, by = "Transcript_id", all.y = TRUE)
    head(DF_merged)
    rm(SUPPA2); rm(res_SUPPA); rm(name)
  }
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  # iCOBRA:
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  # set NA's to 1:
  colSums(is.na(DF_merged))
  
  DF_merged[is.na(DF_merged)] = 1
  
  pval = data.frame(
    DEXSeq_TECs = DF_merged$p_val_DEXSeq_TECs,
    DEXSeq_ECs = DF_merged$p_val_DEXSeq_ECs,
    DRIMSeq = DF_merged$p_val_DRIMSeq,
    DifferentialRegulation = DF_merged$Prob_DifferentialRegulation,
    satuRn = DF_merged$p_val_satuRn,
    SUPPA2 = DF_merged$p_val_SUPPA2,
    eisaR = DF_merged$p_val_eisaR
  )
  
  # AUC calculation:
  library(pROC)
  AUC[[type]] = apply(pval, 2, function(x){
    auc(DF_merged$truth, x)[1]
  })
  #}
  padj = data.frame(
    DEXSeq_TECs = DF_merged$p_adj_DEXSeq_TECs,
    DEXSeq_ECs = DF_merged$p_adj_DEXSeq_ECs,
    DRIMSeq = DF_merged$p_adj_DRIMSeq,
    DifferentialRegulation = DF_merged$Prob_DifferentialRegulation,
    satuRn = DF_merged$p_adj_satuRn,
    SUPPA2 = DF_merged$p_adj_SUPPA2,
    eisaR = DF_merged$p_adj_eisaR
  )
  
  truth = data.frame(status = DF_merged$truth)
  
  rownames(pval) <- DF_merged$Transcript_id
  rownames(padj) <- DF_merged$Transcript_id
  rownames(truth) <- DF_merged$Transcript_id
  
  DF_COBRA <- COBRAData(
    pval,
    padj,
    truth = truth)
  
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  # ROC, FDR plots:
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  perf <- calculate_performance(DF_COBRA, binary_truth = "status", 
                                aspects = c("roc", "fpc", "fdrtpr"),
                                thrs = c(0.01, 0.05, 0.1, 0.2))
  cobra_plot <- prepare_data_for_plot(perf, colorscheme = all_colours, incloverall = FALSE,
                                      facetted = TRUE,conditionalfill = FALSE,
                                      keepmethods = methods_all)
  gg_roc_with_pval <- plot_roc(cobra_plot,linewidth=2) + 
    scale_color_manual(values = all_colours,
                       name = "",
                       breaks=methods_all,
                       labels=methods_all) +
    my_theme
  
  # plot ROC curve
  ROC <- gg_roc_with_pval +
    geom_abline(slope = 1, intercept = 0 )
  
  # top N results:
  topN <- plot_fpc(cobra_plot, maxnfdc = 2000, linewidth = 2 ) +
    theme(legend.position="bottom") + my_theme+
    scale_color_manual(values = all_colours,
                       name = "",
                       breaks=methods_all,
                       labels=methods_all)
  
  ROC_type[[type]] = ROC
  topN_type[[type]] = topN
  
  save_path = file.path("5 - plots")
  ggsave(filename = paste0('ROC_', directories_results[type], '_2.pdf'),
         plot = ROC,
         device = "pdf",
         path = save_path,
         width = 20,
         height = 18,
         units = "in",
         dpi = 300,
         limitsize = TRUE)
  
  ggsave(filename = paste0('TopN_', directories_results[type], '_2.pdf'),
         plot = topN,
         device = "pdf",
         path = save_path,
         width = 20,
         height = 18,
         units = "in",
         dpi = 300,
         limitsize = TRUE)
  
  # for each type, make 1 plot (1 row and 2 columns) with ROC + topN - ggarrange
  # 1 unique legend at the bottom
  ROC_topN <- ggpubr::ggarrange(ROC, topN, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")
  ggsave(filename = paste0('ROC_topN_', directories_results[type], '_bulk1.pdf'),
         plot = ROC_topN,
         device = "pdf",
         path = save_path,
         width = 35,
         height = 20,
         units = "in",
         dpi = 300,
         limitsize = TRUE)
  
  # 
  ggsave(filename = paste0('ROC_topN_', directories_results[type], '_bulk2.pdf'),
         plot = ROC_topN,
         device = "pdf",
         path = save_path,
         width = 25,
         height = 15,
         units = "in",
         dpi = 300,
         limitsize = TRUE)
  
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  # stratify plots by abundance:
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  sel_rows = list()
  
  # split in low, medium and high counts:
  cuts = quantile(DF_merged$counts, c(1/3, 2/3))
  
  sel_rows[[1]] = DF_merged$counts <= cuts[1]
  sel_rows[[2]] = (DF_merged$counts > cuts[1]) & (DF_merged$counts <= cuts[2])
  sel_rows[[3]] = DF_merged$counts > cuts[2]
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  # iCOBRA:
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
  # set NA's to 1:
  colSums(is.na(DF_merged))
  
  DF_merged[is.na(DF_merged)] = 1
  
  plot_ROC_list[[type]] = list()
  plot_topN_list[[type]] = list()
  
  for(sel in 1:3){
    DF_sel = DF_merged[sel_rows[[sel]], ]
    
    pval = data.frame(
      DEXSeq_TECs = DF_sel$p_val_DEXSeq_TECs,
      DEXSeq_ECs = DF_sel$p_val_DEXSeq_ECs,
      DRIMSeq = DF_sel$p_val_DRIMSeq,
      DifferentialRegulation = DF_sel$Prob_DifferentialRegulation,
      satuRn = DF_sel$p_val_satuRn,
      SUPPA2 = DF_sel$p_val_SUPPA2,
      eisaR = DF_sel$p_val_eisaR
    )
    
    truth = data.frame(status = DF_sel$truth)
    
    rownames(pval) <- DF_sel$Transcript_id
    rownames(truth) <- DF_sel$Transcript_id
    
    DF_COBRA <- COBRAData(
      pval,
      truth = truth)
    
    #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
    # ROC, FDR plots:
    #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
    perf <- calculate_performance(DF_COBRA, binary_truth = "status", 
                                  aspects = c("roc", "fpc", "fdrtpr"),
                                  thrs = c(0.01, 0.05, 0.1, 0.2))
    cobra_plot <- prepare_data_for_plot(perf, colorscheme = all_colours, incloverall = FALSE,
                                        facetted = TRUE,conditionalfill = FALSE,
                                        keepmethods = methods_all)
    
    gg_roc_with_pval <- plot_roc(cobra_plot,linewidth=2) + 
      scale_color_manual(values = all_colours,
                         name = "",
                         breaks=methods_all,
                         labels=methods_all) +
      my_theme
    
    # plot ROC curve
    plot_ROC_list[[type]][[sel]] <- gg_roc_with_pval +
      geom_abline(slope = 1, intercept = 0 )
    
    topN <- plot_fpc(cobra_plot, maxnfdc = 700, linewidth = 2 ) +
      theme(legend.position="bottom") + my_theme +
      scale_color_manual(values = all_colours,
                         name = "",
                         breaks=methods_all,
                         labels=methods_all)
    
    plot_topN_list[[type]][[sel]] = topN
  }
}

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# AUC:
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
A = do.call(rbind, AUC)
overall = colMeans(A)

A = A[,order(overall)]
if(FALSE){
  library(xtable)
  xtable(A, digits =3)
  \begin{table}[ht]
  \centering
  \begin{tabular}{rrrrrrrr}
  \hline
  & DEXSeq\_ECs & DRIMSeq & DEXSeq\_TECs & SUPPA2 & satuRn & eisaR & DifferentialRegulation \\ 
  \hline
  1 & 0.411 & 0.425 & 0.864 & 0.908 & 0.905 & 0.938 & 0.951 \\ 
  2 & 0.411 & 0.432 & 0.858 & 0.903 & 0.897 & 0.940 & 0.951 \\ 
  3 & 0.415 & 0.448 & 0.857 & 0.874 & 0.884 & 0.940 & 0.950 \\ 
  \hline
  \end{tabular}
  \end{table}
  
  xtable(A, digits =2)
  \begin{table}[ht]
  \centering
  \begin{tabular}{rrrrrrrr}
  \hline
  & DEXSeq\_ECs & DRIMSeq & DEXSeq\_TECs & SUPPA2 & satuRn & eisaR & DifferentialRegulation \\ 
  \hline
  1 & 0.41 & 0.43 & 0.86 & 0.91 & 0.90 & 0.94 & 0.95 \\ 
  2 & 0.41 & 0.43 & 0.86 & 0.90 & 0.90 & 0.94 & 0.95 \\ 
  3 & 0.42 & 0.45 & 0.86 & 0.87 & 0.88 & 0.94 & 0.95 \\ 
  \hline
  \end{tabular}
  \end{table}
}
# make an overall plot of the 9 images in plot_ROC_list
# 1 unique legend at the bottom
# 3 rows = low, mid and high abundance (sel 1-3)
# 3 cols = type 2-4 (DR, DR + DGE, DR + DAS)

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# 9 ROC curves together:
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
legend <- ggpubr::get_legend(plot_ROC_list[[2]][[1]] +
                               guides(colour = guide_legend(nrow = 2, byrow = FALSE,
                                                            override.aes = list(#shape = shape_fill[-8],
                                                              fill = all_colours[-8]) ) ) +
                               theme(legend.key.height=unit(3,"line"),
                                     legend.text=element_text(size=35),
                                     legend.key.width=unit(6,"line")))

ggpubr::as_ggplot(legend)
my_theme = theme(legend.position = "none", 
                 text = element_text(size = 9),
                 plot.title = element_text(hjust = 0, face = "bold", size=35),
                 plot.subtitle = element_text(hjust = 0.5, face = "bold", size=30))

AA1 = ggpubr::ggarrange(plot_ROC_list[[2]][[1]] + xlab("") + labs(title = "Low abundance", subtitle = "DR") + my_theme,
                        plot_ROC_list[[3]][[1]] + xlab("") + ylab("") + labs(subtitle = "DR + DGE") + my_theme,
                        plot_ROC_list[[4]][[1]] + xlab("") + ylab("") + labs(subtitle = "DR + DAS") + my_theme,
                        ncol = 3, nrow = 1)

AA2 = ggpubr::ggarrange(plot_ROC_list[[2]][[2]] + xlab("") + labs(title = "Mid abundance") + my_theme,
                        plot_ROC_list[[3]][[2]] + xlab("") + ylab("")  + my_theme,
                        plot_ROC_list[[4]][[2]] + xlab("") + ylab("") + my_theme,
                        ncol = 3, nrow = 1)

AA3 = ggpubr::ggarrange(plot_ROC_list[[2]][[3]] + labs(title = "High abundance") + my_theme,
                        plot_ROC_list[[3]][[3]] + ylab("") + my_theme,
                        plot_ROC_list[[4]][[3]] + ylab("") + my_theme,
                        ncol = 3, nrow = 1)

AA = ggpubr::ggarrange( AA1, AA2, AA3, legend, 
                        # font.label = list(color = "black", size = 40, face = "bold"),
                        #label.x = 0, labrl.y = 1,
                        #labels = c("LIBD", "melanoma", "mouse cerebellum"),
                        heights=c(3,3,3,1),
                        widths = c(1,1,1,1),align = "v",
                        ncol = 1, nrow = 4 )

ggsave(filename = paste0('Sup_topN_bulk.pdf'),
       plot = AA,
       device = "pdf",
       path = save_path,
       width = 30,
       height = 33,
       units = "in",
       dpi = 300,
       limitsize = TRUE)

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# 9 topN curves together:
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
AA1 = ggpubr::ggarrange(plot_topN_list[[2]][[1]] + xlab("") + labs(title = "Low abundance", subtitle = "DR") + my_theme,
                        plot_topN_list[[3]][[1]] + xlab("") + ylab("") + labs(subtitle = "DR + DGE") + my_theme,
                        plot_topN_list[[4]][[1]] + xlab("") + ylab("") + labs(subtitle = "DR + DAS") + my_theme,
                        ncol = 3, nrow = 1)

AA2 = ggpubr::ggarrange(plot_topN_list[[2]][[2]] + xlab("") + labs(title = "Mid abundance") + my_theme,
                        plot_topN_list[[3]][[2]] + xlab("") + ylab("")  + my_theme,
                        plot_topN_list[[4]][[2]] + xlab("") + ylab("") + my_theme,
                        ncol = 3, nrow = 1)

AA3 = ggpubr::ggarrange(plot_topN_list[[2]][[3]] + labs(title = "High abundance") + my_theme,
                        plot_topN_list[[3]][[3]] + ylab("") + my_theme,
                        plot_topN_list[[4]][[3]] + ylab("") + my_theme,
                        ncol = 3, nrow = 1)

AA = ggpubr::ggarrange( AA1, AA2, AA3, legend, 
                        # font.label = list(color = "black", size = 40, face = "bold"),
                        #label.x = 0, labrl.y = 1,
                        #labels = c("LIBD", "melanoma", "mouse cerebellum"),
                        heights=c(3,3,3,1),
                        widths = c(1,1,1,1),align = "v",
                        ncol = 1, nrow = 4 )

ggsave(filename = paste0('Sup_ROC_bulk.pdf'),
       plot = AA,
       device = "pdf",
       path = save_path,
       width = 30,
       height = 33,
       units = "in",
       dpi = 300,
       limitsize = TRUE)



#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# 3 ROC together:
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
GG_ROC = ggpubr::ggarrange(ROC_type[[2]] + xlab("") + labs(subtitle = "DR") + my_theme,
                           ROC_type[[3]] + ylab("") + labs(subtitle = "DR + DGE") + my_theme,
                           ROC_type[[4]] + xlab("") + ylab("") + labs(subtitle = "DR + DAS") + my_theme,
                           legend = "bottom", 
                           common.legend = TRUE,
                           ncol = 3, nrow = 1)

ggsave(filename = paste0('3_ROC_bulk.pdf'),
       plot = GG_ROC,
       device = "pdf",
       path = save_path,
       width = 30,
       height = 12,
       units = "in",
       dpi = 300,
       limitsize = TRUE)

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# 3 topN together:
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
GG_topN = ggpubr::ggarrange(topN_type[[2]] + xlab("") + labs(subtitle = "DR") + my_theme,
                            topN_type[[3]] + ylab("") + labs(subtitle = "DR + DGE") + my_theme,
                            topN_type[[4]] + xlab("") + ylab("") + labs(subtitle = "DR + DAS") + my_theme,
                            legend = "bottom", 
                            common.legend = TRUE,
                            ncol = 3, nrow = 1)

ggsave(filename = paste0('3_topN_bulk.pdf'),
       plot = GG_topN,
       device = "pdf",
       path = save_path,
       width = 30,
       height = 12,
       units = "in",
       dpi = 300,
       limitsize = TRUE)

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# 3 ROC and 3 topN together:
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
GG_ROC_topN = ggpubr::ggarrange(ROC_type[[2]] + labs(subtitle = "DR") + my_theme,
                                ROC_type[[3]] + ylab("") + labs(subtitle = "DR + DGE") + my_theme,
                                ROC_type[[4]] + ylab("") + labs(subtitle = "DR + DAS") + my_theme,
                                topN_type[[2]] + my_theme,
                                topN_type[[3]] + ylab("") + my_theme,
                                topN_type[[4]] + ylab("") + my_theme,
                                legend = "bottom", 
                                common.legend = TRUE,
                                ncol = 3, nrow = 2)

ggsave(filename = paste0('3_ROC_3_topN_bulk.pdf'),
       plot = GG_ROC_topN,
       device = "pdf",
       path = save_path,
       width = 30,
       height = 22,
       units = "in",
       dpi = 300,
       limitsize = TRUE)

