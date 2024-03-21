# alias R='R_LIBS=/home/stiber/R/4.3.1 /usr/local/R/R-4.3.1/bin/R'
# cd /home/Shared/simone/DifferentialRegulation

# R

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
directories = c("de0_ds0_dr2000",
                "de0_ds0_dr2000_drbatch2000")

directories_results = c("DR",
                        "batch")

############################################################
# specify type of data:
############################################################
# save the 3 plots below (ROC, FDR and top_N)
plot_ROC_list <- list()
plot_topN_list <- list()

ROC_type = topN_type = list()

AUC = list()

for(type in 1:2){
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
}


my_theme = theme(legend.position = "none", 
                 text = element_text(size = 9),
                 plot.title = element_text(hjust = 0, face = "bold", size=35),
                 plot.subtitle = element_text(hjust = 0.5, face = "bold", size=30))

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# 2 ROC and 2 topN together:
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
GG_ROC_topN = ggpubr::ggarrange(ROC_type[[1]] + labs(subtitle = "DR") + my_theme,
                                ROC_type[[2]] + ylab("") + labs(subtitle = "DR + batch") + my_theme,
                                topN_type[[1]] + my_theme,
                                topN_type[[2]] + ylab("") + my_theme,
                                legend = "bottom", 
                                common.legend = TRUE,
                                ncol = 2, nrow = 2)

save_path = file.path("5 - plots")
ggsave(filename = paste0('2_ROC_2_topN_batch.pdf'),
       plot = GG_ROC_topN,
       device = "pdf",
       path = save_path,
       width = 30,
       height = 22,
       units = "in",
       dpi = 300,
       limitsize = TRUE)

