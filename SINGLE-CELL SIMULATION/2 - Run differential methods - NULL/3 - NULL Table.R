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

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# load ground truth:
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
sim = readRDS("03_data/mouse_simulation_data-NULL-ALEVIN.rds")

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
save_name = paste0("07_results/DifferentialRegulation_NULL.RData")

load(save_name)

results_EC = lapply(results_EC, function(X) X[[1]])
results_EC = do.call(rbind, results_EC)

colnames(results_EC)[1:2] = c("gene_id", "cluster_id")
colnames(results_EC)[3:6] = paste0("DifferentialRegulation_", colnames(results_EC)[3:6])

results_EC$ORIGINAL_PROB = results_EC[,6]

results_EC[,6] = sapply(results_EC[,6], function(x){
  2*min(x, 1-x)
})

DF_merged = merge(DF_merged, results_EC, by = c("gene_id", "cluster_id"))
head(DF_merged); rm(results_EC)

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# eisaR:
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
save_name = paste0("07_results/results_NULL_eisaR.RData")
load(save_name)

colnames(DF_eisaR) = c("gene_id", "cluster_id",
                       "eisaR", "eisaR_FDR")

DF_merged = merge(DF_merged, DF_eisaR, by = c("gene_id", "cluster_id"), all.x = TRUE)
head(DF_merged)

rm(DF_eisaR);

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# DEXSeq - USA:
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
save_name = paste0("07_results/results_NULL_DEXSeq_USA.RData")
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
save_name = paste0("07_results/results_NULL_DRIMSeq.RData")
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
load("07_results/BRIE2_USA_NULL.RData")

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
load("07_results/results_NULL_SatuRn.RData")

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
# satuRn:
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
load("07_results/results_NULL_SatuRnSCsampleDesign.RData")

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

FPs = 100 * data.frame( colMeans(pval <= 0.1),
                        colMeans(pval <= 0.05),
                        colMeans(pval <= 0.01))
round(FPs, 2)

df = data.frame(p = DF_merged$ORIGINAL_PROB)
library(ggplot2)
# Basic histogram
ggplot(df, aes(x=p)) + geom_histogram()
# Change colors
p = ggplot(df, aes(x=p)) + 
  geom_histogram(color="black", fill="white")
p

save_path = file.path("08_plots")
ggsave(filename = paste0('null_hist_sc.pdf'),
       plot = p,
       device = "pdf",
       path = save_path,
       width = 6,
       height = 6,
       units = "in",
       dpi = 300,
       limitsize = TRUE)

library(xtable)
xtable(FPs, digits = 1)

\begin{table}[ht]
\centering
\begin{tabular}{rrrr}
\hline
& colMeans.pval....0.1. & colMeans.pval....0.05. & colMeans.pval....0.01. \\ 
\hline
BRIE2 & 4.0 & 2.5 & 0.9 \\ 
DEXSeq\_TECs & 0.1 & 0.1 & 0.0 \\ 
DifferentialRegulation & 3.5 & 1.7 & 0.4 \\ 
DifferentialRegulation\_Wald & 1.5 & 0.8 & 0.3 \\ 
DRIMSeq & 1.3 & 0.6 & 0.2 \\ 
eisaR & 4.9 & 2.1 & 0.3 \\ 
satuRn & 34.4 & 23.2 & 9.2 \\ 
satuRn\_SC & 21.1 & 12.5 & 4.6 \\ 
\hline
\end{tabular}
\end{table}
