#alias R='R_LIBS=/home/stiber/R/4.3.1 /usr/local/R/R-4.3.1/bin/R'
#cd /home/Shared/simone/DifferentialRegulation
#R

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
                "de0_ds2000_dr2000",
                "de2000_ds0_dr2000_logFC_6",
                "de2000_ds0_dr2000_logFC_9")

directories_results = c("NULL",
                        "DR",
                        "DR + DGE",
                        "DR + DAS",
                        "log2FC_6",
                        "log2FC_9")
type = 1

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
  
  DR$ORIGINAL_PROB = DR[,4]
  
  DR[,4] = sapply(DR[,4], function(x){
    2*min(x, 1-x)
  })
  DF_merged = merge(DR, DF_merged, by = "Transcript_id")
  head(DF_merged)
  rm(DR); rm(name)
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
  DifferentialRegulation_Wald = DF_merged$p_val_DifferentialRegulation,
  satuRn = DF_merged$p_val_satuRn,
  SUPPA2 = DF_merged$p_val_SUPPA2,
  eisaR = DF_merged$p_val_eisaR
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

save_path = file.path("5 - plots")
ggsave(filename = paste0('null_hist_bulk.pdf'),
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
DEXSeq\_TECs & 0.1 & 0.1 & 0.1 \\ 
DEXSeq\_ECs & 0.0 & 0.0 & 0.0 \\ 
DRIMSeq & 0.9 & 0.4 & 0.1 \\ 
DifferentialRegulation & 0.1 & 0.0 & 0.0 \\ 
DifferentialRegulation\_Wald & 0.1 & 0.0 & 0.0 \\ 
satuRn & 5.3 & 2.2 & 0.4 \\ 
SUPPA2 & 1.7 & 0.5 & 0.0 \\ 
eisaR & 1.2 & 0.4 & 0.0 \\ 
\hline
\end{tabular}
\end{table}
