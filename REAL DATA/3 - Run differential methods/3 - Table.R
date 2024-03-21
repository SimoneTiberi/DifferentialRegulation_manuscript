setwd("~/Desktop/Differential Velocity/FINAL SCRIPTS/REAL DATA")
# clean enviroment
rm(list = ls())

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# load associated genes
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
BUZZ_WORDS = list.files("02_results/TERMS")
BUZZ_WORDS = substring(BUZZ_WORDS, 1,  nchar(BUZZ_WORDS) - 5)

#brain only detected            # detected in human brain only
#cerebral cortex elevated       # elevated in human cortex, compared to other regions of the brain
#excitatory neurons             # excitatory neurons - role in development
#PGP1                           # brain organoids we have

library(rjson)
library(doRNG)
GENES <- foreach(i =  1:length(BUZZ_WORDS), .combine = "c") %do% {
  a = unlist(fromJSON(file = paste0("02_results/TERMS/", BUZZ_WORDS[i],".json")))
  sel = names(a) == "Ensembl"
  list(a[sel])
}
str(GENES)
names(GENES) <- BUZZ_WORDS

sapply(GENES, length)

GENES$overall = unique(unlist(GENES))

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# DifferentialRegulation
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
load("02_results/DifferentialRegulation.RData")
results_EC = results_EC[[1]][,c(1:2, 3, 6) ]

colnames(results_EC)[3:4] = c("DifferentialRegulation_Wald", 
                              "DifferentialRegulation")

results_EC[,4] = sapply(results_EC[,4], function(x){
  2*min(x, 1-x)
})

DF_merged = results_EC

head(DF_merged);
rm(results_EC)

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# eisaR:
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# load DEXSeq SUA 02_results:
load("02_results/eisaR.RData")
DF_eisaR = DF_eisaR[,1:3]
colnames(DF_eisaR)[3] = "eisaR"

DF_merged = merge(DF_merged, DF_eisaR, by = c("Gene_id", "Cluster_id"), all.x = TRUE)
head(DF_merged)

rm(DF_eisaR)

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# DRIMSeq:
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# load DRIMSeq SUA 02_results:
load("02_results/DRIMSeq.RData")
DF_DRIMSeq = DF_DRIMSeq[,1:3] 

colnames(DF_DRIMSeq)[3] = c("DRIMSeq")

DF_merged = merge(DF_merged, DF_DRIMSeq, by = c("Gene_id", "Cluster_id"), all.x = TRUE)
head(DF_merged)

rm(DF_DRIMSeq)

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# DEXSeq:
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# load DEXSeq SUA 02_results:
load("02_results/DEXSeq_USA.RData")
colnames(DF_DEXSeq) = c("DEXSeq", "Gene_id", "Cluster_id")

DF_merged = merge(DF_merged, DF_DEXSeq, by = c("Gene_id", "Cluster_id"), all.x = TRUE)
head(DF_merged)

rm(DF_DEXSeq)

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# satuRn:
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# load satuRn 02_results:
load("02_results/SatuRn.RData")

n = length(DF_SatuRn)

satuRn = list()
for(cl in 1:n){
  res = DF_SatuRn[[cl]]
  cluster = res$cluster_id[1]
  
  tr_names = res$transcript_id
  tr_names = strsplit(tr_names, "-")
  tr_names = sapply(tr_names, function(x){ x[[1]]})
  
  res$empirical_pval[is.na(res$empirical_pval)] = 1
  res_by_tr = split(res$empirical_pval, tr_names)
  min_p_val = sapply(res_by_tr, min, na.rm = TRUE)
  
  qval = p.adjust(min_p_val, method = "BH")
  
  satuRn[[cl]] = data.frame(Gene_id = names(min_p_val),
                            Cluster_id = cluster,
                            satuRn = min_p_val)#,
  #satuRn_FDR = qval) 
}
satuRn = do.call(rbind, satuRn)

DF_merged = merge(DF_merged, satuRn, by = c("Gene_id", "Cluster_id"), all.x = TRUE)
head(DF_merged)
tail(DF_merged)

rm(DF_SatuRn); rm(satuRn)

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# satuRn - SC:
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# load satuRn 02_results:
load("02_results/SatuRnSCsampleDesign.RData")

n = length(DF_SatuRn)

satuRn = list()
for(cl in 1:n){
  res = DF_SatuRn[[cl]]
  cluster = res$cluster_id[1]
  
  tr_names = res$transcript_id
  tr_names = strsplit(tr_names, "-")
  tr_names = sapply(tr_names, function(x){ x[[1]]})
  
  res$empirical_pval[is.na(res$empirical_pval)] = 1
  res_by_tr = split(res$empirical_pval, tr_names)
  min_p_val = sapply(res_by_tr, min, na.rm = TRUE)
  
  qval = p.adjust(min_p_val, method = "BH")
  
  satuRn[[cl]] = data.frame(Gene_id = names(min_p_val),
                            Cluster_id = cluster,
                            satuRn_SC = min_p_val)
}
satuRn = do.call(rbind, satuRn)

DF_merged = merge(DF_merged, satuRn, by = c("Gene_id", "Cluster_id"), all.x = TRUE)
head(DF_merged)
tail(DF_merged)

rm(DF_SatuRn); rm(satuRn)

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# BRIE2:
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# load DEXSeq SUA 02_results:
load("02_results/BRIE2_USA.RData")
head(RES)
RES = RES[,-3]
head(RES)
colnames(RES) = c("Gene_id", "BRIE2", "Cluster_id")

RES$Cluster_id = ifelse(RES$Cluster_id == "Immature_CPNs", "Immature CPNs", RES$Cluster_id)
RES$Cluster_id = ifelse(RES$Cluster_id == "Immature_PNs", "Immature PNs", RES$Cluster_id)

table(RES$Cluster_id)
table(DF_merged$Cluster_id)

DF_merged = merge(DF_merged, RES, by = c("Gene_id", "Cluster_id"), all.x = TRUE)
head(DF_merged)

table(DF_merged$Cluster_id)

rm(RES)

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# top results:
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# rm what comes after "." in gene id:
DF_merged$Gene_id = substr(DF_merged$Gene_id, 1, 15)

DF_merged[is.na(DF_merged)] = 1

methods = c("BRIE2",
  "DEXSeq",
  "DifferentialRegulation",
  "DifferentialRegulation_Wald",
  "DRIMSeq",
  "eisaR",
  "satuRn",
  "satuRn_SC")

n_top = 200

clusters = unique(DF_merged$Cluster_id); clusters

top_xx = lapply(clusters, function(cluster){
  DF_one_cluster = DF_merged[DF_merged$Cluster_id == cluster, ]
  
  res = lapply(methods, function(method){
    sel = which(colnames(DF_one_cluster) == method )
    if(method == "DifferentialRegulation"){
      # resolve ties in DifferentialRegulation output, based on the result of DifferentialRegulation Wald test:
      sel_wald = which(colnames(DF_one_cluster) == "DifferentialRegulation_Wald" )
      ordering = order(DF_one_cluster[,sel], DF_one_cluster[,sel_wald])
    }else{
      ordering = order(DF_one_cluster[,sel])
    }
    genes = DF_one_cluster$Gene_id[ ordering ]
    genes[1:n_top]
  })
  res = as.data.frame(do.call(cbind, res))
  colnames(res) = methods
  res$cluster = cluster
  res
})
top_xx = do.call(rbind,top_xx)

found = round(apply(top_xx[,1:length(methods)],2, function(xx){
  sapply(GENES, function(gg){
    sum(xx %in% gg)
  })
}))
found

genes_in_DF = sapply(GENES, function(x){
  sum(sapply(x, function(y){
    sum(y == DF_merged$Gene_id)
  }))
})/length(clusters)
genes_in_DF

# how many we'd expect at random
expected_n = genes_in_DF/nrow(DF_merged) * n_top * length(clusters)
expected_n

found = found[,order(found[5,], decreasing = TRUE)]

found

library(xtable)
xtable(found, digits = 0)

\begin{table}[ht]
\centering
\begin{tabular}{rrrrrrrrr}
\hline
& satuRn\_SC & DifferentialRegulation & DifferentialRegulation\_Wald & BRIE2 & eisaR & satuRn & DEXSeq & DRIMSeq \\ 
\hline
brain only detected & 8 & 3 & 3 & 0 & 0 & 0 & 0 & 0 \\ 
cerebral cortex elevated & 8 & 7 & 7 & 0 & 0 & 1 & 0 & 0 \\ 
excitatory neurons & 3 & 2 & 1 & 1 & 0 & 0 & 0 & 0 \\ 
PGP1 & 0 & 2 & 2 & 0 & 1 & 0 & 0 & 0 \\ 
overall & 17 & 14 & 13 & 1 & 1 & 1 & 0 & 0 \\ 
\hline
\end{tabular}
\end{table}


n_top = 200

clusters = unique(DF_merged$Cluster_id); clusters

top_xx = lapply(clusters, function(cluster){
  DF_one_cluster = DF_merged[DF_merged$Cluster_id == cluster, ]
  
  res = lapply(methods, function(method){
    sel = which(colnames(DF_one_cluster) == method )
    if(method == "DifferentialRegulation"){
      # resolve ties in DifferentialRegulation output, based on the result of DifferentialRegulation Wald test:
      sel_wald = which(colnames(DF_one_cluster) == "DifferentialRegulation_Wald" )
      ordering = order(DF_one_cluster[,sel], DF_one_cluster[,sel_wald])
    }else{
      ordering = order(DF_one_cluster[,sel])
    }
    genes = DF_one_cluster$Gene_id[ ordering ]
    genes[1:n_top]
  })
  res = as.data.frame(do.call(cbind, res))
  colnames(res) = methods
  res$cluster = cluster
  res
})
top_xx = do.call(rbind,top_xx)

found = round(apply(top_xx[,1:length(methods)],2, function(xx){
  sapply(GENES, function(gg){
    sum(xx %in% gg)
  })
}))
found

genes_in_DF = sapply(GENES, function(x){
  sum(sapply(x, function(y){
    sum(y == DF_merged$Gene_id)
  }))
})/length(clusters)
genes_in_DF

# how many we'd expect at random
expected_n = genes_in_DF/nrow(DF_merged) * n_top * length(clusters)
expected_n

found = found[,order(found[5,], decreasing = TRUE)]

found

library(xtable)
xtable(found, digits = 0)

\begin{table}[ht]
\centering
\begin{tabular}{rrrrrrrrr}
\hline
& satuRn\_SC & DifferentialRegulation & DifferentialRegulation\_Wald & BRIE2 & eisaR & satuRn & DEXSeq & DRIMSeq \\ 
\hline
brain only detected & 8 & 3 & 3 & 0 & 0 & 0 & 0 & 0 \\ 
cerebral cortex elevated & 8 & 7 & 7 & 0 & 0 & 1 & 0 & 0 \\ 
excitatory neurons & 3 & 2 & 1 & 1 & 0 & 0 & 0 & 0 \\ 
PGP1 & 0 & 2 & 2 & 0 & 1 & 0 & 0 & 0 \\ 
overall & 17 & 14 & 13 & 1 & 1 & 1 & 0 & 0 \\ 
\hline
\end{tabular}
\end{table}
