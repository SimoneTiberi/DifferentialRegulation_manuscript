setwd("~/Desktop/Differential Velocity/FINAL SCRIPTS/REAL DATA - BULK/04_results")

# clean enviroment
rm(list = ls())
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# load associated genes
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
BUZZ_WORDS = list.files("TERMS")
# remove ".json" at the end
BUZZ_WORDS = substring(BUZZ_WORDS, 1,  nchar(BUZZ_WORDS) - 5)

# mouse pancreas development
# ZNF800

BUZZ_WORDS = BUZZ_WORDS

library(rjson)
library(doRNG)
GENES <- foreach(i =  1:length(BUZZ_WORDS), .combine = "c") %do% {
  a = unlist(fromJSON(file = paste0("TERMS/", BUZZ_WORDS[i],".json")))
  sel = names(a) == "Gene"
  list(a[sel])
}
str(GENES)
names(GENES) <- BUZZ_WORDS

GENES = lapply(GENES, toupper)
GENES = lapply(GENES, unique)

sapply(GENES, length)

GENES$overall = unique(unlist(GENES))

head(GENES$overall)

############################################################
# load results:
############################################################
# DifferentialRegulation:
name = paste("DifferentialRegulation.RData")
load(name)
DR = res_EC$Differential_results
DR = DR[,c(1,2,4)]

names(DR)[2:3]  = c("DifferentialRegulation_Wald", "DifferentialRegulation")

DR[,3] = sapply(DR[,3], function(x){
  2*min(x, 1-x)
})

DF_merged = DR

rm(DR); rm(res_EC); rm(name)

# eisaR:
name = paste("eisaR.RData")
load(name)
DF_eisaR = data.frame(Transcript_id = DF_eisaR$Transcript_id,
                      eisaR = DF_eisaR$p_eisaR)
DF_merged = merge(DF_eisaR, DF_merged, by = "Transcript_id", all.y = TRUE)
head(DF_merged)
rm(DF_eisaR); rm(name)

#library(DEXSeq)
# DEXSeq - TECs:
name = paste("DEXSeq_TECs.RData")
load(name)
DEXSeq_TECs = data.frame(Transcript_id = res_gene$gene,
                         DEXSeq_TECs = res_gene$qval)
DF_merged = merge(DEXSeq_TECs, DF_merged, by = "Transcript_id", all.y = TRUE)
head(DF_merged)
rm(DEXSeq_TECs); rm(res); rm(name)

# DEXSeq - ECs:
name = paste("DEXSeq_ECCs.RData")
load(name)
DEXSeq_ECs = data.frame(Transcript_id = res_gene$gene,
                        DEXSeq_ECs = res_gene$qval)
DF_merged = merge(DEXSeq_ECs, DF_merged, by = "Transcript_id", all.y = TRUE)
head(DF_merged)
rm(DEXSeq_ECs); rm(res_gene); rm(name)

# DRIMSeq:
name = paste("DRIMSeq.RData")
load(name)
DRIMSeq = data.frame(Transcript_id = results_gene$gene_id,
                     DRIMSeq = results_gene$pvalue)
DF_merged = merge(DRIMSeq, DF_merged, by = "Transcript_id", all.y = TRUE)
head(DF_merged)
rm(DRIMSeq); rm(results_gene); rm(name)

# satuRn:
name = paste("satuRn.RData")
load(name)

tr_names = rownames(res)
tr_names = strsplit(tr_names, "-")
tr_names = sapply(tr_names, function(x){ x[[1]]})

res$empirical_pval[is.na(res$empirical_pval)] = 1
res_by_tr = split(res$empirical_pval, tr_names)
min_p_val = sapply(res_by_tr, min, na.rm = TRUE)

satuRn = data.frame(Transcript_id = names(min_p_val),
                    satuRn = min_p_val)
DF_merged = merge(satuRn, DF_merged, by = "Transcript_id", all.y = TRUE)
head(DF_merged)
rm(satuRn); rm(res); rm(name);
rm(min_p_val); rm(tr_names)

# SUPPA2:
name = paste("SUPPA2.RData")
load(name)

SUPPA2 = data.frame(Transcript_id = names(res_SUPPA),
                    SUPPA2 = res_SUPPA)

DF_merged = merge(SUPPA2, DF_merged, by = "Transcript_id", all.y = TRUE)
head(DF_merged)
rm(SUPPA2); rm(res_SUPPA); rm(name)

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# move from Transcript to Gene ids:
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
library(data.table)
tr_to_gene = fread("gencode.vM33.annotation.gz.expanded.tx2gene.tsv", header = FALSE)
colnames(tr_to_gene) = c("transcript_id", "gene_id")

matches = match(DF_merged$Transcript_id, tr_to_gene$transcript_id)
head(DF_merged$Transcript_id); head(tr_to_gene$transcript_id[matches])

DF_merged$Gene_id = tr_to_gene$gene_id[matches]

# rm what comes after "." in Transcript_id:
DF_merged$Gene_id = substr(DF_merged$Gene_id, 1, 18)

head(DF_merged$Gene_id)

library(gprofiler2)
a = gconvert( unique(DF_merged$Gene_id),
              organism="mmusculus",
              target="ENTREZGENE",
              filter_na = FALSE)
# a
length(unique(DF_merged$Gene_id)); dim(a)

DF = data.frame(Gene_id = a$input, Gene_name = a$name)
DF = unique(DF)

dim(DF_merged)
DF_merged = merge(DF_merged, DF, all.x = TRUE)
dim(DF_merged)

# remove genes with no match:
dim(DF_merged)
DF_merged = DF_merged[ !is.na(DF_merged$Gene_name), ]
dim(DF_merged)

# ~ 5% of results are removed

# all UPPER-CASE
DF_merged$Gene_name = toupper(DF_merged$Gene_name)

length(unique(DF_merged$Gene_name))
# 16,347

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# top results:
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
DF_merged[is.na(DF_merged)] = 1

methods = c("DEXSeq_ECs",
            "DEXSeq_TECs",
            "DifferentialRegulation",
            "DifferentialRegulation_Wald",
            "DRIMSeq",
            "eisaR",
            "satuRn",
            "SUPPA2")

n_top = 1000

res = sapply(methods, function(method){
  sel = which(colnames(DF_merged) == method )
  if(method == "DifferentialRegulation"){
    sel_wald = which(colnames(DF_merged) == "DifferentialRegulation_Wald" )
    ordering = order(DF_merged[,sel], DF_merged[,sel_wald])
  }else{
    ordering = order(DF_merged[,sel])
  }
  genes = unique(DF_merged$Gene_name[ ordering ])
  genes[1:n_top]
})
head(res)

found = round(apply(res[,1:length(methods)],2, function(xx){
  sapply(GENES, function(gg){
    sum(xx %in% gg)
  })
}))
found
sum(found)

genes_in_DF = sapply(GENES, function(x){
  y = 0
  if(length(x) > 0){
    y = sum(sapply(x, function(y){
      sum(y == DF_merged$Gene_name)
    }))
  }
  y
})
genes_in_DF

# how many we'd expect at random
expected_n = genes_in_DF/nrow(DF_merged) * n_top
expected_n

found = found[,order(found[nrow(found),], decreasing = TRUE)]

# found = cbind(found, round(expected_n,2))

found
round(found/expected_n, 1)

# Fig. 7 - Zfp800 knockout mice have defects in pancreatic endocrine and acinar tissue development
# and disturbances in pancreatic exocrine cell morphology

                              DifferentialRegulation DifferentialRegulation_Wald DEXSeq_TECs eisaR SUPPA2 DRIMSeq DEXSeq_ECs satuRn
mouse pancreas development                     12                          12          12    11     11      10          9      6
ZNF800                                          1                           1           0     0      0       0          0      0
overall                                        13                          13          12    11     11      10          9      6

unique_genes = unique(DF_merged$Gene_name)
unique_genes_in_DF = sapply(GENES, function(x){
  y = 0
  if(length(x) > 0){
    y = sum(sapply(x, function(y){
      sum(y == unique_genes)
    }))
  }
  y
})
# fraction of GENES that actually appear among the final results,
# and hence are potentially identifiable.
unique_genes_in_DF/sapply(GENES, length)




library(xtable)
xtable(found, digits = 0)
\begin{table}[ht]
\centering
\begin{tabular}{rrrrrrrrr}
\hline
& DifferentialRegulation & DifferentialRegulation\_Wald & DEXSeq\_TECs & eisaR & SUPPA2 & DRIMSeq & DEXSeq\_ECs & satuRn \\ 
\hline
mouse pancreas development & 12 & 12 & 12 & 11 & 11 & 10 & 9 & 6 \\ 
ZNF800 & 1 & 1 & 0 & 0 & 0 & 0 & 0 & 0 \\ 
overall & 13 & 13 & 12 & 11 & 11 & 10 & 9 & 6 \\ 
\hline
\end{tabular}
\end{table}

