alias R='R_LIBS=/home/stiber/R/x86_64-pc-linux-gnu-library/4.3 /home/Shared/simone/DifferentialRegulation/software/R-4.3.0/bin/R'

cd /home/Shared/simone/DifferentialRegulation

R

rm(list =ls())
setwd("/home/Shared/simone/DifferentialRegulation")
############################################################
# choose directories:
############################################################
# TRANSCRIPT ESTIMATED COUNTS:
directories = c("de0_ds0_dr2000_drbatch2000")

directories_results = c("batch")

############################################################
# load data:
############################################################
library(tximport)
library(DEXSeq)

set.seed(169612)

data_dir = file.path("1 - data", directories,
                     "2_quants/salmon" )

sample_names = paste0("sample", seq_len(6))
quant_files = file.path(data_dir, sample_names, "quant.sf")
file.exists(quant_files)

txi <- tximport(files = quant_files, type = "salmon", txOut = TRUE, dropInfReps = TRUE)
counts <- txi$counts

# get transcript id:
transcript_id = sapply(rownames(counts), function(x){
  strsplit(x, "-")[[1]][[1]]
})
names(transcript_id) = NULL

############################################################
# FILTER LOWLY ABUNDANT TRANSCRIPTSa:
############################################################
# keep transcripts with >= 10 counts per condition (across all samples):
counts_split = split(counts, transcript_id)
sel = sapply(counts_split, function(x){
  # sel samples in group A:
  sel_a = 1:6
  (sum(x[sel_a]) >= 10) & (sum(x[-sel_a]) >= 10)
})
transcripts_kept = names(sel[sel])
length(transcripts_kept)
# 46647
sel = transcript_id %in% transcripts_kept
counts = counts[sel,]
transcript_id = transcript_id[sel]

############################################################
# DEXSeq:
############################################################

# Load truth table:
colnames(counts) = sample_names
group_names = rep(c("A", "B"), each = 3)

batch = c("A", "A", "B", "A", "B", "B")

design = data.frame( row.names = sample_names,
                     condition = group_names,
                     batch = batch)

fullModel = ~ sample + exon + batch:exon + condition:exon
reducedModel = ~ sample + exon + batch:exon

dxd = DEXSeqDataSet(countData = round( counts ),
                    sampleData=design,
                    design=fullModel,
                    featureID = rownames(counts),
                    groupID = transcript_id, 
                    transcripts = rownames(counts))

dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd)
dxd <- testForDEU(dxd, 
                  reducedModel=reducedModel, 
                  fullModel=fullModel)

res <- DEXSeqResults(dxd, independentFiltering = FALSE)

qval = perGeneQValue(res)
res_gene = data.frame(gene = names(qval), qval)

############################################################
# save results:
############################################################
data_dir = file.path("3 - results",
                     directories_results)

name = paste("DEXSeq_TECs")

name = paste0(name, ".RData")

save(res, res_gene, file = file.path(data_dir,name))
