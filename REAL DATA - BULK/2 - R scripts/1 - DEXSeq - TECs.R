alias R='R_LIBS=/home/stiber/R/x86_64-pc-linux-gnu-library/4.3 /home/Shared/simone/DifferentialRegulation/software/R-4.3.0/bin/R'

cd /home/Shared/simone/Diff_Velo/discovery_bulk

R

rm(list =ls())
############################################################
# load data:
############################################################
library(tximport)
library(DEXSeq)
TIME = system.time({
  set.seed(169612)
  data_dir = "03_Salmon"
  
  sample_names = paste0("ERR45792", 67:72)
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
  # FILTER LOWLY ABUNDANT TRANSCRIPTS:
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
  # 97539
  sel = transcript_id %in% transcripts_kept
  counts = counts[sel,]
  transcript_id = transcript_id[sel]
  
  ############################################################
  # DEXSeq:
  ############################################################
  # Load truth table:
  colnames(counts) = sample_names
  group_names = rep(c("A", "B"), each = 3)
  
  design = data.frame( row.names = sample_names,
                       condition = group_names )
  
  dxd = DEXSeqDataSet(countData = round( counts ),
                      sampleData=design,
                      design= ~ sample + exon + condition:exon,
                      featureID = rownames(counts),
                      groupID = transcript_id, 
                      transcripts = rownames(counts))
  
  dxd <- estimateSizeFactors(dxd)
  dxd <- estimateDispersions(dxd) #, BPPARAM = BPPARAM)
  dxd <- testForDEU(dxd, reducedModel = ~sample + exon) #, BPPARAM = BPPARAM)
  
  res <- DEXSeqResults(dxd, independentFiltering = FALSE)
  
  qval = perGeneQValue(res)
  res_gene = data.frame(gene = names(qval), qval)
}) 
############################################################
# save results:
############################################################
data_dir = "04_results"

name = paste("DEXSeq_TECs.RData")

save(res, res_gene, file = file.path(data_dir,name))

#OUTSIDE THE LOOP, save time:
name = paste0("TIME_", name)
save(TIME, file = file.path(data_dir, name))
