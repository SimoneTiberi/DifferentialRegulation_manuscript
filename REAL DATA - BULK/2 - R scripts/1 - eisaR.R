alias R='R_LIBS=/home/stiber/R/x86_64-pc-linux-gnu-library/4.3 /home/Shared/simone/DifferentialRegulation/software/R-4.3.0/bin/R'

#BiocManager::install("Rbowtie")
#BiocManager::install("ShortRead")
#BiocManager::install("GenomicFiles")
# install_github("cran/deldir")

cd /home/Shared/simone/Diff_Velo/discovery_bulk

R

rm(list =ls())
############################################################
# load data:
############################################################
# BiocManager::install("tximport")

library(eisaR)
library(tximport)
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
  # eisaR:
  ############################################################
  # Load truth table:
  colnames(counts) = sample_names
  group_names = factor(rep(c("A", "B"), each = 3))
  
  sel_S = rownames(counts) %in% transcript_id
  # run eisaR
  R_ex = counts[sel_S,] # Spliced counts
  R_in = counts[!sel_S,] # Unsplced counts
  
  rownames(R_in) = rownames(R_ex)
  
  if(!all(dim(R_ex) == dim(R_in))){
    print('DIMENTIONS DO NOT MATCH!')
    break
  }
  
  set.seed(169612)
  tab <- runEISA(cntEx  = R_ex, cntIn = R_in, 
                 cond = group_names, 
                 geneSelection = "none")$tab.ExIn
  
  # initialize results data frame
  DF_eisaR = data.frame(Transcript_id = rownames(tab),
                        p_eisaR = tab$PValue,
                        p_eisaR_adj = tab$FDR)
})
############################################################
# save results:
############################################################
data_dir = "04_results"

name = paste("eisaR.RData")

save(DF_eisaR, file = file.path(data_dir,name))

#OUTSIDE THE LOOP, save time:
name = paste0("TIME_", name)
save(TIME, file = file.path(data_dir, name))
