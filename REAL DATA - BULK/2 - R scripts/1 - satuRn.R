alias R='R_LIBS=/home/stiber/R/x86_64-pc-linux-gnu-library/4.3 /home/Shared/simone/DifferentialRegulation/software/R-4.3.0/bin/R'
cd /home/Shared/simone/Diff_Velo/discovery_bulk

R

rm(list =ls())
############################################################
# load data:
############################################################
library(satuRn)
library(tximport)
library(SummarizedExperiment)

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
  
  sel = transcript_id %in% transcripts_kept
  counts = counts[sel,]
  transcript_id = transcript_id[sel]
  
  ############################################################
  # satuRn:
  ############################################################
  
  # Load truth table:
  colnames(counts) = sample_names
  group_names = rep(c("A", "B"), each = 3)
  
  design <- data.frame(row.names = sample_names,
                       sample_id = sample_names,
                       group = group_names)
  
  # the columns with transcript identifiers is names isoform_id, 
  # while the column containing gene identifiers should be named gene_id
  txInfo = data.frame(isoform_id = rownames(counts),
                      gene_id = transcript_id)
  
  # Create a SummarizedExperiment object
  sumExp <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    colData = design,
    rowData = txInfo
  )
  
  # for sake of completeness: specify design formula from colData
  metadata(sumExp)$formula <- ~ 0 + as.factor(colData(sumExp)$group)
  sumExp
  
  # We fit the model
  # ?satuRn::fitDTU
  sumExp <- satuRn::fitDTU(
    object = sumExp,
    formula = ~ 0 + group,
    parallel = FALSE,
    BPPARAM = MulticoreParam(1),
    verbose = TRUE
  )
  
  # construct design matrix
  group <- as.factor(design$group)
  design_full <- model.matrix(~ 0 + group)
  colnames(design_full) <- levels(group)
  
  # initialize contrast matrix
  L <- limma::makeContrasts(
    Contrast1 = A - B,
    levels = design_full
  )
  
  # We test the genes
  # ?satuRn::testDTU
  
  sumExp <- satuRn::testDTU(
    object = sumExp,
    contrasts = L,
    diagplot1 = TRUE,
    diagplot2 = TRUE,
    sort = FALSE
  )
  res <- rowData(sumExp)[["fitDTUResult_Contrast1"]]
})  
############################################################
# save results:
############################################################
data_dir = "04_results"

name = paste("satuRn.RData")

save(res, file = file.path(data_dir,name))

#OUTSIDE THE LOOP, save time:
name = paste0("TIME_", name)
save(TIME, file = file.path(data_dir, name))
