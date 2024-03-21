alias R='R_LIBS=/home/stiber/R/x86_64-pc-linux-gnu-library/4.3 /home/Shared/simone/DifferentialRegulation/software/R-4.3.0/bin/R'
cd /home/Shared/simone/DifferentialRegulation

R

rm(list =ls())

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
library(DRIMSeq)

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
# DRIMSeq:
############################################################

# Load truth table:
colnames(counts) = sample_names
group_names = rep(c("A", "B"), each = 3)

batch = c("A", "A", "B", "A", "B", "B")

design <- data.frame(sample_id = sample_names,
                     group = group_names,
                     batch = batch)

# save transcripts as "gene_id"
counts_df = data.frame(counts, 
                       gene_id = transcript_id,
                       feature_id = rownames(counts))

# Create a dmDSdata object
d <- dmDSdata(counts = counts_df, samples = design)
# with 46647 genes and 6 samples

design_full <- model.matrix(~ group + batch, data = samples(d))
design_full

# infer the precision parameters:
d <- dmPrecision(d, genewise_precision = TRUE, 
                 design = design_full)

# We fit the model
# ?dmFit
d <- dmFit(d, design = design_full, 
           verbose = 1)

# We test the genes
# ?dmTest
d <- dmTest(d, coef = "groupB", verbose = 1)

results_gene = results(d, level = "gene")

############################################################
# save results:
############################################################
data_dir = file.path("3 - results",
                     directories_results)

name = paste("DRIMSeq.RData")

save(results_gene, file = file.path(data_dir,name))
