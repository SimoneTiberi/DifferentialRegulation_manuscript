gtf="gencode.v41.annotation.protein_coding.expanded.gtf"

suppa="/home/Shared/simone/DifferentialRegulation/software/SUPPA"

export LC_ALL=C

########## need to compute an expression_file: TPM with 1 col per sample, 1 row per transcript
# in R:
alias R='R_LIBS=/home/stiber/R/x86_64-pc-linux-gnu-library/4.3 /home/Shared/simone/DifferentialRegulation/software/R-4.3.0/bin/R'

R

rm(list =ls())

############################################################
# choose directories:
############################################################
# TRANSCRIPT ESTIMATED COUNTS:
# TRANSCRIPT ESTIMATED COUNTS:
directories = c("de0_ds0_dr2000_drbatch2000")

directories_results = c("batch")

setwd("/home/Shared/simone/DifferentialRegulation")
############################################################
# load data:
############################################################
library(data.table); library(tximport)
set.seed(169612)
data_dir = file.path("1 - data", directories,
                     "2_quants/salmon" )

sample_names = paste0("sample", seq_len(6))
quant_files = file.path(data_dir, sample_names, "quant.sf")
file.exists(quant_files)

txi <- tximport(files = quant_files, type = "salmon", txOut = TRUE, dropInfReps = TRUE)
counts <- txi$counts
TPM <- txi$abundance
colnames(TPM) = paste0("sample", 1:6)
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
# length(transcript_id): 329476
# length(transcripts_kept): 50740
sel = transcript_id %in% transcripts_kept
TPM = TPM[sel,]

data_dir = file.path("3 - results",
                     directories_results)

write.table(TPM[,1:3], file = file.path(data_dir, "expression_file_Cond_A.csv"), row.names=TRUE, col.names=TRUE, sep="\t")
write.table(TPM[,4:6], file = file.path(data_dir, "expression_file_Cond_B.csv"), row.names=TRUE, col.names=TRUE, sep="\t")
# \t for tab separated.

# OPEN WITH TEXTEDIT AND then manually remove quotes.

############################################################
# 'NULL':
############################################################
########## At the moment the PSI per transcript isoform is calculated in the following way:
time python3 /home/Shared/simone/DifferentialRegulation/software/SUPPA/suppa.py \
psiPerIsoform -g /home/Shared/simone/DifferentialRegulation/'1 - data'/gencode.v41.annotation.protein_coding.expanded.gtf \
-e /home/Shared/simone/DifferentialRegulation/'3 - results'/batch/expression_file_Cond_A.csv \
-o /home/Shared/simone/DifferentialRegulation/'3 - results'/batch/psi_A

time python3 /home/Shared/simone/DifferentialRegulation/software/SUPPA/suppa.py \
psiPerIsoform -g /home/Shared/simone/DifferentialRegulation/'1 - data'/gencode.v41.annotation.protein_coding.expanded.gtf \
-e /home/Shared/simone/DifferentialRegulation/'3 - results'/batch/expression_file_Cond_B.csv \
-o /home/Shared/simone/DifferentialRegulation/'3 - results'/batch/psi_B

########## Differential splicing:
time python3 /home/Shared/simone/DifferentialRegulation/software/SUPPA/suppa.py \
diffSplice -m empirical -gc -i /home/Shared/simone/DifferentialRegulation/'3 - results'/events.ioi \
-p /home/Shared/simone/DifferentialRegulation/'3 - results'/batch/psi_A_isoform.psi /home/Shared/simone/DifferentialRegulation/'3 - results'/batch/psi_B_isoform.psi \
-e /home/Shared/simone/DifferentialRegulation/'3 - results'/batch/expression_file_Cond_A.csv /home/Shared/simone/DifferentialRegulation/'3 - results'/batch/expression_file_Cond_B.csv \
-o /home/Shared/simone/DifferentialRegulation/'3 - results'/batch/DTU_OUTPUT

########## Sort results in R:
# DTU with transcripts:
# in R:
alias R='R_LIBS=/home/stiber/R/x86_64-pc-linux-gnu-library/4.3 /home/Shared/simone/DifferentialRegulation/software/R-4.3.0/bin/R'

R
rm(list =ls())

############################################################
# choose directories:
############################################################
# TRANSCRIPT ESTIMATED COUNTS:
directories = c("de0_ds0_dr2000_drbatch2000")

directories_results = c("batch")

setwd("/home/Shared/simone/DifferentialRegulation")
############################################################
# load data:
############################################################
library(data.table); library(tximport)
set.seed(169612)
data_dir = file.path("3 - results", directories_results,
                     "DTU_OUTPUT.dpsi")
res_SUPPA = read.table(data_dir, header = TRUE)
head(res_SUPPA); dim(res_SUPPA)
tail(res_SUPPA); dim(res_SUPPA)

# get the gene names:              
gene_names = rownames(res_SUPPA)
# extract transcript name:
tr_names = strsplit(gene_names, ";")
tr_names = sapply(tr_names, function(x){ x[[2]]})

tr_names = strsplit(tr_names, "-")
tr_names = sapply(tr_names, function(x){ x[[1]]})

res_by_tr = split(res_SUPPA[,2], tr_names)
min_p_val = lapply(res_by_tr, min, na.rm = T)
min_p_val = unlist(min_p_val)
res_SUPPA = min_p_val
############################################################
# save results:
############################################################
data_dir = file.path("3 - results",
                     directories_results)

name = paste("SUPPA2.RData")

save(res_SUPPA, file = file.path(data_dir,name))
