gtf="/home/Shared_sherborne/simone/Diff_Velo/discovery_bulk/reference/gencode.vM33.annotation.gz.expanded.gtf"
suppa="/home/Shared/simone/DifferentialRegulation/software/SUPPA"
output="/home/Shared/simone/Diff_Velo/discovery_bulk/04_results"

export LC_ALL=C
#python3 -m pip install --upgrade pandas --user
#python3 -m pip install --upgrade scikit-learn --user
#python3 -m pip install --upgrade statsmodels --user
#python3 -m pip install --upgrade numexpr --user
########## To generate the events from the GTF file one has to run the following command:
time python3 $suppa/suppa.py generateEvents -i $gtf -o $output/events -f ioi

########## need to compute an expression_file: TPM with 1 col per sample, 1 row per transcript
# in R:
alias R='R_LIBS=/home/stiber/R/x86_64-pc-linux-gnu-library/4.3 /home/Shared/simone/DifferentialRegulation/software/R-4.3.0/bin/R'

cd /home/Shared/simone/Diff_Velo/discovery_bulk/
R

rm(list =ls())
############################################################
# load data:
############################################################
library(data.table); library(tximport)
TIME = system.time({
  set.seed(169612)
  data_dir = "03_Salmon"
  
  sample_names = paste0("ERR45792", 67:72)
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
})
data_dir = "04_results"

write.table(TPM[,1:3], file = file.path(data_dir, "expression_file_Cond_A.csv"), row.names=TRUE, col.names=TRUE, sep="\t")
write.table(TPM[,4:6], file = file.path(data_dir, "expression_file_Cond_B.csv"), row.names=TRUE, col.names=TRUE, sep="\t")
# \t for tab separated.

# OPEN WITH TEXTEDIT AND then manually remove quotes.

#OUTSIDE THE LOOP, save time:
name = paste("SUPPA2.RData")
name = paste0("TIME_R1_", name)
save(TIME, file = file.path("04_results",name))

# remove quotes " " from files.

############################################################
# DR:
############################################################
########## At the moment the PSI per transcript isoform is calculated in the following way:
time python3 /home/Shared/simone/DifferentialRegulation/software/SUPPA/suppa.py \
psiPerIsoform -g $gtf \
-e $output/expression_file_Cond_A.csv \
-o $output/psi_A

time python3 /home/Shared/simone/DifferentialRegulation/software/SUPPA/suppa.py \
psiPerIsoform -g $gtf \
-e $output/expression_file_Cond_B.csv \
-o $output/psi_B

########## Differential splicing:
time python3 /home/Shared/simone/DifferentialRegulation/software/SUPPA/suppa.py \
diffSplice -m empirical -gc -i $output/events.ioi \
-p $output/psi_A_isoform.psi $output/psi_B_isoform.psi \
-e $output/expression_file_Cond_A.csv $output/expression_file_Cond_B.csv \
-o $output/DTU_OUTPUT

########## Sort results in R:
# DTU with transcripts:
# in R:
alias R='R_LIBS=/home/stiber/R/x86_64-pc-linux-gnu-library/4.3 /home/Shared/simone/DifferentialRegulation/software/R-4.3.0/bin/R'

cd /home/Shared/simone/Diff_Velo/discovery_bulk

R
rm(list =ls())
############################################################
# load data:
############################################################
library(data.table); library(tximport)
set.seed(169612)
data_dir = file.path("04_results", "DTU_OUTPUT.dpsi")

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
data_dir = "04_results"

name = paste("SUPPA2.RData")
save(res_SUPPA, file = file.path(data_dir,name))
