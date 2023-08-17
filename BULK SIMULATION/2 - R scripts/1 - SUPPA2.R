gtf="gencode.v41.annotation.protein_coding.expanded.gtf"

suppa="/home/Shared/simone/DifferentialRegulation/software/SUPPA"

export LC_ALL=C
# python3.11 -m pip install scikit-learn
# python3.11 -m pip install statsmodels
########## To generate the events from the GTF file one has to run the following command:
time python3.11 /home/Shared/simone/DifferentialRegulation/software/SUPPA/suppa.py \
generateEvents -i /home/Shared/simone/DifferentialRegulation/'1 - data'/gencode.v41.annotation.protein_coding.expanded.gtf \
-o /home/Shared/simone/DifferentialRegulation/'3 - results'/events -f ioi
########## need to compute an expression_file: TPM with 1 col per sample, 1 row per transcript
# in R:
alias R='R_LIBS=/home/stiber/R/x86_64-pc-linux-gnu-library/4.3 /home/Shared/simone/DifferentialRegulation/software/R-4.3.0/bin/R'

R

rm(list =ls())

############################################################
# choose directories:
############################################################
# TRANSCRIPT ESTIMATED COUNTS:
directories = c("de0_ds0_dr0",
                "de0_ds0_dr2000",
                "de2000_ds0_dr2000",
                "de0_ds2000_dr2000")

directories_results = c("NULL",
                        "DR",
                        "DR + DGE",
                        "DR + DAS")
setwd("/home/Shared/simone/DifferentialRegulation")
############################################################
# load data:
############################################################
library(data.table); library(tximport)
TIME = list()
for(type in 1:4){
  TIME[[type]] = system.time({
    set.seed(169612)
    data_dir = file.path("1 - data", directories[type],
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
  })
  data_dir = file.path("3 - results",
                       directories_results[type])
  
  write.table(TPM[,1:3], file = file.path(data_dir, "expression_file_Cond_A.csv"), row.names=TRUE, col.names=TRUE, sep="\t")
  write.table(TPM[,4:6], file = file.path(data_dir, "expression_file_Cond_B.csv"), row.names=TRUE, col.names=TRUE, sep="\t")
  # \t for tab separated.
  
  # OPEN WITH TEXTEDIT AND then manually remove quotes.
}

#OUTSIDE THE LOOP, save time:
name = paste("SUPPA2")

name = paste0(name, ".RData")
name = paste0("TIME_R1_", name)
save(TIME, file = file.path("3 - results",name))

scp -r  /Users/peiyingcai/Desktop/"3 - results"/DR/expression_file_Cond_A.csv  peicai@imlssherborne:/home/Shared/simone/DifferentialRegulation/"3 - results"/DR/expression_file_Cond_A.csv
scp -r  /Users/peiyingcai/Desktop/"3 - results"/DR/expression_file_Cond_B.csv  peicai@imlssherborne:/home/Shared/simone/DifferentialRegulation/"3 - results"/DR/expression_file_Cond_B.csv
scp -r  /Users/peiyingcai/Desktop/"3 - results"/NULL/expression_file_Cond_A.csv  peicai@imlssherborne:/home/Shared/simone/DifferentialRegulation/"3 - results"/NULL/expression_file_Cond_A.csv
scp -r  /Users/peiyingcai/Desktop/"3 - results"/NULL/expression_file_Cond_B.csv  peicai@imlssherborne:/home/Shared/simone/DifferentialRegulation/"3 - results"/NULL/expression_file_Cond_B.csv
scp -r  /Users/peiyingcai/Desktop/"3 - results"/"DR + DAS"/expression_file_Cond_B.csv  peicai@imlssherborne:/home/Shared/simone/DifferentialRegulation/"3 - results"/"DR + DAS"/expression_file_Cond_B.csv
scp -r  /Users/peiyingcai/Desktop/"3 - results"/"DR + DAS"/expression_file_Cond_A.csv  peicai@imlssherborne:/home/Shared/simone/DifferentialRegulation/"3 - results"/"DR + DAS"/expression_file_Cond_A.csv
scp -r  /Users/peiyingcai/Desktop/"3 - results"/"DR + DGE"/expression_file_Cond_A.csv  peicai@imlssherborne:/home/Shared/simone/DifferentialRegulation/"3 - results"/"DR + DGE"/expression_file_Cond_A.csv
scp -r  /Users/peiyingcai/Desktop/"3 - results"/"DR + DGE"/expression_file_Cond_B.csv  peicai@imlssherborne:/home/Shared/simone/DifferentialRegulation/"3 - results"/"DR + DGE"/expression_file_Cond_B.csv


########## At the moment the PSI per transcript isoform is calculated in the following way:
python3.11 /home/Shared/simone/DifferentialRegulation/software/SUPPA-2.3/suppa.py \
psiPerIsoform -g $gtf -e expression_file_Cond_A.csv -o psi_A

python3.11 /home/Shared/simone/DifferentialRegulation/software/SUPPA-2.3/suppa.py \
psiPerIsoform -g $gtf -e expression_file_Cond_B.csv -o psi_B

########## Differential splicing:
python3.11 $suppa/suppa.py \
diffSplice -m empirical -gc -i events.ioi \
-p psi_A_isoform.psi psi_B_isoform.psi \
-e expression_file_Cond_A.csv expression_file_Cond_B.csv \
-o DTU_OUTPUT

############################################################
# 'NULL':
############################################################
########## At the moment the PSI per transcript isoform is calculated in the following way:
time python3.11 /home/Shared/simone/DifferentialRegulation/software/SUPPA/suppa.py \
psiPerIsoform -g /home/Shared/simone/DifferentialRegulation/'1 - data'/gencode.v41.annotation.protein_coding.expanded.gtf \
-e /home/Shared/simone/DifferentialRegulation/'3 - results'/NULL/expression_file_Cond_A.csv \
-o /home/Shared/simone/DifferentialRegulation/'3 - results'/NULL/psi_A

# real    1m13.010s
# user    0m57.590s
# sys     0m11.485s

time python3.11 /home/Shared/simone/DifferentialRegulation/software/SUPPA/suppa.py \
psiPerIsoform -g /home/Shared/simone/DifferentialRegulation/'1 - data'/gencode.v41.annotation.protein_coding.expanded.gtf \
-e /home/Shared/simone/DifferentialRegulation/'3 - results'/NULL/expression_file_Cond_B.csv \
-o /home/Shared/simone/DifferentialRegulation/'3 - results'/NULL/psi_B

# real    1m52.006s
# user    0m56.477s
# sys     0m11.782s

########## Differential splicing:
time python3.11 /home/Shared/simone/DifferentialRegulation/software/SUPPA/suppa.py \
diffSplice -m empirical -gc -i /home/Shared/simone/DifferentialRegulation/'3 - results'/events.ioi \
-p /home/Shared/simone/DifferentialRegulation/'3 - results'/NULL/psi_A_isoform.psi /home/Shared/simone/DifferentialRegulation/'3 - results'/NULL/psi_B_isoform.psi \
-e /home/Shared/simone/DifferentialRegulation/'3 - results'/NULL/expression_file_Cond_A.csv /home/Shared/simone/DifferentialRegulation/'3 - results'/NULL/expression_file_Cond_B.csv \
-o /home/Shared/simone/DifferentialRegulation/'3 - results'/NULL/DTU_OUTPUT

# real    28m55.760s
# user    28m48.702s
# sys     0m15.718s

############################################################
# DR:
############################################################
########## At the moment the PSI per transcript isoform is calculated in the following way:
time python3.11 /home/Shared/simone/DifferentialRegulation/software/SUPPA/suppa.py \
psiPerIsoform -g /home/Shared/simone/DifferentialRegulation/'1 - data'/gencode.v41.annotation.protein_coding.expanded.gtf \
-e /home/Shared/simone/DifferentialRegulation/'3 - results'/DR/expression_file_Cond_A.csv \
-o /home/Shared/simone/DifferentialRegulation/'3 - results'/DR/psi_A

# real    1m6.424s
# user    1m4.804s
# sys     0m11.747s

time python3.11 /home/Shared/simone/DifferentialRegulation/software/SUPPA/suppa.py \
psiPerIsoform -g /home/Shared/simone/DifferentialRegulation/'1 - data'/gencode.v41.annotation.protein_coding.expanded.gtf \
-e /home/Shared/simone/DifferentialRegulation/'3 - results'/DR/expression_file_Cond_B.csv \
-o /home/Shared/simone/DifferentialRegulation/'3 - results'/DR/psi_B

# real    0m57.938s
# user    0m56.617s
# sys     0m11.560s

########## Differential splicing:
time python3.11 /home/Shared/simone/DifferentialRegulation/software/SUPPA/suppa.py \
diffSplice -m empirical -gc -i /home/Shared/simone/DifferentialRegulation/'3 - results'/events.ioi \
-p /home/Shared/simone/DifferentialRegulation/'3 - results'/DR/psi_A_isoform.psi /home/Shared/simone/DifferentialRegulation/'3 - results'/DR/psi_B_isoform.psi \
-e /home/Shared/simone/DifferentialRegulation/'3 - results'/DR/expression_file_Cond_A.csv /home/Shared/simone/DifferentialRegulation/'3 - results'/DR/expression_file_Cond_B.csv \
-o /home/Shared/simone/DifferentialRegulation/'3 - results'/DR/DTU_OUTPUT

# real    24m48.912s
# user    24m40.931s
# sys     0m16.843s

############################################################
# 'DR + DAS':
############################################################
########## At the moment the PSI per transcript isoform is calculated in the following way:
time python3.11 /home/Shared/simone/DifferentialRegulation/software/SUPPA/suppa.py \
psiPerIsoform -g /home/Shared/simone/DifferentialRegulation/'1 - data'/gencode.v41.annotation.protein_coding.expanded.gtf \
-e /home/Shared/simone/DifferentialRegulation/'3 - results'/'DR + DAS'/expression_file_Cond_A.csv \
-o /home/Shared/simone/DifferentialRegulation/'3 - results'/'DR + DAS'/psi_A

# real    0m58.567s
# user    0m57.040s
# sys     0m11.317s

time python3.11 /home/Shared/simone/DifferentialRegulation/software/SUPPA/suppa.py \
psiPerIsoform -g /home/Shared/simone/DifferentialRegulation/'1 - data'/gencode.v41.annotation.protein_coding.expanded.gtf \
-e /home/Shared/simone/DifferentialRegulation/'3 - results'/'DR + DAS'/expression_file_Cond_B.csv \
-o /home/Shared/simone/DifferentialRegulation/'3 - results'/'DR + DAS'/psi_B

# real    0m59.295s
# user    0m57.667s
# sys     0m11.187s

########## Differential splicing:
time python3.11 /home/Shared/simone/DifferentialRegulation/software/SUPPA/suppa.py \
diffSplice -m empirical -gc -i /home/Shared/simone/DifferentialRegulation/'3 - results'/events.ioi \
-p /home/Shared/simone/DifferentialRegulation/'3 - results'/'DR + DAS'/psi_A_isoform.psi /home/Shared/simone/DifferentialRegulation/'3 - results'/'DR + DAS'/psi_B_isoform.psi \
-e /home/Shared/simone/DifferentialRegulation/'3 - results'/'DR + DAS'/expression_file_Cond_A.csv /home/Shared/simone/DifferentialRegulation/'3 - results'/'DR + DAS'/expression_file_Cond_B.csv \
-o /home/Shared/simone/DifferentialRegulation/'3 - results'/'DR + DAS'/DTU_OUTPUT

# real    23m36.729s
# user    23m30.922s
# sys     0m13.082s

############################################################
# 'DR + DGE':
############################################################
########## At the moment the PSI per transcript isoform is calculated in the following way:
time python3.11 /home/Shared/simone/DifferentialRegulation/software/SUPPA/suppa.py \
psiPerIsoform -g /home/Shared/simone/DifferentialRegulation/'1 - data'/gencode.v41.annotation.protein_coding.expanded.gtf \
-e /home/Shared/simone/DifferentialRegulation/'3 - results'/'DR + DGE'/expression_file_Cond_A.csv \
-o /home/Shared/simone/DifferentialRegulation/'3 - results'/'DR + DGE'/psi_A

# real    0m58.810s
# user    0m56.798s
# sys     0m11.557s

time python3.11 /home/Shared/simone/DifferentialRegulation/software/SUPPA/suppa.py \
psiPerIsoform -g /home/Shared/simone/DifferentialRegulation/'1 - data'/gencode.v41.annotation.protein_coding.expanded.gtf \
-e /home/Shared/simone/DifferentialRegulation/'3 - results'/'DR + DGE'/expression_file_Cond_B.csv \
-o /home/Shared/simone/DifferentialRegulation/'3 - results'/'DR + DGE'/psi_B

# real    1m10.020s
# user    0m56.923s
# sys     0m10.756s

########## Differential splicing:
time python3.11 /home/Shared/simone/DifferentialRegulation/software/SUPPA/suppa.py \
diffSplice -m empirical -gc -i /home/Shared/simone/DifferentialRegulation/'3 - results'/events.ioi \
-p /home/Shared/simone/DifferentialRegulation/'3 - results'/'DR + DGE'/psi_A_isoform.psi /home/Shared/simone/DifferentialRegulation/'3 - results'/'DR + DGE'/psi_B_isoform.psi \
-e /home/Shared/simone/DifferentialRegulation/'3 - results'/'DR + DGE'/expression_file_Cond_A.csv /home/Shared/simone/DifferentialRegulation/'3 - results'/'DR + DGE'/expression_file_Cond_B.csv \
-o /home/Shared/simone/DifferentialRegulation/'3 - results'/'DR + DGE'/DTU_OUTPUT

# real    18m13.924s
# user    18m6.487s
# sys     0m18.098s

# real    26m7.168s
# user    25m57.346s
# sys     0m15.216s
# SUPPA2 includes an option to correct for multiple testing using the Benjamini-Hochberg method across all events
# from the same gene, as they cannot be considered to be entirely independent of each other, for which the false discovery rate (FDR) cut-off can be given as input.

# -gc, --gene-correction
# Boolean. If True, SUPPA correct the p-values by gene. (Default: False).
########## Sort results in R:
# DTU with transcripts:
# in R:
alias R=‘R_LIBS=/home/stiber/R/x86_64-pc-linux-gnu-library/4.3 /home/Shared/simone/DifferentialRegulation/software/R-4.3.0/bin/R’

R
rm(list =ls())

############################################################
# choose directories:
############################################################
# TRANSCRIPT ESTIMATED COUNTS:
directories = c("de0_ds0_dr0",
                "de0_ds0_dr2000",
                "de2000_ds0_dr2000",
                "de0_ds2000_dr2000")

directories_results = c("NULL",
                        "DR",
                        "DR + DGE",
                        "DR + DAS")
setwd("/home/Shared/simone/DifferentialRegulation")
############################################################
# load data:
############################################################
library(data.table); library(tximport)
for(type in 1:4){
  set.seed(169612)
  data_dir = file.path("3 - results", directories_results[type],
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
                       directories_results[type])
  
  name = paste("SUPPA2.RData")
  
  save(res_SUPPA, file = file.path(data_dir,name))
}


TIME = list()

## NULL
TIME[[1]] = c(39.662, 66.241, 73.558, 1286.571)
## DR
TIME[[2]] = c(39.662, 54.453, 55.122, 1030.645)
## DR + DGE
TIME[[3]] = c(39.662, 52.251, 52.291, 1093.924)
## DR + DAS
TIME[[4]] = c(39.662, 58.838, 53.862, 953.522)
save(TIME, file = file.path("3 - results","TIME_Py2_SUPPA2.RData"))