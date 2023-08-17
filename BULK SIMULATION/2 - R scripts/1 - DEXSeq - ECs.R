export LC_ALL=C
# upgrade pip locally (user):
# pip install --upgrade pip --user

# update python pandas:
# python3 -m pip install --upgrade pandas --user
# git clone git@github.com:Oshlack/ec-dtu-pipe.git

ec="/home/Shared/simone/DifferentialRegulation/"
############################################################
# 'NULL':
############################################################
time python3 /home/Shared/simone/DifferentialRegulation/software/ec-dtu-pipe/create_salmon_ec_count_matrix.py \
/home/Shared/simone/DifferentialRegulation/"1 - data"/de0_ds0_dr0/2_quants/salmon/sample1/aux_info/eq_classes.txt.gz \
/home/Shared/simone/DifferentialRegulation/"1 - data"/de0_ds0_dr0/2_quants/salmon/sample2/aux_info/eq_classes.txt.gz \
/home/Shared/simone/DifferentialRegulation/"1 - data"/de0_ds0_dr0/2_quants/salmon/sample3/aux_info/eq_classes.txt.gz \
/home/Shared/simone/DifferentialRegulation/"1 - data"/de0_ds0_dr0/2_quants/salmon/sample4/aux_info/eq_classes.txt.gz \
/home/Shared/simone/DifferentialRegulation/"1 - data"/de0_ds0_dr0/2_quants/salmon/sample5/aux_info/eq_classes.txt.gz \
/home/Shared/simone/DifferentialRegulation/"1 - data"/de0_ds0_dr0/2_quants/salmon/sample6/aux_info/eq_classes.txt.gz \
sample1,sample2,sample3,sample4,sample5,sample6 \
/home/Shared/simone/DifferentialRegulation/"3 - results"/NULL/ec_matrix.txt

# real    2m29.880s
# user    2m20.468s
# sys     0m13.159s
############################################################
# 'DR':
############################################################
time python3 /home/Shared/simone/DifferentialRegulation/software/ec-dtu-pipe/create_salmon_ec_count_matrix.py \
/home/Shared/simone/DifferentialRegulation/"1 - data"/de0_ds0_dr2000/2_quants/salmon/sample1/aux_info/eq_classes.txt.gz \
/home/Shared/simone/DifferentialRegulation/"1 - data"/de0_ds0_dr2000/2_quants/salmon/sample2/aux_info/eq_classes.txt.gz \
/home/Shared/simone/DifferentialRegulation/"1 - data"/de0_ds0_dr2000/2_quants/salmon/sample3/aux_info/eq_classes.txt.gz \
/home/Shared/simone/DifferentialRegulation/"1 - data"/de0_ds0_dr2000/2_quants/salmon/sample4/aux_info/eq_classes.txt.gz \
/home/Shared/simone/DifferentialRegulation/"1 - data"/de0_ds0_dr2000/2_quants/salmon/sample5/aux_info/eq_classes.txt.gz \
/home/Shared/simone/DifferentialRegulation/"1 - data"/de0_ds0_dr2000/2_quants/salmon/sample6/aux_info/eq_classes.txt.gz \
sample1,sample2,sample3,sample4,sample5,sample6 \
/home/Shared/simone/DifferentialRegulation/"3 - results"/DR/ec_matrix.txt

# real    3m10.513s
# user    3m0.615s
# sys     0m14.199s
############################################################
# 'DR + DGE':
############################################################
time python3 /home/Shared/simone/DifferentialRegulation/software/ec-dtu-pipe/create_salmon_ec_count_matrix.py \
/home/Shared/simone/DifferentialRegulation/"1 - data"/de2000_ds0_dr2000/2_quants/salmon/sample1/aux_info/eq_classes.txt.gz \
/home/Shared/simone/DifferentialRegulation/"1 - data"/de2000_ds0_dr2000/2_quants/salmon/sample2/aux_info/eq_classes.txt.gz \
/home/Shared/simone/DifferentialRegulation/"1 - data"/de2000_ds0_dr2000/2_quants/salmon/sample3/aux_info/eq_classes.txt.gz \
/home/Shared/simone/DifferentialRegulation/"1 - data"/de2000_ds0_dr2000/2_quants/salmon/sample4/aux_info/eq_classes.txt.gz \
/home/Shared/simone/DifferentialRegulation/"1 - data"/de2000_ds0_dr2000/2_quants/salmon/sample5/aux_info/eq_classes.txt.gz \
/home/Shared/simone/DifferentialRegulation/"1 - data"/de2000_ds0_dr2000/2_quants/salmon/sample6/aux_info/eq_classes.txt.gz \
sample1,sample2,sample3,sample4,sample5,sample6 \
/home/Shared/simone/DifferentialRegulation/"3 - results"/'DR + DGE'/ec_matrix.txt
# real    2m19.612s
# user    2m11.024s
# sys     0m13.296s
############################################################
# 'DR + DAS':
############################################################
time python3 /home/Shared/simone/DifferentialRegulation/software/ec-dtu-pipe/create_salmon_ec_count_matrix.py \
/home/Shared/simone/DifferentialRegulation/"1 - data"/de0_ds2000_dr2000/2_quants/salmon/sample1/aux_info/eq_classes.txt.gz \
/home/Shared/simone/DifferentialRegulation/"1 - data"/de0_ds2000_dr2000/2_quants/salmon/sample2/aux_info/eq_classes.txt.gz \
/home/Shared/simone/DifferentialRegulation/"1 - data"/de0_ds2000_dr2000/2_quants/salmon/sample3/aux_info/eq_classes.txt.gz \
/home/Shared/simone/DifferentialRegulation/"1 - data"/de0_ds2000_dr2000/2_quants/salmon/sample4/aux_info/eq_classes.txt.gz \
/home/Shared/simone/DifferentialRegulation/"1 - data"/de0_ds2000_dr2000/2_quants/salmon/sample5/aux_info/eq_classes.txt.gz \
/home/Shared/simone/DifferentialRegulation/"1 - data"/de0_ds2000_dr2000/2_quants/salmon/sample6/aux_info/eq_classes.txt.gz \
sample1,sample2,sample3,sample4,sample5,sample6 \
/home/Shared/simone/DifferentialRegulation/"3 - results"/'DR + DAS'/ec_matrix.txt
# real    2m5.829s
# user    1m58.065s
# sys     0m11.426s

## in R:
alias R='R_LIBS=/home/stiber/R/x86_64-pc-linux-gnu-library/4.3 /home/Shared/simone/DifferentialRegulation/software/R-4.3.0/bin/R'

cd /home/Shared/simone/DifferentialRegulation

R

rm(list =ls())
setwd("/home/Shared/simone/DifferentialRegulation")
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

############################################################
# load data:
############################################################
library(tximport)
library(dplyr)
library(data.table)
library(DEXSeq)
library(DRIMSeq)

TIME = list()
for(type in 1:4){
  TIME[[type]] = system.time({
    set.seed(169612)
    
    data_dir = file.path("1 - data", directories[type],
                         "3_truth/simulation_details.txt")
    
    truth = fread(data_dir)
    ensg <- data.frame(gene_name = truth$gene_id, gene_id = truth$gene_id, transcript = truth$transcript_id,
                       stringsAsFactors = FALSE)
    ############################################################
    # DEXSeq - ECCs:
    ############################################################
    # We are now ready to load our data:
    ecm <- read.delim(file.path('3 - results', directories_results[type], 
                                'ec_matrix.txt'), sep='\t')
    ecm <- inner_join(ecm, ensg, by='transcript')
    
    # create transcript id without S or U version
    tr = ecm$transcript
    tr = strsplit(tr, "-")
    tr = sapply(tr, function(x){
      x[[1]]
    })
    ecm$gene_id = tr
    ecm$gene_name = tr
    
    # Next we remove equivalence classes that map to multiple different genes:
    multi_ecs <- data.table(ecm)[, length(unique(gene_id)), keyby=ec_names]
    multi_ecs <- multi_ecs[multi_ecs$V1 > 1]
    multi_ecs <- multi_ecs$ec_names
    df <- ecm[!ecm$ec_names %in% multi_ecs,]
    
    # We now define our sample groups:
    samples = paste0("sample", seq_len(6))
    group_names = rep(c("A", "B"), each = 3)
    
    sampleTable = data.frame(sample_id = samples,
                             condition = group_names )
    
    # Our EC data frame still potentially contains repeat entries 
    # as multiple transcripts mapping to the same EC will be divided into multiple rows.
    # Therefore, we will remove that transcript ID column, 
    # and consider all distinct equivalence classes per gene. 
    # As we will use DRIMSeq to prepare and filter our data,
    # we will also change create a 'feature_id' column which holds the equivalence class IDs:
    # Our EC data frame still potentially contains repeat entries as multiple transcripts mapping to the same EC will be divided into multiple rows. Therefore, we will remove that transcript ID column, and consider all distinct equivalence classes per gene. As we will use DRIMSeq to prepare and filter our data, we will also change create a 'feature_id' column which holds the equivalence class IDs:
    df$feature_id <- df$ec_names
    df <- distinct(df[,c(samples, 'gene_id', 'feature_id')])
    
    # Our data is now ready to load into a DRIMSeq object.
    d <- dmDSdata(counts = df, samples = sampleTable)
    
    sample.data <- DRIMSeq::samples(d)
    count.data <- round(as.matrix(counts(d)[,-c(1:2)]))
    
    # Our data is now ready to load into DEXSeq to perform differential transcript usage analysis:
    dxd <- DEXSeqDataSet(countData=count.data,
                         sampleData=sample.data,
                         design=~sample + exon + condition:exon,
                         featureID=counts(d)$feature_id,
                         groupID=counts(d)$gene_id)
    
    # We will use MulticoreParam to allow DEXSeq to take advantage of multiple for cores when performing analysis. This greatly speeds up analysis. First, we estimate the library sizes, then estimate dispersions, which we can visualise:
    BPPARAM = MulticoreParam(workers=1)
    dxd <- estimateSizeFactors(dxd)
    dxd <- estimateDispersions(dxd, BPPARAM=BPPARAM)
    
    dxd <- testForDEU(dxd, BPPARAM=BPPARAM)
    res <- DEXSeqResults(dxd)
    
    qval = perGeneQValue(res)
    res_gene = data.frame(gene = names(qval), qval)
  })
  ############################################################
  # save results:
  ############################################################
  data_dir = file.path("3 - results",
                       directories_results[type])
  
  name = paste("DEXSeq_ECCs")
  
  name = paste0(name, ".RData")
  
  save(res, res_gene, file = file.path(data_dir,name))
}

#OUTSIDE THE LOOP, save time:
name = paste0("TIME_", name)
save(TIME, file = file.path("3 - results", name))
