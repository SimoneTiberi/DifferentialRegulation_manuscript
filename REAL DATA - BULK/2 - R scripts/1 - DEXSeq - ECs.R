export LC_ALL=C
# upgrade pip locally (user):
# pip install --upgrade pip --user

# update python pandas:
# python3 -m pip install --upgrade pandas --user
# git clone git@github.com:Oshlack/ec-dtu-pipe.git

############################################################
# 'discovery':
############################################################
time python3 /home/Shared/simone/DifferentialRegulation/software/ec-dtu-pipe/create_salmon_ec_count_matrix.py \
/home/Shared_sherborne/simone/Diff_Velo/discovery_bulk/03_Salmon/ERR4579267/aux_info/eq_classes.txt.gz \
/home/Shared_sherborne/simone/Diff_Velo/discovery_bulk/03_Salmon/ERR4579268/aux_info/eq_classes.txt.gz \
/home/Shared_sherborne/simone/Diff_Velo/discovery_bulk/03_Salmon/ERR4579269/aux_info/eq_classes.txt.gz \
/home/Shared_sherborne/simone/Diff_Velo/discovery_bulk/03_Salmon/ERR4579270/aux_info/eq_classes.txt.gz \
/home/Shared_sherborne/simone/Diff_Velo/discovery_bulk/03_Salmon/ERR4579271/aux_info/eq_classes.txt.gz \
/home/Shared_sherborne/simone/Diff_Velo/discovery_bulk/03_Salmon/ERR4579272/aux_info/eq_classes.txt.gz \
ERR4579267,ERR4579268,ERR4579269,ERR4579270,ERR4579271,ERR4579272 \
/home/Shared_sherborne/simone/Diff_Velo/discovery_bulk/04_results/ec_matrix.txt

############################################################
#R:
############################################################
## in R:
alias R='R_LIBS=/home/stiber/R/x86_64-pc-linux-gnu-library/4.3 /home/Shared/simone/DifferentialRegulation/software/R-4.3.0/bin/R'

cd /home/Shared/simone/Diff_Velo/discovery_bulk

R

rm(list =ls())
############################################################
# load data:
############################################################
library(tximport)
library(dplyr)
library(data.table)
library(DEXSeq)
library(DRIMSeq)

TIME = system.time({
  set.seed(169612)
  truth = fread("reference/gencode.vM33.annotation.gz.expanded.tx2gene.tsv", header = FALSE)
  colnames(truth) = c("transcript_id", "gene_id")
  
  # TODO -> create this one!
  ensg <- data.frame(gene_name = truth$gene_id, gene_id = truth$gene_id, transcript = truth$transcript_id,
                     stringsAsFactors = FALSE)
  ############################################################
  # DEXSeq - ECCs:
  ############################################################
  # We are now ready to load our data:
  ecm <- read.delim("04_results/ec_matrix.txt", sep='\t')
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
  samples = paste0("ERR45792", 67:72)
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
data_dir = file.path("04_results")

name = paste("DEXSeq_ECCs.RData")

save(res, res_gene, file = file.path(data_dir,name))

#OUTSIDE THE LOOP, save time:
name = paste0("TIME_", name)
save(TIME, file = file.path(data_dir, name))
