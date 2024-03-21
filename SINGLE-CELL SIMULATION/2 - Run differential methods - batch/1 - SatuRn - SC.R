alias R='R_LIBS=/home/stiber/R/x86_64-pc-linux-gnu-library/4.3 /home/Shared/simone/DifferentialRegulation/software/R-4.3.0/bin/R'

cd /home/Shared_sherborne/simone/Diff_Velo/NEW_simulation/06_alevin_fry_USA

R
rm(list = ls())

library(DifferentialRegulation)
library(DEXSeq)
library(SingleCellExperiment)
library(BiocParallel)

# specify 4 samples ids:
sample_ids = paste0("normal", seq_len(4))

min_counts_per_group = 10

clusters = c("Adipocytes", "Epithelial_cells", "Hepatocytes")

data_dir = "simulation_batch"

TIMES = list()  

DF_SatuRn = list()

for(cl in seq_along(clusters)){
  
  TIMES[[cl]] = system.time({
    cluster = clusters[cl]
    
    # set directories of each sample input data (obtained via alevin-fry):
    base_dir = file.path(data_dir, paste0(sample_ids, "_", cluster) )
    file.exists(base_dir)
    
    # set paths to objects:
    path_to_counts = file.path(base_dir,"/alevin/quants_mat.mtx")
    path_to_cell_id = file.path(base_dir,"/alevin/quants_mat_rows.txt")
    path_to_gene_id = file.path(base_dir,"/alevin/quants_mat_cols.txt")
    
    file.exists(path_to_counts)
    file.exists(path_to_cell_id)
    file.exists(path_to_gene_id)
    
    # load USA counts:
    sce = load_USA(path_to_counts,
                   path_to_cell_id,
                   path_to_gene_id,
                   sample_ids)
    
    sce$group = factor(ifelse(sce$sample_id %in% c("normal1", "normal3"),
                              "A", "B"))
    
    ############################################################
    # sort data:
    ############################################################
    n_samples = nlevels(sce$sample_id); n_samples
    
    SUA = list()
    for(i in 1:n_samples){
      sel_cells = sce$sample_id == levels(sce$sample_id)[i]
      sce_one_sample = sce[, sel_cells]
      S = rowSums(assays(sce_one_sample)$spliced)
      U = rowSums(assays(sce_one_sample)$unspliced)
      A = rowSums(assays(sce_one_sample)$ambiguous)
      
      SUA[[i]] = cbind(S, U, A)
    }
    x = do.call(cbind, SUA)
    rm(sce_one_sample)
    
    # start PI_counts from relative counts estimates
    # start PI_SU from 50:50 or relative abundance from SU counts estimates
    genes = rownames(sce)
    n_genes = length(genes)
    
    group = factor(c("A", "B", "A", "B"))
    sel_A = c(1:3, 7:9) # select columns of group A
    sel = (rowSums(x[,sel_A]) >= min_counts_per_group) & (rowSums(x[,-sel_A]) >= min_counts_per_group)
    sel_genes = genes[sel]
    n_genes = length(sel_genes)
    sum(sel)
    
    # filter transcripts based on sel
    sce = sce[sel,]
    
    SUA = rbind(as.matrix(assays(sce)$spliced),
                as.matrix(assays(sce)$unspliced),
                as.matrix(assays(sce)$ambiguous))
    
    # prepare design matrix
    design <- data.frame(condition = factor(group))
    spliced_unspliced_ambiguous <- c(rep("S", n_genes), rep("U", n_genes),  rep("A", n_genes))
    
    rownames(SUA) = paste( rep(sel_genes, 3),
                           spliced_unspliced_ambiguous,
                           sep  = "-")
    ############################################################
    # satuRn:
    ############################################################
    # Load truth table:
    # prepare design matrix
    batch = factor(c("1", "1", "2", "2"))
    
    design <- data.frame(row.names = colnames(SUA),
                         group = sce$group,
                         batch = factor(ifelse(sce$sample_id %in% c("normal1", "normal2"), "1", "2")))
    
    # the columns with transcript identifiers is names isoform_id, 
    # while the column containing gene identifiers should be named gene_id
    txInfo = data.frame(isoform_id = rownames(SUA),
                        gene_id = rep(sel_genes, 3))
    
    # Create a SummarizedExperiment object
    sumExp <- SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = SUA),
      colData = design,
      rowData = txInfo
    )
    
    # for sake of completeness: specify design formula from colData
    metadata(sumExp)$formula <- ~ 0 + 
      as.factor(colData(sumExp)$group) +
      as.factor(colData(sumExp)$batch)
    sumExp
    
    # We fit the model
    # ?satuRn::fitDTU
    sumExp <- satuRn::fitDTU(
      object = sumExp,
      formula = ~ 0 + group + batch,
      parallel = FALSE,
      BPPARAM = MulticoreParam(1),
      verbose = TRUE
    )
    
    # construct design matrix
    group <- as.factor(design$group)
    batch <- as.factor(design$batch)
    design_full <- model.matrix(~ 0 + group + batch)
    colnames(design_full)[1:2] <- levels(group)
    
    # initialize contrast matrix
    L <- limma::makeContrasts(
      Contrast1 = A - B,
      levels = design_full
    )
    
    # We test the genes
    # ?satuRn::testDTU
    set.seed(169612)
    sumExp <- satuRn::testDTU(
      object = sumExp,
      contrasts = L,
      diagplot1 = TRUE,
      diagplot2 = TRUE,
      sort = FALSE
    )
    res <- rowData(sumExp)[["fitDTUResult_Contrast1"]]
    DF_SatuRn[[cl]] = data.frame(empirical_pval = res$empirical_pval, 
                                 transcript_id = rownames(res),
                                 cluster_id = clusters[cl])
  })
  ############################################################
  # save results:
  ############################################################
}

name = paste0("results_batch_SatuRnSC.RData")

full_name = file.path("../07_results", name)

save(DF_SatuRn, TIMES, file = full_name)
