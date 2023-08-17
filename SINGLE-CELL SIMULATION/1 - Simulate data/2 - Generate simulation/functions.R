# load libraries
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(doParallel)
  library(doRNG)
  library(Matrix)
})

# simulation of mouse data
mouse_simulation <- function (sce, CLUSTERS, 
                              p_genes) {
  SCES = list()
  TRUTH = list()
  for(i in 1:length(CLUSTERS)){
    # select cluster
    temp <- sce[, sce$cell_type == CLUSTERS[[i]]]
    
    # remove lowly expressed genes: at least 10 non-zero cells:
    filter <- (rowSums(assays(temp)$TOT_counts[, temp$group == "A"] > 0) >= 10) & (rowSums(assays(temp)$TOT_counts[, temp$group == "B"] > 0) >= 10)
    
    # filter min distince in pi_S between groups (0.2 minimum):
    sel_A = temp$group == "A"
    tot_counts_A = rowSums(assays(temp)$spliced[,sel_A] + assays(temp)$unspliced[,sel_A])
    tot_counts_B = rowSums(assays(temp)$spliced[,!sel_A] + assays(temp)$unspliced[,!sel_A])
    
    tot_S_A = rowSums(assays(temp)$spliced[,sel_A])
    tot_S_B = rowSums(assays(temp)$spliced[,!sel_A])
    
    pi_S_A = tot_S_A/tot_counts_A
    pi_U_B = 1 - tot_S_B/tot_counts_B
    difference = abs( pi_S_A - pi_U_B )
    sel = (difference > 0.2)
    sel[is.na(sel)] = FALSE
    
    filter = filter & sel
    
    # randomly draw proportion of genes and store ground truth
    genes <- sample(rownames(temp)[filter], size = p_genes*length(filter), replace = FALSE)
    truth <- data.frame(Gene_id = rownames(temp), Cell_type = CLUSTERS[[i]],
                        truth = ifelse(rownames(temp) %in% genes, 1, 0))
    
    for(gene in genes){
      # sample the group to swap (50:50 prob for the 2 groups):
      if(rbinom(1,1,0.5) == 1){
        sel_cells = temp$group == "A"
      }else{
        sel_cells = temp$group == "B"
      }
      sel_gene = which(rownames(temp) == gene)
      
      # swap S and U counts:
      spliced = assays(temp)$spliced[sel_gene, sel_cells]
      unspliced = assays(temp)$unspliced[sel_gene, sel_cells]
      assays(temp)$spliced[sel_gene, sel_cells] = unspliced
      assays(temp)$unspliced[sel_gene, sel_cells] = spliced
    }
    SCES[[i]] = temp
    TRUTH[[i]] = truth
    print(i)
  }
  
  SCE <- do.call("cbind", SCES)
  metadata(SCE)$truth <- do.call("rbind", TRUTH)
  
  return(SCE)
}

add_DGE <- function (sce, CLUSTERS,  
                     p_genes_DGE) {
  # store truth from first simulation
  truth_DA <- metadata(sce)$truth
  metadata(sce)$truth <- NULL
  
  SCES = list()
  TRUTH = list()
  
  for(i in 1:length(CLUSTERS)){
    # select cluster
    temp <- sce[, sce$cell_type == CLUSTERS[[i]]]
    
    # remove lowly expressed genes: at least 10 non-zero cells:
    filter <- (rowSums(assays(temp)$TOT_counts[, temp$group == "A"] > 0) >= 10) & (rowSums(assays(temp)$TOT_counts[, temp$group == "B"] > 0) >= 10)
    
    # randomly draw proportion of genes and store ground truth
    genes <- sample(rownames(temp)[filter], size = p_genes_DGE*length(filter), replace = FALSE)
    truth <- data.frame(Gene_id = rownames(temp), Cell_type = CLUSTERS[[i]],
                        truth_DGE = ifelse(rownames(temp) %in% genes, 1, 0),
                        FC = 1)
    
    for(gene in genes){
      # sample the group to swap (50:50 prob for the 2 groups):
      if(rbinom(1,1,0.5) == 1){
        sel_cells = temp$group == "A"
      }else{
        sel_cells = temp$group == "B"
      }
      sel_gene = which(rownames(temp) == gene)
      
      # simulate a random FC, with mean 3 and variance 1:
      FC = 2 + rexp(1, rate = 1)

      # store FC in truth matrix:
      truth$FC[sel_gene] = FC 
      
      # swap S and U counts:
      assays(temp)$spliced[sel_gene, sel_cells] = FC * assays(temp)$spliced[sel_gene, sel_cells]
      assays(temp)$unspliced[sel_gene, sel_cells] = FC * assays(temp)$unspliced[sel_gene, sel_cells]
    }
    SCES[[i]] = temp
    TRUTH[[i]] = truth
    print(i)
  }

  SCE <- do.call("cbind", SCES)
  truth_DGE <- do.call("rbind", TRUTH)
  metadata(SCE)$truth <- merge(truth_DA, truth_DGE, by = c("Gene_id", "Cell_type"), all.x = TRUE)
  
  return(SCE)
}