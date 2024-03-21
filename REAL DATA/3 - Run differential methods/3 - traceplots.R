setwd("~/Desktop/Differential Velocity/FINAL SCRIPTS/REAL DATA")
# clean enviroment
rm(list = ls())

library(DifferentialRegulation)

############################################################
# load results:
############################################################
load("02_results/DifferentialRegulation_WITH_trace.RData")
DR = results_EC$Differential_results

n_res = nrow(DR); n_res
res_names = DR$Gene_id
cl_names = DR$Cluster_id

# plot randomly selected results:
if(FALSE){
  set.seed(169612)
  sel = sample.int(n_res, 20)
  for(i in sel){
    plot_traceplot(results_EC, res_names[i], cl_names[i])
  }
}

library(ggplot2)
# plot top results (according to p-val):
save_path = file.path("5 - plots")
A = list()
for(i in 1:20){
  A[[i]] = MY_plot_traceplot(results_EC, res_names[i], cl_names[i])
}

AA_1 = ggpubr::ggarrange( A[[1]], A[[2]], A[[3]], A[[4]], A[[5]], A[[6]], A[[7]], A[[8]], A[[9]], A[[10]],
                          ncol = 2, nrow = 5 )
AA_2 = ggpubr::ggarrange( A[[11]], A[[12]], A[[13]], A[[14]], A[[15]], A[[16]], A[[17]], A[[18]], A[[19]], A[[20]], 
                          ncol = 2, nrow = 5 )

ggsave(filename = "trace_sc_1.pdf",
       plot = AA_1,
       device = "pdf",
       path = save_path,
       width = 15,
       height = 21,
       units = "in",
       dpi = 300,
       limitsize = TRUE)

ggsave(filename = "trace_sc_2.pdf",
       plot = AA_2,
       device = "pdf",
       path = save_path,
       width = 15,
       height = 21,
       units = "in",
       dpi = 300,
       limitsize = TRUE)




# Heidelberg and Welch stationarity test:
CONV_A = CONV_B = list()
for(cl in 1:length(results_EC$MCMC_U)){
  res_one_cl_A = results_EC$MCMC_U[[cl]][[1]]
  res_one_cl_B = results_EC$MCMC_U[[cl]][[2]]
  R = nrow(res_one_cl_A)
  by. = 100
  n_genes = ncol(res_one_cl_A)
  conv_A = rep(NA, n_genes)
  conv_B = rep(NA, n_genes)
  for(gene in 1:n_genes){
    conv_A[gene] = heidel.diag(res_one_cl_A[,gene], pvalue=0.05)
    #DifferentialRegulation:::my_heidel_diag(res_one_cl_A[,gene], R = R, by. = by., pvalue = 0.05)[1]
    conv_B[gene] = heidel.diag(res_one_cl_B[,gene], pvalue=0.05)
    #DifferentialRegulation:::my_heidel_diag(res_one_cl_B[,gene], R = R, by. = by., pvalue = 0.05)[1]
  }
  CONV_A[[cl]] = conv_A
  CONV_B[[cl]] = conv_B
  
  print(cl)
  print(mean(conv_A))
  print(mean(conv_B))
}

mean(c(unlist(CONV_A),unlist(CONV_B)))
# 0.9830507
# 0.9695782 via heidel.diag


library(coda)
# Heidelberg and Welch stationarity test:
CONV_A = CONV_B = list()
for(cl in 1:length(results_EC$MCMC_U)){
  cluster = names(results_EC$MCMC_U)[cl]
  sel_cluster = results_EC$Differential_results$Cluster_id == cluster
  DR = results_EC$Differential_results[sel_cluster,]
  
  names = DR$Gene_id[ DR$`Prob-6mon-UP` > 0.99 | DR$`Prob-6mon-UP` < 0.01 ]
  
  SEL = results_EC$MCMC_U[[cl]]$Gene_id %in% names
  
  res_one_cl_A = results_EC$MCMC_U[[cl]][[1]][,SEL]
  res_one_cl_B = results_EC$MCMC_U[[cl]][[2]][,SEL]
  R = nrow(res_one_cl_A)
  by. = 100
  n_genes = ncol(res_one_cl_A)
  conv_A = rep(NA, n_genes)
  conv_B = rep(NA, n_genes)
  for(gene in 1:n_genes){
    conv_A[gene] = #heidel.diag(res_one_cl_A[,gene], pvalue=0.05)[1]
      DifferentialRegulation:::my_heidel_diag(res_one_cl_A[,gene], R = R, by. = by., pvalue = 0.05)[1]
    conv_B[gene] =  #heidel.diag(res_one_cl_B[,gene], pvalue=0.05)[1]
      DifferentialRegulation:::my_heidel_diag(res_one_cl_B[,gene], R = R, by. = by., pvalue = 0.05)[1]
  }
  CONV_A[[cl]] = conv_A
  CONV_B[[cl]] = conv_B
  
  print(cl)
  print(mean(conv_A))
  print(mean(conv_B))
}

mean(c(unlist(CONV_A),unlist(CONV_B)))
# 0.9868184
# 0.975172 via heidel.diag

############################################################
# Coherency across runs:
############################################################
# DifferentialRegulation:
name = paste("02_results/DifferentialRegulation_WITH_trace.RData")
load(name)
res_EC_1 = results_EC$US_results; rm(results_EC)

name = paste("02_results/DifferentialRegulation_WITH_trace_seed123.RData")
load(name)
res_EC_2 = results_EC$US_results; rm(results_EC)

# group A:
DR_1 = res_EC_1[,c(1,2,8)]
colnames(DR_1)[3] = "run_1"
DR_2 = res_EC_2[,c(1,2,8)]
colnames(DR_2)[3] = "run_2"
DR_A = merge(DR_1, DR_2, by = c("Gene_id", "Cluster_id"))

# group B:
DR_1 = res_EC_1[,c(1,2,10)]
colnames(DR_1)[3] = "run_1"
DR_2 = res_EC_2[,c(1,2,10)]
colnames(DR_2)[3] = "run_2"
DR_B = merge(DR_1, DR_2, by = c("Gene_id", "Cluster_id"))

# combine groups A and B:
DR = rbind(DR_A, DR_B)
plot(DR$run_1, DR$run_2)
cor(DR$run_1, DR$run_2)
# 0.9959115
