setwd("~/Desktop/Differential Velocity/FINAL SCRIPTS/REAL DATA - BULK")

# clean enviroment
rm(list = ls())

library(DifferentialRegulation)

############################################################
# load results:
############################################################
# DifferentialRegulation:
name = paste("04_results/DifferentialRegulation_WITH_trace.RData")
load(name)
DR = res_EC$Differential_results

n_res = nrow(DR); n_res
res_names = DR$Transcript_id

if(FALSE){
  # plot randomly selected results:
  set.seed(169612)
  sel = sample.int(n_res, 20)
  for(i in sel){
    plot_bulk_traceplot(res_EC, res_names[i])
  }
}

library(ggplot2)
# plot top results (according to p-val):
save_path = file.path("5 - plots")
A = list()
for(i in 1:20){
  A[[i]] = MY_plot_bulk_traceplot(results =  res_EC, transcript_id = res_names[i])
}

AA_1 = ggpubr::ggarrange( A[[1]], A[[2]], A[[3]], A[[4]], A[[5]], A[[6]], A[[7]], A[[8]], A[[9]], A[[10]],
                        ncol = 2, nrow = 5 )
AA_2 = ggpubr::ggarrange( A[[11]], A[[12]], A[[13]], A[[14]], A[[15]], A[[16]], A[[17]], A[[18]], A[[19]], A[[20]], 
                          ncol = 2, nrow = 5 )

ggsave(filename = "trace_1.pdf",
       plot = AA_1,
       device = "pdf",
       path = save_path,
       width = 15,
       height = 21,
       units = "in",
       dpi = 300,
       limitsize = TRUE)

ggsave(filename = "trace_2.pdf",
       plot = AA_2,
       device = "pdf",
       path = save_path,
       width = 15,
       height = 21,
       units = "in",
       dpi = 300,
       limitsize = TRUE)

# Heidelberg and Welch stationarity test:
res_A = res_EC$MCMC_U[[1]]
res_B = res_EC$MCMC_U[[2]]
R = nrow(res_A)
by. = 100
n_genes = ncol(res_A)
conv_A = rep(NA, n_genes)
conv_B = rep(NA, n_genes)
for(gene in 1:n_genes){
  conv_A[gene] = #heidel.diag(res_A[,gene], pvalue=0.05)[1]
    DifferentialRegulation:::my_heidel_diag(res_A[,gene], R = R, by. = by., pvalue = 0.05)[1]
  conv_B[gene] = #heidel.diag(res_B[,gene], pvalue=0.05)[1]
    DifferentialRegulation:::my_heidel_diag(res_B[,gene], R = R, by. = by., pvalue = 0.05)[1]
}

mean(conv_A)
mean(conv_B)
mean( c(conv_A, conv_B) )
# 0.9583182
# 0.9397038 via heidel.diag

# try for significant results only!
# Heidelberg and Welch stationarity test:
names = DR$Transcript_id[ DR$`Prob-B-UP` > 0.99 | DR$`Prob-B-UP` < 0.01 ]

library(coda)
SEL = res_EC$MCMC_U$Transcript_id %in% names
res_A = res_EC$MCMC_U[[1]][,SEL]
res_B = res_EC$MCMC_U[[2]][,SEL]
R = nrow(res_A)
by. = 100
n_genes = ncol(res_A)
conv_A = rep(NA, n_genes)
conv_B = rep(NA, n_genes)
for(gene in 1:n_genes){
  conv_A[gene] = #heidel.diag(res_A[,gene], pvalue=0.05)[1]
    DifferentialRegulation:::my_heidel_diag(res_A[,gene], R = R, by. = by., pvalue = 0.05)[1]
  conv_B[gene] = #heidel.diag(res_B[,gene], pvalue=0.05)[1]
    DifferentialRegulation:::my_heidel_diag(res_B[,gene], R = R, by. = by., pvalue = 0.05)[1]
}

mean(conv_A)
mean(conv_B)
mean( c(conv_A, conv_B) )
# 0.9568966
# 0.9482759 via heidel.diag

############################################################
# Coherency across runs:
############################################################
# DifferentialRegulation:
name = paste("04_results/DifferentialRegulation_WITH_trace.RData")
load(name)
res_EC_1 = res_EC; rm(res_EC)

name = paste("04_results/DifferentialRegulation_WITH_trace_seed123.RData")
load(name)
res_EC_2 = res_EC; rm(res_EC)

# group A:
DR_1 = res_EC_1$Differential_results[,c(1,6)]
colnames(DR_1)[2] = "run_1"
DR_2 = res_EC_2$Differential_results[,c(1,6)]
colnames(DR_2)[2] = "run_2"
DR_A = merge(DR_1, DR_2, by = "Transcript_id")

# group B:
DR_1 = res_EC_1$Differential_results[,c(1,8)]
colnames(DR_1)[2] = "run_1"
DR_2 = res_EC_2$Differential_results[,c(1,8)]
colnames(DR_2)[2] = "run_2"
DR_B = merge(DR_1, DR_2, by = "Transcript_id")

# combine groups A and B:
DR = rbind(DR_A, DR_B)
plot(DR$run_1, DR$run_2)
cor(DR$run_1, DR$run_2)
# 0.9870339
