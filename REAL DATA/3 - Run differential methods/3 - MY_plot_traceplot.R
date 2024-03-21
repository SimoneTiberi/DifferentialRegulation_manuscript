MY_plot_traceplot = function(results,
                          gene_id,
                          cluster_id){
  sel_cluster = which( names(results$MCMC_U) == cluster_id)

  sel_gene = which( results$MCMC_U[[sel_cluster]]$Gene_id == gene_id)
  
  group_names = names(results$MCMC_U[[sel_cluster]])[1:2]
  n_iter = nrow(results$MCMC_U[[sel_cluster]][[1]])
  
  DF = data.frame(pi_U_A = results$MCMC_U[[sel_cluster]][[1]][,sel_gene], 
                  pi_U_B = results$MCMC_U[[sel_cluster]][[2]][,sel_gene], 
                  MCMC_iterations = 1:n_iter)
  
  # plot vertical line at burn-in
  sel_cluster_convergence = which( results$Convergence_results$Cluster_id == cluster_id)
  burn_in = results$Convergence_results$burn_in[sel_cluster_convergence]
  
  # Plot the estimated average proportions of each groups:
  ggp1 = ggplot() +
    geom_line(data = DF, aes_string(x = "MCMC_iterations", y = "pi_U_A")) +
    theme_bw() + 
    ylim(0, 1) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text = element_text(size=16), 
          axis.title = element_text(size=14, face="bold"), 
          plot.title = element_text(size=16)) +
    ggtitle(paste(gene_id, "-", cluster_id, "-", group_names[1])) +
    ylab(expression(pi[U])) + 
    geom_vline(xintercept = burn_in, linetype="dashed", 
               colour = "darkgrey")
  
  ggp2 = ggplot() +
    geom_line(data = DF, aes_string(x = "MCMC_iterations", y = "pi_U_B")) +
    theme_bw() + 
    ylim(0, 1) +
    theme(axis.text.x = element_text(vjust = 0.5), 
          axis.text = element_text(size=16), 
          axis.title = element_text(size=14, face="bold"), 
          plot.title = element_text(size=16)) +
    ggtitle(paste(gene_id, "-", cluster_id, "-", group_names[2])) +
    xlab("MCMC iteration") +
    ylab(expression(pi[U])) + 
    geom_vline(xintercept = burn_in, linetype="dashed", 
               colour = "darkgrey")
  
  ggpubr::ggarrange(ggp1, ggp2, ncol = 1, nrow = 2)
}