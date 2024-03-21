MY_plot_bulk_traceplot = function(results,
                                  transcript_id){
  sel = which(results$MCMC_U$Transcript_id == transcript_id)
  
  group_names = names(results$MCMC_U)[1:2]
  n_iter = nrow(results$MCMC_U[[1]])
  
  DF = data.frame(pi_U_A = results$MCMC_U[[1]][,sel], 
                  pi_U_B = results$MCMC_U[[2]][,sel], 
                  MCMC_iterations = 1:n_iter)
  
  # plot vertical line at burn-in
  burn_in = results$Convergence_results$burn_in
  
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
    ggtitle(paste(transcript_id, " - ", group_names[1])) +
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
    ggtitle(paste(transcript_id, " - ", group_names[2])) +
    xlab("MCMC iteration") +
    ylab(expression(pi[U])) + 
    geom_vline(xintercept = burn_in, linetype="dashed", 
               colour = "darkgrey")
  
  ggpubr::ggarrange(ggp1, ggp2, ncol = 1, nrow = 2)
}