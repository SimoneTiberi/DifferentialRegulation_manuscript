alias R='R_LIBS=/home/stiber/R/x86_64-pc-linux-gnu-library/4.3 /home/Shared/simone/DifferentialRegulation/software/R-4.3.0/bin/R'

R

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# Single-cell SIMULATION:
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
rm(list = ls())
# setwd("/home/Shared/simone/Diff_Velo/NEW_simulation")
setwd("~/Desktop/Differential Velocity/FINAL SCRIPTS/SINGLE-CELL SIMULATION")

DF = c()

# DifferentialRegulation
save_name = paste0("07_results/DifferentialRegulation_DGE_2000_10.RData")
load(save_name)
TIME = do.call(rbind, TIMES)
TIME_DGE = sum(TIME[,3])
rm(TIME); rm(TIMES)

save_name = paste0("07_results/DifferentialRegulation_NO_DGE_2000_10.RData")
load(save_name)
TIME = do.call(rbind, TIMES)
TIME = sum(TIME[,3])
rm(TIMES)

DF$DifferentialRegulation = (TIME + TIME_DGE)/2

rm(TIME); rm(TIME_DGE)

# eisaR:
save_name = paste0("07_results/results_DGE_eisaR.RData")
load(save_name)
TIME = do.call(rbind, TIMES)
TIME_DGE = sum(TIME[,3])
rm(TIME); rm(TIMES)

save_name = paste0("07_results/results_NO_DGE_eisaR.RData")
load(save_name)
TIME = do.call(rbind, TIMES)
TIME = sum(TIME[,3])
rm(TIMES)

DF$eisaR = (TIME + TIME_DGE)/2

rm(TIME); rm(TIME_DGE)

# BRIE2 50:50 without sample_id
save_name = paste0("07_results/BRIE2_USA_DGE.RData")
load(save_name)
TIME_DGE = TIME[3]
rm(TIME);

save_name = paste0("07_results/BRIE2_USA.RData")
load(save_name)
TIME = TIME[3]

DF$BRIE2 = (TIME + TIME_DGE)/2

rm(TIME); rm(TIME_DGE)

# DEXSeq
save_name = paste0("07_results/results_DGE_DEXSeq_USA.RData")
load(save_name)
TIME = do.call(rbind, TIMES)
TIME_DGE = sum(TIME[,3])
rm(TIME); rm(TIMES)

save_name = paste0("07_results/results_NO_DGE_DEXSeq_USA.RData")
load(save_name)
TIME = do.call(rbind, TIMES)
TIME = sum(TIME[,3])
rm(TIMES)

DF$DEXSeq = (TIME + TIME_DGE)/2

rm(TIME); rm(TIME_DGE)

# SatuRn
load("07_results/results_DGE_SatuRn.RData")
TIME = do.call(rbind, TIMES)
TIME_DGE = sum(TIME[,3])
rm(TIME); rm(TIMES)

load("07_results/results_NO_DGE_SatuRn.RData")
TIME = do.call(rbind, TIMES)
TIME = sum(TIME[,3])
rm(TIMES)

DF$satuRn = (TIME + TIME_DGE)/2

rm(TIME); rm(TIME_DGE)

# DRIMSeq
load("07_results/results_DGE_DRIMSeq.RData")
TIME = do.call(rbind, TIMES)
TIME_DGE = sum(TIME[,3])
rm(TIME); rm(TIMES)

load("07_results/results_NO_DGE_DRIMSeq.RData")
TIME = do.call(rbind, TIMES)
TIME = sum(TIME[,3])
rm(TIMES)

DF$DRIMSeq = (TIME + TIME_DGE)/2

rm(TIME); rm(TIME_DGE)

# average time:
DF = unlist(DF)

# change from seconds to minutes:
DF = DF/60
DF

# plot time in DF
gg_data_sc = data.frame(Method_sc = names(DF),
                        Time_Minutes = c(DF),
                        Colour_sc = c("#1F78B4", # DifferentialRegulation
                                      "#A6761D",  # eisaR
                                      "#DE77AE", # BRIE2
                                      "#A63603",  # "DEXSeq"
                                      "#404040", # "satuRn"
                                      "#66A61E" # "DRIMSeq"
                        ))

gg_data_sc$Time_sc <- round(gg_data_sc$Time_Minutes)
gg_data_sc$Method_sc <- factor(gg_data_sc$Method_sc)
## Main plots
library(ggplot2)
library(dplyr)

Methods <- arrange(gg_data_sc, Time_Minutes)$Method_sc
gg_sc = ggplot(aes(x = Method_sc, Time_Minutes, fill = Colour_sc), data = gg_data_sc) +
  geom_bar( stat = "identity",position = "stack",width = 0.7) + 
  theme_bw() +   xlab("") +  ylab("Time (Minutes)") + 
  geom_text(aes(y = Time_Minutes, label=Time_sc), vjust=-0.5, 
            color="black", size=3.5) +
  coord_trans(y = "sqrt") +
  scale_x_discrete(limits = Methods) +
  scale_y_continuous(breaks = c(0,15,60,120,500,1000,1500),expand = c(0,0),
                     limits = c(0, 2000))+
  scale_fill_identity() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=rel(1)),   
        axis.text.y= element_text(size=rel(1),angle = 0, hjust = 1),  
        axis.title = element_text(size=rel(1)),      
        panel.grid.minor = element_blank(),         
        panel.grid.major.x = element_blank(),       
        panel.grid.major.y = element_line(size = 0.2, color = "grey", linetype = 2),        
        aspect.ratio = 1,        
        legend.position = c(0.2, 0.85),
        legend.text=element_text(size=10)
  ) 
gg_sc

save(gg_sc, file = "08_plots/gg_sc.RData")

ggsave(filename = "Barplot_Computational_Cost_sc_Simulation.pdf",
       plot = gg_sc,
       path = file.path("08_plots"),
       device = "pdf",
       width = 5,
       height = 5,
       units = "in",
       dpi = 300,
       limitsize = TRUE)

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# BULK SIMULATION:
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
rm(list = ls())

library(data.table)
# library(iCOBRA)
library(ggplot2)
library(dplyr)

# setwd("/home/Shared/simone/DifferentialRegulation")
setwd("~/Desktop/Differential Velocity/FINAL SCRIPTS/BULK SIMULATION")
source("2 - R scripts/2 - plot theme.R")

# directories:
directories = c("de0_ds0_dr0",
                "de0_ds0_dr2000",
                "de2000_ds0_dr2000",
                "de0_ds2000_dr2000")

directories_results = c("NULL",
                        "DR",
                        "DR + DGE",
                        "DR + DAS")

# load TIME:
# DEXSeq - ECCs:
name = paste0("TIME_DEXSeq_ECCs.RData")
file = file.path("3 - results",name)
if(file.exists(file)){
  load(file)
  DEXSeq_ECs = do.call(cbind, lapply(TIME,"[",3))
  TIME_merged = DEXSeq_ECs
  rm(DEXSeq_ECs); rm(TIME)
}

DEXSeq_ECs_py <- c(2*60+29.880, 3*60+10.513, 2*60+19.612, 2*60+5.829)
TIME_merged = rbind(TIME_merged, DEXSeq_ECs_py)

# DEXSeq - TECs:
name = paste0("TIME_DEXSeq_TECs.RData")
file = file.path("3 - results",name)
if(file.exists(file)){
  load(file)
  DEXSeq_TECs = do.call(cbind, lapply(TIME,"[",3))
  TIME_merged = rbind(TIME_merged, DEXSeq_TECs)
  rm(DEXSeq_TECs); rm(TIME)
}

# eisaR:
name = paste0("TIME_eisaR.RData")
file = file.path("3 - results",name)
if(file.exists(file)){
  load(file)
  eisaR = do.call(cbind, lapply(TIME,"[",3))
  TIME_merged = rbind(TIME_merged, eisaR)
  rm(eisaR); rm(TIME)
}

# DRIMSeq:
name = paste0("TIME_DRIMSeq.RData")
file = file.path("3 - results",name)
if(file.exists(file)){
  load(file)
  DRIMSeq = do.call(cbind, lapply(TIME,"[",3))
  TIME_merged = rbind(TIME_merged, DRIMSeq)
  rm(DRIMSeq); rm(TIME)
}

# DifferentialRegulation:
name = paste0("TIME_DifferentialRegulation.RData")
file = file.path("3 - results",name)
if(file.exists(file)){
  load(file)
  DR = do.call(cbind, lapply(TIME,"[",3))
  TIME_merged = rbind(TIME_merged, DR)
  rm(DR); rm(TIME)
}

# satuRn:
name = paste0("TIME_satuRn.RData")
file = file.path("3 - results",name)
if(file.exists(file)){
  load(file)
  satuRn = do.call(cbind, lapply(TIME,"[",3))
  TIME_merged = rbind(TIME_merged, satuRn)
  rm(satuRn); rm(TIME)
}

# SUPPA2:
SUPPA2_py1 <- rep(39.662, 4)
TIME_merged = rbind(TIME_merged, SUPPA2_py1)

name = paste0("TIME_R1_SUPPA2.RData")
file = file.path("3 - results",name)
if(file.exists(file)){
  load(file)
  SUPPA2 = do.call(cbind, lapply(TIME,"[",3))
  TIME_merged = rbind(TIME_merged, SUPPA2)
  rm(SUPPA2); rm(TIME)
}

t1 <- c(1*60+3.701 + 	0*1+59.536 + 28*60+33.114)
t2 <- c(0*60+55.931 + 0*60+51.392 + 20*60+25.862)
t3 <- c(0*60+53.825 + 0*60+58.392  + 17*60+46.103)
t4 <- c(0*60+55.115 + 1*60+1.240 + 19*60+53.248)
SUPPA2_py2 <- c(t1, t2, t3, t4)
TIME_merged = rbind(TIME_merged, SUPPA2_py2)

name = paste0("TIME_R2_SUPPA2.RData")
file = file.path("3 - results",name)
if(file.exists(file)){
  load(file)
  SUPPA2 = do.call(cbind, lapply(TIME,"[",3))
  TIME_merged = rbind(TIME_merged, SUPPA2)
  rm(SUPPA2); rm(TIME)
}

rownames(TIME_merged) <- c("DEXSeq_ECs_r", "DEXSeq_ECs_py",
                           "DEXSeq_TECs", "eisaR",
                           "DRIMSeq",
                           "DifferentialRegulation", "satuRn",
                           "SUPPA2_py1", "SUPPA2_r1",
                           "SUPPA2_py2", "SUPPA2_r2")
colnames(TIME_merged) <- directories_results
DEXSeq_ECs <- colSums(TIME_merged[c("DEXSeq_ECs_r","DEXSeq_ECs_py"),])
SUPPA2 <- colSums(TIME_merged[c("SUPPA2_py1", "SUPPA2_r1",
                                "SUPPA2_py2", "SUPPA2_r2"),])
TIME_merged = rbind(TIME_merged, DEXSeq_ECs, SUPPA2)
Average <- rowSums(TIME_merged[,2:4])/3
df <- t(as.data.frame(Average))
df <- reshape::melt(df)
df$X2 <- as.factor(df$X2)
df$value <- df$value/60
gg_data <- df[-c(1)]
colnames(gg_data) <- c("Method", "Time_Minutes")
Time <- round(gg_data$Time_Minutes)
gg_data$Method <- factor(gg_data$Method, levels=c(methods_all), labels=c(methods_all))
gg_data <- na.omit(gg_data)
## Main plots
colors_method <- as.data.frame(colors_method)
colors_method$Method <- rownames(colors_method)
gg_data <- merge(gg_data,colors_method,by = "Method")
gg_data$Time_bulk <- round(gg_data$Time_Minutes)
gg_data$label_ypos = gg_data$Time_Minutes

`%notin%` <- Negate(`%in%`)
Methods <- arrange(gg_data, Time_Minutes)$Method
(gg_bulk <- gg_data %>%
    ggplot(aes(x = Method, Time_Minutes, fill = colors_method)) +
    geom_bar( stat = "identity",position = "stack",width = 0.7) + 
    theme_bw() +   xlab("") +  ylab("Time (Minutes)") + 
    geom_text(aes(y = label_ypos, label=Time_bulk), vjust=-0.5, 
              color="black", size=3.5)+
    coord_trans(y = "sqrt") +
    scale_x_discrete(limits = Methods) +
    scale_y_continuous(breaks = c(5,15,30,60),
                       limits = c(0, 70))+
    scale_fill_identity() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size=rel(1)),   
          axis.text.y= element_text(size=rel(1),angle = 0, hjust = 1),  
          axis.title = element_text(size=rel(1)),      
          panel.grid.minor = element_blank(),         
          panel.grid.major.x = element_blank(),       
          panel.grid.major.y = element_line(size = 0.2, color = "grey", linetype = 2),        
          aspect.ratio = 1,        
          legend.position = c(0.2, 0.85),
          legend.text=element_text(size=10)
    ) )
gg_bulk

save(gg_bulk, file = "5 - plots/gg_bulk.RData")

ggsave(filename = "Barplot_Computational_Cost.pdf",
       plot = gg_bulk,
       path = file.path("5 - plots"),
       device = "pdf",
       width = 5,
       height = 5,
       units = "in",
       dpi = 300,
       limitsize = TRUE)

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# ARRANGE 2 plots:
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
rm(list =ls())

setwd("~/Desktop/Differential Velocity/FINAL SCRIPTS")

load("SINGLE-CELL SIMULATION/08_plots/gg_sc.RData")
load("BULK SIMULATION/5 - plots/gg_bulk.RData")

AA = egg::ggarrange( plots = list(gg_bulk + labs(title = "bulk"),
                                   gg_sc + ylab("") + labs(title = "single-cell")),
                      ncol = 2, nrow = 1)
# sometimes error with gg_sc -> just re-run it.

ggsave(filename = "Runtime_BOTH.pdf",
       plot = AA,
       device = "pdf",
       width = 9,
       height = 5,
       units = "in",
       dpi = 300,
       limitsize = TRUE)

