#This script shows the distribution of the exact statistic under the null hypothesis for multiple covariates
rm(list = ls())
set.seed(123)
options(scipen=999)

#Working directory
if(rstudioapi::isAvailable()) setwd(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))

#Libraries
library(mvtnorm)
library(pbmcapply)
library(gridExtra)
library(dplyr)

#Import functions
source("Package/MSEmin_package v1.0.R")
source("Functions/Functions Simulations.R")
source("Functions/Functions Analyze MSE.R")
source("Functions/Functions Analyze Multivariate Exact Test.R")
source("Functions/Functions Simulations PostSelection.R")
#Modify the effect size to change the situation
EffectSizeMSE=1
source("DGP/DGP_4.R")

#Number of simulations
S=150000
alpha=0.10
R=200 #Number of draws of priors

#Note: We need very large S to compare distribution when the Rejection Rate is very low. 
#It is true when EffectSizeMSE is far from the boundary=1 

set.seed(123)
simRes=pbmclapply(1:S,function(i){
  set.seed(i)
  main_data_funct=generateDatasetMulti(N_funct=N, EStreatment_funct = beta,
                                       EScontrol_funct = ESvec, CorrMatrix_funct = CorrMat, 
                                       sigma2_funct=sigma2,
                                       SD_vec_funct = SD_vec)
  pilot_data_funct=generateDatasetMulti(N_funct=N, 
                                        EStreatment_funct = beta,
                                        EScontrol_funct = ESvec, 
                                        CorrMatrix_funct = diag(4), 
                                        sigma2_funct=sigma2,
                                        SD_vec_funct = SD_vec)
  cbind(i,bind_rows(computeStatsSelectionRules(main_data_funct=main_data_funct,
                                                               #selectionRule=c("Final Sample","Pilot","Point Prior","Distrib Prior"),
                                                               selectionRule=c("Final Sample","Pilot","Point Prior"),
                                                               pilot_df=pilot_data_funct, 
                                                               priors_mu=ESvec,
                                                               priors_var=rep(0.00001,3),
                                                               alpha_funct=alpha,
                                                               R_funct=R)))
},  mc.cores = min(detectCores()-1,20))
simRes_DF=do.call(rbind.data.frame, simRes)
head(simRes_DF)

mean(simRes_DF[simRes_DF$cm == "y ~ treatment",]$correction_unscaled)
mean(simRes_DF[simRes_DF$cm == "y ~ treatment",]$nu)
mean(simRes_DF[simRes_DF$cm == "y ~ treatment + X1 + X2 + X3",]$nu)

#Select Model 0 and Full Model  
simRes_DF = simRes_DF %>% filter(cm == "y ~ treatment" | cm == "y ~ treatment + X1 + X2 + X3")

#Select Model 0 only
simRes_DF_NullModelOnly = simRes_DF %>% filter(cm == "y ~ treatment")

hist(simRes_DF_NullModelOnly$pval_finalSampleBased)
mean(simRes_DF_NullModelOnly$pval_finalSampleBased<0.10)

hist(simRes_DF_NullModelOnly$pval_pilotBased)
mean(simRes_DF_NullModelOnly$pval_pilotBased<0.10)

hist(simRes_DF_NullModelOnly$pval_pointPriorBased)
mean(simRes_DF_NullModelOnly$pval_pointPriorBased<0.10)

#Matrix of Rejection Rate
RR=matrix(data=NA,nrow=3,ncol=3)
RR[,1]=c(mean(simRes_DF_NullModelOnly$pval_finalSampleBased<0.01),
          mean(simRes_DF_NullModelOnly$pval_finalSampleBased<0.05),
          mean(simRes_DF_NullModelOnly$pval_finalSampleBased<0.10))
RR[,2]=c(mean(simRes_DF_NullModelOnly$pval_pilotBased<0.01),
          mean(simRes_DF_NullModelOnly$pval_pilotBased<0.05),
          mean(simRes_DF_NullModelOnly$pval_pilotBased<0.10))
RR[,3]=c(mean(simRes_DF_NullModelOnly$pval_pointPriorBased<0.01),
          mean(simRes_DF_NullModelOnly$pval_pointPriorBased<0.05),
          mean(simRes_DF_NullModelOnly$pval_pointPriorBased<0.10))
RR=round(RR,3)
RR

table(simRes_DF_NullModelOnly$pval_pointPriorBased<0.10)

#Generate selection decisions
#Select model with final data: Full model as default
simRes_DF$selected_finalSample = ifelse(simRes_DF$pval_finalSampleBased<=alpha,1,0)
simRes_DF[simRes_DF$cm == "y ~ treatment + X1 + X2 + X3",]$selected_finalSample=1-simRes_DF[simRes_DF$cm == "y ~ treatment",]$selected_finalSample
#Select model with pilot data: Full model as default
simRes_DF$selected_pilotSample = ifelse(simRes_DF$pval_pilotBased<=alpha,1,0)
simRes_DF[simRes_DF$cm == "y ~ treatment + X1 + X2 + X3",]$selected_pilotSample=1-simRes_DF[simRes_DF$cm == "y ~ treatment",]$selected_pilotSample
#Select model with registered priors: Full model as default
simRes_DF$selected_pointPrior = ifelse(simRes_DF$pval_pointPriorBased<=alpha,1,0)
simRes_DF[simRes_DF$cm == "y ~ treatment + X1 + X2 + X3",]$selected_pointPrior=1-simRes_DF[simRes_DF$cm == "y ~ treatment",]$selected_pointPrior


##Show densities

cov_val_gr1 = cov(
  simRes_DF[simRes_DF$cm == "y ~ treatment", ]$estimate,
  simRes_DF[simRes_DF$cm == "y ~ treatment", ]$selected_finalSample
)

cov_lab_gr1 <- paste0(
  "atop(",
  "Cov(widehat(beta)[m], 1[m]), ",
  "'=' * ", round(cov_val_gr1, 5),
  ")"
)


gr1=ggplot(simRes_DF[simRes_DF$cm == "y ~ treatment",], 
          aes(x = estimate, fill = as.factor(selected_finalSample))) +
  geom_density(
    data = subset(simRes_DF[simRes_DF$cm == "y ~ treatment",], selected_finalSample == 1),
    alpha=0.4
  ) +
  geom_density(
    data = subset(simRes_DF[simRes_DF$cm == "y ~ treatment",], selected_finalSample == 0),
    aes(y = after_stat(density)),
    alpha=0.4
  ) +
  geom_hline(yintercept = 0, color = "black") +
  annotate(
    "text",
    x = -Inf, y = Inf,
    label = cov_lab_gr1,
    hjust = -0.05, vjust = 1.5,
    size = 4,
    parse = TRUE
  )  +
  labs(
    #title = "Selection with final data estimates",
    x = NULL,
    y = NULL,
    fill = NULL
  ) + 
  scale_fill_discrete(
    labels = c("Not selected", "Selected")
  ) +
  theme_minimal() + 
  theme( axis.text.y = element_blank(),      # drop y-axis numerical values
         axis.ticks.y = element_blank(),     # optional: drop y-axis ticks
         legend.position = "none",          # legend below the graph
         panel.grid = element_blank()
  ) 
gr1

cov_val_gr2 = cov(
  simRes_DF[simRes_DF$cm == "y ~ treatment + X1 + X2 + X3", ]$estimate,
  simRes_DF[simRes_DF$cm == "y ~ treatment + X1 + X2 + X3", ]$selected_finalSample
)

cov_lab_gr2 = paste0(
  "atop(",
  "Cov(widehat(beta)[m], 1[m]), ",
  "'=' * ", round(cov_val_gr2, 5),
  ")"
)


gr2=ggplot(simRes_DF[simRes_DF$cm == "y ~ treatment + X1 + X2 + X3",], 
       aes(x = estimate, fill = as.factor(selected_finalSample))) +
  geom_density(
    data = subset(simRes_DF[simRes_DF$cm == "y ~ treatment + X1 + X2 + X3",], selected_finalSample == 1),
    alpha=0.4
  ) +
  geom_density(
    data = subset(simRes_DF[simRes_DF$cm == "y ~ treatment + X1 + X2 + X3",], selected_finalSample == 0),
    aes(y = after_stat(density)),
    alpha=0.4
  ) +
  geom_hline(yintercept = 0, color = "black") +
  annotate(
    "text",
    x = -Inf, y = Inf,
    label = cov_lab_gr2,
    hjust = -0.05, vjust = 1.5,
    size = 4,
    parse = TRUE
  ) +
  labs(
    #title = "Selection with final data estimates",
    x = NULL,
    y = NULL,
    fill = NULL
  ) + 
  scale_fill_discrete(
    labels = c("Not selected", "Selected")
  ) +
  theme_minimal() + 
  theme( axis.text.y = element_blank(),      # drop y-axis numerical values
          axis.ticks.y = element_blank(),     # optional: drop y-axis ticks
          legend.position = "none",          # legend below the graph
         panel.grid = element_blank()
   ) 
  

cov_val_gr3 = cov(
  simRes_DF[simRes_DF$cm == "y ~ treatment", ]$estimate,
  simRes_DF[simRes_DF$cm == "y ~ treatment", ]$selected_pilotSample
)

cov_lab_gr3 = paste0(
  "atop(",
  "Cov(widehat(beta)[m], 1[m]), ",
  "'=' * ", round(cov_val_gr3, 5),
  ")"
)


gr3=ggplot(simRes_DF[simRes_DF$cm == "y ~ treatment",], 
           aes(x = estimate, fill = as.factor(selected_pilotSample))) +
  geom_density(
    data = subset(simRes_DF[simRes_DF$cm == "y ~ treatment",], selected_pilotSample == 1),
    alpha=0.4
  ) +
  geom_density(
    data = subset(simRes_DF[simRes_DF$cm == "y ~ treatment",], selected_pilotSample == 0),
    aes(y = after_stat(density)),
    alpha=0.4
  ) +
  geom_hline(yintercept = 0, color = "black") +
  annotate(
    "text",
    x = -Inf, y = Inf,
    label = cov_lab_gr3,
    hjust = -0.05, vjust = 1.5,
    size = 4,
    parse = TRUE
  )  +
  labs(
    #title = "Selection with final data estimates",
    x = NULL,
    y = NULL,
    fill = NULL
  ) + 
  scale_fill_discrete(
    labels = c("Not selected", "Selected")
  ) +
  theme_minimal() + 
  theme( axis.text.y = element_blank(),      # drop y-axis numerical values
         axis.ticks.y = element_blank(),     # optional: drop y-axis ticks
         legend.position = "none",          # legend below the graph
         panel.grid = element_blank()
  )  

cov_val_gr4 = cov(
  simRes_DF[simRes_DF$cm == "y ~ treatment + X1 + X2 + X3", ]$estimate,
  simRes_DF[simRes_DF$cm == "y ~ treatment + X1 + X2 + X3", ]$selected_pilotSample
)
#cov_lab_gr4 = paste0("Cov(widehat(beta)[m], 1[m]) == ", round(cov_val_gr4, 5))

cov_lab_gr4 = paste0(
  "atop(",
  "Cov(widehat(beta)[m], 1[m]), ",
  "'=' * ", round(cov_val_gr4, 5),
  ")"
)


gr4=ggplot(simRes_DF[simRes_DF$cm == "y ~ treatment + X1 + X2 + X3",], 
           aes(x = estimate, fill = as.factor(selected_pilotSample))) +
  geom_density(
    data = subset(simRes_DF[simRes_DF$cm == "y ~ treatment + X1 + X2 + X3",], selected_pilotSample == 1),
    alpha=0.4
  ) +
  geom_density(
    data = subset(simRes_DF[simRes_DF$cm == "y ~ treatment + X1 + X2 + X3",], selected_pilotSample == 0),
    aes(y = after_stat(density)),
    alpha=0.4
  ) +
  geom_hline(yintercept = 0, color = "black") +
  annotate(
    "text",
    x = -Inf, y = Inf,
    label = cov_lab_gr4,
    hjust = -0.05, vjust = 1.5,
    size = 4,
    parse = TRUE
  )  +
  labs(
    #title = "Selection with final data estimates",
    x = NULL,
    y = NULL,
    fill = NULL
  ) + 
  scale_fill_discrete(
    labels = c("Not selected", "Selected")
  ) +
  theme_minimal() + 
  theme( axis.text.y = element_blank(),      # drop y-axis numerical values
         axis.ticks.y = element_blank(),     # optional: drop y-axis ticks
         legend.position = "none",          # legend below the graph
         panel.grid = element_blank()
  )  


cov_val_gr5 = cov(
  simRes_DF[simRes_DF$cm == "y ~ treatment", ]$estimate,
  simRes_DF[simRes_DF$cm == "y ~ treatment", ]$selected_pointPrior
)
#cov_lab_gr5 = paste0("Cov(widehat(beta)[m], 1[m]) == ", round(cov_val_gr5, 5))

cov_lab_gr5 = paste0(
  "atop(",
  "Cov(widehat(beta)[m], 1[m]), ",
  "'=' * ", round(cov_val_gr5, 5),
  ")"
)

gr5=ggplot(simRes_DF[simRes_DF$cm == "y ~ treatment",], 
           aes(x = estimate, fill = as.factor(selected_pointPrior))) +
  geom_density(
    data = subset(simRes_DF[simRes_DF$cm == "y ~ treatment",], selected_pointPrior == 1),
    alpha=0.4
  ) +
  geom_density(
    data = subset(simRes_DF[simRes_DF$cm == "y ~ treatment",], selected_pointPrior == 0),
    aes(y = after_stat(density)),
    alpha=0.4
  ) +
  geom_hline(yintercept = 0, color = "black") +
  annotate(
    "text",
    x = -Inf, y = Inf,
    label = cov_lab_gr5,
    hjust = -0.05, vjust = 1.5,
    size = 4,
    parse = TRUE
  )  +
  labs(
    #title = "Selection with final data estimates",
    x = NULL,
    y = NULL,
    fill = NULL
  ) + 
  scale_fill_discrete(
    labels = c("Not selected", "Selected")
  ) +
  theme_minimal() + 
  theme( axis.text.y = element_blank(),      # drop y-axis numerical values
         axis.ticks.y = element_blank(),     # optional: drop y-axis ticks
         legend.position = "none",          # legend below the graph
         panel.grid = element_blank()
  ) 

cov_val_gr6 = cov(
  simRes_DF[simRes_DF$cm == "y ~ treatment + X1 + X2 + X3", ]$estimate,
  simRes_DF[simRes_DF$cm == "y ~ treatment + X1 + X2 + X3", ]$selected_pointPrior
)
#cov_lab_gr6 = paste0("Cov(widehat(beta)[m], 1[m]) == ", round(cov_val_gr6, 5))

cov_lab_gr6 = paste0(
  "atop(",
  "Cov(widehat(beta)[m], 1[m]), ",
  "'=' * ", round(cov_val_gr6, 5),
  ")"
)


gr6=ggplot(simRes_DF[simRes_DF$cm == "y ~ treatment + X1 + X2 + X3",], 
           aes(x = estimate, fill = as.factor(selected_pointPrior))) +
  geom_density(
    data = subset(simRes_DF[simRes_DF$cm == "y ~ treatment + X1 + X2 + X3",], selected_pointPrior == 1),
    alpha=0.4
  ) +
  geom_density(
    data = subset(simRes_DF[simRes_DF$cm == "y ~ treatment + X1 + X2 + X3",], selected_pointPrior == 0),
    aes(y = after_stat(density)),
    alpha=0.4
  ) +
  geom_hline(yintercept = 0, color = "black") +
  annotate(
    "text",
    x = -Inf, y = Inf,
    label = cov_lab_gr6,
    hjust = -0.05, vjust = 1.5,
    size = 4,
    parse = TRUE
  )  +
  labs(
    #title = "Selection with final data estimates",
    x = NULL,
    y = NULL,
    fill = NULL
  ) + 
  scale_fill_discrete(
    labels = c("Not selected", "Selected")
  ) +
  theme_minimal() + 
  theme( axis.text.y = element_blank(),      # drop y-axis numerical values
         axis.ticks.y = element_blank(),     # optional: drop y-axis ticks
         legend.position = "none",          # legend below the graph
         panel.grid = element_blank()
  ) 


gr=grid.arrange(gr1, gr3, gr5, gr2, gr4, gr6, ncol = 3, nrow=2)
gr

ggsave(
  paste0("Output/DistributionSelectionRules_ESMSE=",EffectSizeMSE,".pdf"),
  gr,
  width = 10,
  height = 5,
  dpi = 1200,
  bg='white'
)



