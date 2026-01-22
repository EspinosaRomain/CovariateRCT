rm(list = ls())
set.seed(123)
options(scipen=999)

if(rstudioapi::isAvailable()) setwd(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))

#Libraries
library(mvtnorm)
library(pbmcapply)
library(goftest)
library(nloptr)

#Import functions
source("Package/MSEmin_package v1.0.R")
source("Functions/Functions Simulations.R")
source("Functions/Functions Analyze MSE.R")

#Vector of effect size multipliers
EffectSizeMSE_vec=seq(1,5,1)
#Results
ResultsTable=matrix(nrow=length(EffectSizeMSE_vec),ncol=4,data=NA)

#Set up parameters
S=10000
B=2000

#Loop over effect sizes
for(i in 1:length(EffectSizeMSE_vec)){
  print(i)
  EffectSizeMSE=EffectSizeMSE_vec[i]
  source("DGP/DGP_2.R")
  
  #Careful:
  #The bootstrap process returns: 0.01 if rejected at 1%, 0.05 if rejected at 5%, 0.10 if rejected at 10%, 1 otherwise
  simRes=pbmclapply(1:S,function(i){
    set.seed(i)
    data_funct=generateDatasetMulti(N_funct=N, EStreatment_funct = beta,
                                    EScontrol_funct = ESvec, CorrMatrix_funct = CorrMat, sigma2_funct=sigma2,
                                    SD_vec_funct = SD_vec)
    pvals_test(fm="y ~ treatment + X1 + X2 + X3", 
               df=data_funct,
               model_j="y ~ treatment + X1",
               model_k="y ~ treatment + X1 + X2",
               tv="treatment",
               method_test=c("LR","LM","Wald","Bootstrap"),
               pval_only=TRUE,
               report_comment=FALSE,
               B_funct=B)
  }, mc.cores = min(detectCores()-1,20))
  simRes_DF=do.call(rbind.data.frame, simRes)
  
  ResultsTable[i,]=c(mean(simRes_DF$pval_Wald<=0.05),
  mean(simRes_DF$pval_LR<=0.05),
  mean(simRes_DF$pval_LM<=0.05),
  mean(simRes_DF$reject_boot<=0.05))
  
  ResultsTable
}
ResultsTable=round(ResultsTable,4)
ResultsTable


