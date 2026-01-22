rm(list = ls())
set.seed(123)
options(scipen=999)

#Working directory
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
source("DGP/DGP_2.R")

#Set up parameters
S=30000
B=2000

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
head(simRes_DF)

#Results 
ResultsTable=matrix(nrow=3,ncol=4,data=NA)
colnames(ResultsTable)=c("Wald","LR","LM","Bootstrap")
rownames(ResultsTable)=c("1%","5%","10%")
ResultsTable[1:3,1]=c(mean(simRes_DF$pval_Wald<=0.01),
                   mean(simRes_DF$pval_Wald<=0.05),
                   mean(simRes_DF$pval_Wald<=0.10))
ResultsTable[1:3,2]=c(mean(simRes_DF$pval_LR<=0.01),
                   mean(simRes_DF$pval_LR<=0.05),
                   mean(simRes_DF$pval_LR<=0.10))
ResultsTable[1:3,3]=c(mean(simRes_DF$pval_LM<=0.01),
                   mean(simRes_DF$pval_LM<=0.05),
                   mean(simRes_DF$pval_LM<=0.10))
ResultsTable[1:3,4]=c(mean(simRes_DF$reject_boot<=0.01),
                    mean(simRes_DF$reject_boot<=0.05),
                    mean(simRes_DF$reject_boot<=0.10))

ResultsTable=round(ResultsTable,4)
ResultsTable
