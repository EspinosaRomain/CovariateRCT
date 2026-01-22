rm(list = ls())
set.seed(1234)
options(scipen=999)

if(rstudioapi::isAvailable()) setwd(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))

#library
library(pbmcapply)
library(gridExtra)

#Import functions
source("Package/MSEmin_package v1.0.R")
source("Functions/Functions Simulations.R")
source("Functions/Functions Analyze Univariate Simulations.R")
source("DGP/DGP_1.R")

set.seed(123)
S=30000

simRes=pbmclapply(1:S,function(i){
  set.seed(i)
  data_funct=generateDatasetMulti(N_funct=N, 
                                  EStreatment_funct = beta,
                                  EScontrol_funct = ESvec, 
                                  CorrMatrix_funct = CorrMat, 
                                  sigma2_funct=sigma2,
                                  SD_vec_funct = c(1,1))
  reg_res_funct=summary(lm("y ~ treatment + X1 - 1", data=data_funct))
  list(s=i,gamma=reg_res_funct$coefficients["X1","Estimate"],
       gamma_se=reg_res_funct$coefficients["X1","Std. Error"])
}, mc.cores = detectCores()-1)
simResDF=do.call(rbind.data.frame, simRes)
head(simResDF)

#Graph results
graphStatUni(simResDF)

#False rejection rates
TableFRR=round(rejectionRateUni(simResDF),3)
TableFRR
