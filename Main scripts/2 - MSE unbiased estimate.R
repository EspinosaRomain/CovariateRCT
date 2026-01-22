rm(list = ls())
set.seed(123)
options(scipen=999)

if(rstudioapi::isAvailable()) setwd(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))

#library
library(pbmcapply)
library(gridExtra)

#Import functions
source("Package/MSEmin_package v1.0.R")
source("Functions/Functions Simulations.R")
source("Functions/Functions Analyze MSE.R")
source("DGP/DGP_2.R")

#Number of simulations
S=30000

#Launch simulations
simRes=pbmclapply(1:S,function(i){
  set.seed(i)
  data_funct=generateDatasetMulti(N_funct=N, EStreatment_funct = beta,
                                  EScontrol_funct = ESvec, CorrMatrix_funct = CorrMat, sigma2_funct=sigma2,
                                  SD_vec_funct = SD_vec)
  cbind(i,do.call(rbind.data.frame, computeStats(data_funct)))
}, mc.cores = min(detectCores()-1,20))
simRes_DF=do.call(rbind.data.frame, simRes)
simRes_DF

#Results Table
AnalysisTable=analysisTable(simRes_DF)
AnalysisTable