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

#Import functions
source("Package/MSEmin_package v1.0.R")
source("Functions/Functions Simulations.R")
source("Functions/Functions Analyze MSE.R")
source("Functions/Functions Analyze Multivariate Exact Test.R")
source("DGP/DGP_3.R")

#Number of simulations
#S=30000
S=1000

#Compute empirical MSE
set.seed(123)
simRes=pbmclapply(1:S,function(i){
  set.seed(i)
  data_funct=generateDatasetMulti(N_funct=N, EStreatment_funct = beta,
                                  EScontrol_funct = ESvec, CorrMatrix_funct = CorrMat, sigma2_funct=sigma2,
                                  SD_vec_funct = SD_vec)
  cbind(i,do.call(rbind.data.frame, computeStats(data_funct)))
},  mc.cores = min(detectCores()-1,20))
simRes_DF=do.call(rbind.data.frame, simRes)
head(simRes_DF)

#Check estimates
AnalysisTable=analysisTable(simRes_DF)
AnalysisTable

#Show that the statistic follows the appropriate distribution
graphStatMultiExact(simRes_DF)

#False rejection rates
TableFRR=round(rejectionRateMultiExact(simRes_DF),3)
TableFRR

