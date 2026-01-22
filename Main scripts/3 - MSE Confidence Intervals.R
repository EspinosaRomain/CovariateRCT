#This script shows that the CI formula are correct: good coverage probability
rm(list = ls())
set.seed(123)
options(scipen=999)

#Working directory
if(rstudioapi::isAvailable()) setwd(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))

#library
library(pbmcapply)
library(gridExtra)

#Import functions
source("Package/MSEmin_package v1.0.R")
source("Functions/Functions Simulations.R")
source("Functions/Functions Analyze MSE.R")
source("Functions/Functions CIs.R")
source("DGP/DGP_2.R")

#Number of simulations
S=30000

#Specify models
models=list("","X1",c("X1","X2"),c("X1","X2","X3"))
list_cm=sapply(models, function(k){
  if(paste0(k, collapse="")==""){
    return("y ~ treatment")
  }else{
    return(paste0("y ~ treatment + ", paste0(k, collapse=" + ")))
  }
})

#ComputeCIs
set.seed(123)
simRes=pbmclapply(1:S,function(i){
  set.seed(i)
  data_funct=generateDatasetMulti(N_funct=N, EStreatment_funct = beta,
                                  EScontrol_funct = ESvec, CorrMatrix_funct = CorrMat, sigma2_funct=sigma2,
                                  SD_vec_funct = SD_vec)
  cbind(i,do.call(rbind.data.frame, lapply(list_cm, function(k){
                                            append(append(list(cm=k),
                                            returnSeveralCIs(df=data_funct, 
                                              fm="y ~ treatment + X1 + X2 + X3",
                                              cm=k,
                                              tv="treatment")),
                                            list(estimate=summary(lm(as.formula(k), data=data_funct))$coefficients[2,1])
                                            )
                                            })
                  )
        )
},  mc.cores = min(detectCores()-1,20))
simResDF=do.call(rbind.data.frame, simRes)
head(simResDF)

#Show the confidence intervals
matResults=analyzeSimCIs(simResDF,cm_list_funct=list_cm)
matResults[,1:2]=round(matResults[,1:2],5)
matResults


