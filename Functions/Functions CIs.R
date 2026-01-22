analyzeSimCIs=function(simCI_DF_funct,cm_list_funct){
  matResults_funct=matrix(data=NA,ncol=5,nrow=4)
  
  colnames(matResults_funct)=c("TrueMSE_star","EmpiricalMSE_star","covProba90","covProba95","covProba99")
  rownames(matResults_funct)=c("No control","1 covariate","2 covariates","All covariates")
  
  matResults_funct[,1]=c(theoreticalMSE_0,theoreticalMSE_1,theoreticalMSE_2,theoreticalMSE_3)/sigma2
  
  matResults_funct[,2]=sapply(cm_list_funct, function(k_funct) mean((simCI_DF_funct[simCI_DF_funct$cm==k_funct,]$estimate-beta)^2)/sigma2)
  
  matResults_funct[,3]=sapply(1:length(cm_list_funct), function(k_funct) 1-(round(mean(matResults_funct[k_funct,1]<=simCI_DF_funct[simCI_DF_funct$cm==cm_list_funct[k_funct],]$lb_90)+mean(simCI_DF_funct[simCI_DF_funct$cm==cm_list_funct[k_funct],]$ub_90<=matResults_funct[k_funct,1]),3)))
  
  matResults_funct[,4]=sapply(1:length(cm_list_funct), function(k_funct) 1-(round(mean(matResults_funct[k_funct,1]<=simCI_DF_funct[simCI_DF_funct$cm==cm_list_funct[k_funct],]$lb_95)+mean(simCI_DF_funct[simCI_DF_funct$cm==cm_list_funct[k_funct] ,]$ub_95<=matResults_funct[k_funct,1]),3)))
  
  matResults_funct[,5]=sapply(1:length(cm_list_funct), function(k_funct) 1-(round(mean(matResults_funct[k_funct,1]<=simCI_DF_funct[simCI_DF_funct$cm==cm_list_funct[k_funct],]$lb_99)+mean(simCI_DF_funct[simCI_DF_funct$cm==cm_list_funct[k_funct],]$ub_99<=matResults_funct[k_funct,1]),3)))
  
  matResults_funct
}


