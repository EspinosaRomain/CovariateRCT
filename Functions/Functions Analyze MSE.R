#Compute stats
computeStats=function(df_funct){
  controlvar_vec_funct=colnames(df_funct)[-c(1,2,(dim(df_funct)[2]))]
  #Full regression
  formula_string_funct=paste0("y ~ treatment + ",paste0(controlvar_vec_funct,collapse=" + "))
  formula_funct=as.formula(paste0("y ~ treatment + ",paste0(controlvar_vec_funct,collapse=" + ")))
  full_reg_funct=summary(lm(formula_funct, data=df_funct))
  
  all_combinations=unlist(
    sapply(1:length(controlvar_vec_funct), function(k) {
      if(k>0) combn(controlvar_vec_funct, k, simplify = FALSE)
    }),
    recursive = FALSE
  )
  
  res_list=append(list(MSEstatReport(df=df_funct,
                                     fm=formula_string_funct,
                                     cm="y ~ treatment",
                                     tv="treatment")),
                  lapply(all_combinations, function(k) {
                    MSEstatReport(df=df_funct,
                                  fm=formula_string_funct,
                                  cm=paste0("y ~ treatment + ", paste0(unlist(k), collapse=" + ")),
                                  tv="treatment")
                  }))
  names(res_list)=paste0("elem_", seq_along(res_list))
  return(res_list)
}


#Analysis Results Table
analysisTable=function(resSim_funct){
  analTab_funct=matrix(data=NA,nrow=12,ncol=4)
  colnames(analTab_funct)=c("Model_0","Model_1","Model_2","Model_3")
  rownames(analTab_funct)=c("Theoretical bias","Empirical bias","Estimated bias",
                            "Theoretical sqrd bias","(Empirical bias)^2","Estimated sqrd bias corrected",
                            "Theoretical Var","Empirical Var","Estimated Var",
                            "Theoretical MSE", "Empirical MSE","Estimated MSE")
  listSpec_funct=c("y ~ treatment","y ~ treatment + X1","y ~ treatment + X1 + X2", "y ~ treatment + X1 + X2 + X3")
  #Theoretical bias
  analTab_funct[1,]=c((g0%*%ESvec),(g1%*%ESvec),(g2%*%ESvec),(g3%*%ESvec))
  #Empirical Bias
  analTab_funct[2,]=unname(sapply(listSpec_funct, function(i) mean(resSim_funct[resSim_funct$cm==i,"estimate"])))-beta
  #Estimated bias
  analTab_funct[3,]=unname(sapply(listSpec_funct, function(i) mean(resSim_funct[resSim_funct$cm==i,"bias"])))
  #Theoretical Squared Bias
  analTab_funct[4,]=analTab_funct[1,]^2
  #Empirical Squared Bias
  analTab_funct[5,]=analTab_funct[2,]^2
  #Estimated Squared Bias
  analTab_funct[6,]=unname(sapply(listSpec_funct, function(i) mean(resSim_funct[resSim_funct$cm==i,"bias_sqrd_corr"])))
  #Theoretical Variance
  analTab_funct[7,]=c(nu_0,nu_1,nu_2,nu_3)*sigma2
  #Empirical Variance
  analTab_funct[8,]=unname(sapply(listSpec_funct, function(i) var(resSim_funct[resSim_funct$cm==i,"estimate"]))) 
  #Estimated variance
  analTab_funct[9,]=unname(sapply(listSpec_funct, function(i) mean(resSim_funct[resSim_funct$cm==i,"var"])))
  #Theoretical MSE
  analTab_funct[10,]=analTab_funct[4,]+analTab_funct[7,]
  #Empirical MSE
  analTab_funct[11,]=unname(sapply(listSpec_funct, function(i) mean((resSim_funct[resSim_funct$cm==i,"estimate"]-beta)^2))) 
  #Estimated MSE
  analTab_funct[12,]=analTab_funct[6,]+analTab_funct[9,]
  analTab_funct=round(analTab_funct,6)
  analTab_funct
}
