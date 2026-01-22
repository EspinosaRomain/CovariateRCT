#Libraries
library(tidyverse)

#####

diagnosis=function(full_model_formula, 
                   default_model_formula, 
                   treatment_variable,
                   list_candidate_model_formulas, 
                   list_excluded_models_formulas,
                   data,
                   pilot_data,
                   test_method=c("LR"),
                   B_boot_test=1000, 
                   ci_level=0.95){
  #full_model_formula: The formula for the full model (including all relevant covariates)
  #Ex: "y ~ treat + x1 + x2 + x3"
  #default_model_formula: The formula for the default model. If missing, it's the full model.
  #Ex: "y ~ x1 + x2"
  #treatment_variable: name of the treatment variable
  #Ex: "treat"
  #list_candidate_model_formulas: List of formulas for the set of candidate models.
  #If missing: all possible models are investigated
  #Ex: list("y ~ x1 + x2", "y ~ x1 + x3")
  #test_method: Specify the type of tests to compare models
  #Can take values: LR (Likelihood ratio), LM (lagrange multiplier), Wald (Wald test), Bootstrap (Bootstrap test)
  #B_boot_test: Number of draws for the bootstrap test if it is selected. Default=1000 draws.
  #data: dataset on which to estimate the ATE
  #pilotdata: Optional. If provided, it is used to calibrate the test parameters.
  #ci_level: size of the confidence interval
  
  #Retrieve user's call
  call = match.call()
  
  if(missing(treatment_variable)) stop("You must specify treatment_variable")
  tv=call[["treatment_variable"]]
  
  if(missing(full_model_formula)) stop("You must specify full_model_formula")
  fm=call[["full_model_formula"]]
  
  if(!has_linear_treat(fm,tv)) stop("You must include the treatment variable in the full model.")
  fm=reorder_treat_first(fm,tv)
  
  if (!"default_model_formula" %in% names(call)){
    dm=fm
  }else{
    dm=call[["default_model_formula"]]
    if(!has_linear_treat(dm,tv)) stop("You must include the treatment variable in the default model.")
  }
  dm=reorder_treat_first(dm,tv)
  if (!"list_candidate_model_formulas" %in% names(call)){
    cm=all_formulas_with_treat(fm, tv)
  }else{
    cm=as.list(call[["list_candidate_model_formulas"]])[-1]
    cm=reorder_treat_first(cm,tv)
    if(!has_linear_treat(cm,tv)) stop("You must include the treatment variable in all specified candidate models.")
  }
  
  if(!extra_variables(fm,cm)) stop("There are some explanatory variables included in the candidate models that are not included in the full model.")
  
  cat("The full model is: ", fm, "\n")
  cat("The default model is: ", dm, "\n")
  cat("The candidate models are: \n", paste(unlist(cm), collapse = "\n "), "\n")
  
  if(missing(pilot_data)){
    cat("No pilot data provided. \n")
  }else{
    cat("Pilot data provided. \n")
    pd=call[["pilot_data"]]
  }
  
  #Estimated MSEs
  allModels=order_models(unique(unlist(list(fm,dm,cm))), tv)
  nbModels=length(allModels)
  cat("There are ", nbModels, "unique models. \n")
  
  estimatedMSEs=data.frame(model=allModels) 
  estimatedMSEs$estimMSE=NA_real_
  estimatedMSEs$estimSqrdBias=NA_real_
  estimatedMSEs$estimVar=NA_real_
  estimatedMSEs$estimMSEstar=NA_real_
  estimatedMSEs[, paste0("CI_",ci_level*100,"%_lowerBound_MSEstar")]=NA_real_
  estimatedMSEs[, paste0("CI_",ci_level*100,"%_upperBound_MSEstar")]=NA_real_
  if(full_model_formula==default_model_formula) estimatedMSEs$ExactTest_stat=NA_real_
  rownames(estimatedMSEs)=1:nbModels
  #print(estimatedMSEs)
  for(i in 1:nbModels){
    #print(i)
    tmp_res=MSEstatReport(df=data,
                          fm=fm, 
                          cm=estimatedMSEs[i,]$model, 
                          tv=tv)
    estimatedMSEs[i,]$estimMSE=tmp_res$MSE
    estimatedMSEs[i,]$estimSqrdBias=tmp_res$bias_sqrd_corr
    estimatedMSEs[i,]$estimVar=tmp_res$var
    if(full_model_formula==default_model_formula) estimatedMSEs[i,"ExactTest_stat"]=tmp_res$stat
    tmp_res=returnCIs(df=data,
                      fm=fm, 
                      cm=estimatedMSEs[i,]$model, 
                      tv=tv,
                      alpha_funct=1-ci_level)
    estimatedMSEs[i, ]$estimMSEstar=tmp_res$MSEstar
    estimatedMSEs[i, paste0("CI_",ci_level*100,"%_lowerBound_MSEstar")]=tmp_res$lb
    estimatedMSEs[i, paste0("CI_",ci_level*100,"%_upperBound_MSEstar")]=tmp_res$ub
  }
  
  #Model comparisons
  MSEcomparison = as.data.frame(t(combn(allModels, 2)))
  names(MSEcomparison) = c("Model1", "Model2")
  MSEcomparison$estimMSE_Model1=NA_real_
  MSEcomparison$estimMSE_Model2=NA_real_
  #Careful: These are two-sided tests
  if("LR" %in% test_method) MSEcomparison$LR_pval=NA_real_
  if("LM" %in% test_method) MSEcomparison$LM_pval=NA_real_
  if("Wald" %in% test_method) MSEcomparison$Wald_pval=NA_real_
  if("Bootstrap" %in% test_method) MSEcomparison$reject_boot=NA_real_
  if("Bootstrap" %in% test_method) MSEcomparison$Bootstrap_reject=NA_real_
  #This is directional test
  if("Bootstrap" %in% test_method) MSEcomparison$Bootstrap_pval_uni=NA_real_
  #Note: the exact test is one-sided
  #H0: MSE_full<=MSE_reduced
  if(full_model_formula==default_model_formula) MSEcomparison$ExactTest_pval=NA_real_
  #print(MSEcomparison)
  
  #Get pvalues
  for(i in 1:nrow(MSEcomparison)){
    
    #cat("Comparison in progress:", round(i/nrow(MSEcomparison)*100,1), "%")
    MSEcomparison[i,"estimMSE_Model1"]=estimatedMSEs[estimatedMSEs$model==MSEcomparison[i,"Model1"],"estimMSE"]
    MSEcomparison[i,"estimMSE_Model2"]=estimatedMSEs[estimatedMSEs$model==MSEcomparison[i,"Model2"],"estimMSE"]
    tmp_res=pvals_test(fm=fm, 
                       df=data,
                       model_j=MSEcomparison[i,"Model1"],
                       model_k=MSEcomparison[i,"Model2"],
                       tv=tv,
                       method_test=test_method,
                       pval_only=TRUE,
                       B_funct=B_boot_test)
    if("LR" %in% test_method) MSEcomparison[i,"LR_pval"]=tmp_res$pval_LR
    if("LM" %in% test_method) MSEcomparison[i,"LM_pval"]=tmp_res$pval_LM
    if("Wald" %in% test_method) MSEcomparison[i,"Wald_pval"]=tmp_res$pval_Wald
    if("Bootstrap" %in% test_method) MSEcomparison[i,"reject_boot"]=tmp_res$reject_boot
    if("Bootstrap" %in% test_method) MSEcomparison[i,"Bootstrap_reject"]=ifelse(tmp_res$reject_boot<=1-ci_level,"Reject","Fail to reject")
    if("Bootstrap" %in% test_method) MSEcomparison[i,"Bootstrap_pval_uni"]=tmp_res$pval_bootstrap_uni
    if(full_model_formula==default_model_formula & MSEcomparison[i,"Model2"]==fm){
      MSEcomparison[i,"ExactTest_pval"]= pf(estimatedMSEs[estimatedMSEs$model==MSEcomparison[i,"Model1"],"ExactTest_stat"], df1=1, df2=nrow(data)-(sapply(MSEcomparison[i,"Model2"], count_rhs_terms)+1), lower.tail = TRUE, ncp=1)
    }
  }
  
  #Show better models than the default models
  M_tmp=MSEcomparison[MSEcomparison$Model1==dm | MSEcomparison$Model2==dm,]
  alt_model=ifelse(M_tmp$Model1 == dm, M_tmp$Model2, M_tmp$Model1)
  MSE_dm  = ifelse(M_tmp$Model1 == dm, M_tmp$estimMSE_Model1, M_tmp$estimMSE_Model2)
  MSE_alt = ifelse(M_tmp$Model1 == dm, M_tmp$estimMSE_Model2, M_tmp$estimMSE_Model1)
  pval_cols = setdiff(
    names(M_tmp),
    c("Model1", "Model2", "estimMSE_Model1", "estimMSE_Model2")
  )
  MSEcomparison_defaultModel=data.frame(
    Model_dm  = dm,
    Model_alt = alt_model,
    MSE_dm    = MSE_dm,
    MSE_alt   = MSE_alt,
    M_tmp[, pval_cols, drop = FALSE],
    row.names = NULL
  )
  MSEcomparison_defaultModel=MSEcomparison_defaultModel[MSEcomparison_defaultModel$Model_alt %in% cm,]
  
  if("LR" %in% test_method) MSEcomparison_defaultModel$LR_pval_uni=ifelse(MSEcomparison_defaultModel$MSE_alt<MSEcomparison_defaultModel$MSE_dm,0.5*MSEcomparison_defaultModel$LR_pval,1)
  if("LM" %in% test_method) MSEcomparison_defaultModel$LM_pval_uni=ifelse(MSEcomparison_defaultModel$MSE_alt<MSEcomparison_defaultModel$MSE_dm,0.5*MSEcomparison_defaultModel$LM_pval,1)
  if("Wald" %in% test_method) MSEcomparison_defaultModel$Wald_pval_uni=ifelse(MSEcomparison_defaultModel$MSE_alt<MSEcomparison_defaultModel$MSE_dm,0.5*MSEcomparison_defaultModel$Wald_pval,1)

  return(list(estimatedMSEs=estimatedMSEs, 
              MSEcomparison=MSEcomparison, 
              MSEcomparison_defaultModel=MSEcomparison_defaultModel))

}

##########

# Function to compute lower bound of lambda
find_lambda_upper=function(y_funct, df1_funct, df2_funct, alpha_funct) {
  if(pf(y_funct, df1_funct, df2_funct, ncp=0)-alpha_funct/2>0){
    uniroot(function(lambda_funct) pf(y_funct, df1_funct, df2_funct, ncp=lambda_funct) - alpha_funct / 2,
            lower = 0, upper = 1000000, tol=.Machine$double.eps)$root
  }else{
    0
  }
}

# Function to compute upper bound of lambda
find_lambda_lower=function(y_funct, df1_funct, df2_funct, alpha_funct) {
  if(pf(y_funct, df1_funct, df2_funct, ncp=0)-(1-alpha_funct/2)>0){
    uniroot(function(lambda_funct) pf(y_funct, df1_funct, df2_funct, ncp=lambda_funct) - (1 - alpha_funct / 2),
            lower = 0, upper = 100000000, tol=.Machine$double.eps)$root
  }else{
    0
  }
}


########

#Function to return the CIs
returnCIs=function(df, fm, cm, tv, alpha_funct=0.05){
  
  #Prepare data
  y_var=get_lhs(fm)
  x_var=vars_in_formula(fm)
  
  if(has_rhs(gsub(tv, "", cm))){
    x_cm_var=vars_in_formula(gsub(tv, "", cm))
  }else{
    x_cm_var=NULL
  }
  
  #Full model
  x_funct=as.matrix(cbind(1,df[,x_var])) #add intercept
  N_funct=nrow(x_funct)
  K_funct=ncol(x_funct)
  y_funct=as.matrix(df[,y_var]) 
  xtxm1_funct=solve(crossprod(x_funct))
  full_reg_funct=summary(lm(as.formula(fm), data=df))
  
  #Other matrices and vectors
  x_incl_funct=as.matrix(x_funct[,c(1, tv, x_cm_var)]) 
  x_excl_funct=as.matrix(x_funct[ , setdiff(colnames(x_funct), c(1, tv, x_cm_var))])
  g_star_funct=genGstarVec(x_funct,genVecInclusionControls(colnames(x_funct)[-c(1,2)],x_cm_var))
  
  #Get MSEstat
  res_MSE_funct=MSEstatReport(df=df,
                              fm=fm, 
                              cm=cm, 
                              tv=tv)
  
  #Lower bound estimate
  lb = ifelse(fm==cm,NA,(find_lambda_lower(res_MSE_funct$stat, df1_funct = 1, df2_funct = N_funct - K_funct, alpha_funct = alpha_funct) - 1) *  as.numeric(res_MSE_funct$correction_unscaled) + res_MSE_funct$var / full_reg_funct$sigma^2)
  
  #Upper bound estimate
  ub = ifelse(fm==cm,NA,(find_lambda_upper(res_MSE_funct$stat, df1_funct = 1, df2_funct = N_funct - K_funct, alpha_funct = alpha_funct) - 1) *  as.numeric(res_MSE_funct$correction_unscaled) + res_MSE_funct$var / full_reg_funct$sigma^2)
  
  #Return results
  return(list(MSEstar= res_MSE_funct$stat * as.numeric(res_MSE_funct$correction_unscaled) + res_MSE_funct$var / full_reg_funct$sigma^2, lb=lb, ub=ub))
  
  return(out_funct)
}


#Function to return the CIs
returnSeveralCIs=function(df, fm, cm, tv, alpha_funct=c(0.10,0.05,0.01)){
  unlist(lapply(1:length(alpha_funct), function(k){
    res_tmp=returnCIs(df, fm, cm, tv, alpha_funct=alpha_funct[k])
    if(k==1) names(res_tmp)=c("MSEstar",paste0("lb_",(1-alpha_funct[k])*100),paste0("ub_",(1-alpha_funct[k])*100))
    if(k>1){
      res_tmp=res_tmp[-1]
      names(res_tmp)=c(paste0("lb_",(1-alpha_funct[k])*100),paste0("ub_",(1-alpha_funct[k])*100))
    }
    res_tmp
  }))
}

########

#Define log-likelihood function for the unrestricted model
loglik_unrestricted=function(par_funct, X_funct, y_funct) {
  length_funct=length(par_funct)
  beta_funct = par_funct[1:(length_funct-1)]
  sigma2_funct = exp(par_funct[length_funct])
  res_funct = y_funct - X_funct %*% beta_funct
  n_funct = length(y_funct)
  -0.5 * n_funct * log(2 * pi * sigma2_funct) - 0.5 * sum(res_funct^2) / sigma2_funct
}


########

#Function to check if an argument is an integer
is.wholenumber =
  function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

#Define log-likelihood function for the unrestricted model
loglik_restricted = function(par_funct, x_funct, y_funct, model_j_controls_funct, model_k_controls_funct, gj_star_funct, gk_star_funct, nu_j_funct, nu_k_funct) {
  
  if(all(genVecInclusionControls(colnames(x_funct)[-c(1,2)],model_j_controls_funct)==genVecInclusionControls(colnames(x_funct)[-c(1,2)],model_k_controls_funct))) stop("The two models are identical")
  
  if(missing(gj_star_funct)){
    gj_star_funct=genGstarVec(x_funct,genVecInclusionControls(colnames(x_funct)[-c(1,2)],model_j_controls_funct))  
  }
  if(missing(gk_star_funct)){
    gk_star_funct=genGstarVec(x_funct,genVecInclusionControls(colnames(x_funct)[-c(1,2)],model_k_controls_funct))
  }
  
  if(missing(nu_j_funct)){
    nu_j_funct=genNuVec(x_funct,genVecInclusionControls(colnames(x_funct)[-c(1,2)],model_j_controls_funct))
  }
  
  if(missing(nu_k_funct)){
    nu_k_funct=genNuVec(x_funct,genVecInclusionControls(colnames(x_funct)[-c(1,2)],model_k_controls_funct))
  }
  
  sigma2_funct=(as.double((gj_star_funct%*%par_funct)^2)-as.double((gk_star_funct%*%par_funct)^2))/(nu_k_funct-nu_j_funct)
  
  #if (sigma2_funct <= 0) return(-Inf)
  if (sigma2_funct <= 0) {
    return(-1e10 - 1e8 * abs(sigma2_funct))  # large, finite, decreasing
  }
  
  res_funct = y_funct - x_funct %*% par_funct
  n_funct = length(y_funct)
  -0.5 * n_funct * log(2 * pi * sigma2_funct) - 0.5 * sum(res_funct^2) / sigma2_funct
}

##########


MSEstatReport=function(df, fm, cm, tv){
  
  fm=unname(fm)
  
  #Prepare data
  y_var=get_lhs(fm)
  x_var=vars_in_formula(fm)
  
  if(has_rhs(gsub(tv, "", cm))){
    x_cm_var=vars_in_formula(gsub(tv, "", cm))
  }else{
    x_cm_var=NULL
  }
  
  #Full model
  x_funct=as.matrix(cbind(1,df[,x_var])) #add intercept
  y_funct=as.matrix(df[,y_var]) 
  xtxm1_funct=solve(crossprod(x_funct))
  full_reg_funct=summary(lm(as.formula(fm), data=df))
  
  #Other matrices and vectors
  x_incl_funct=as.matrix(x_funct[,c(1, tv, x_cm_var)]) 
  x_excl_funct=as.matrix(x_funct[ , setdiff(colnames(x_funct), c(1, tv, x_cm_var))])
  g_star_funct=genGstarVec(x_funct,genVecInclusionControls(colnames(x_funct)[-c(1,2)],x_cm_var))
  
  #Estimate
  estimate_funct=unname((solve(crossprod(x_incl_funct))%*%crossprod(x_incl_funct,y_funct))[2,])
  #Bias
  bias_funct=as.numeric(stats::coef(full_reg_funct)[,1]%*%g_star_funct)
  #Correction factor unscaled
  correction_funct=g_star_funct%*%xtxm1_funct%*%g_star_funct
  #Squared bias corrected
  bias_sqrd_corrected_funct=as.numeric((bias_funct)^2-correction_funct*full_reg_funct$sigma^2)
  #Variance
  nu_funct=genNuVec(x_funct,genVecInclusionControls(colnames(x_funct)[-c(1,2)],x_cm_var))
  #var_funct=genNuVec(x_funct,genVecInclusionControls(colnames(x_funct)[-c(1,2)],x_cm_var))*full_reg_funct$sigma^2
  var_funct=nu_funct*full_reg_funct$sigma^2
  #Stat
  stat_funct=as.numeric((bias_funct/full_reg_funct$sigma)^2/(g_star_funct%*%xtxm1_funct%*%g_star_funct))
  
  
  #Return
  list(cm=cm,
       estimate=estimate_funct, 
       bias=bias_funct,
       bias_sqrd_corr=bias_sqrd_corrected_funct,
       var=var_funct, 
       MSE=bias_sqrd_corrected_funct+var_funct, 
       correction_unscaled=as.numeric(correction_funct), 
       stat=stat_funct,
       nu=nu_funct,
       full_sigma2=full_reg_funct$sigma^2)
}

########

genVecInclusionControls=function(controlvar_vec_funct,controls_funct){
  if(is.null(controls_funct)){
    vec_logic_inclusion_all=c(TRUE,TRUE,rep(FALSE,length(controlvar_vec_funct)))  
  }else if(length(controls_funct)==1){
    if(controls_funct=="" | controls_funct=="None"){
      vec_logic_inclusion_all=c(TRUE,TRUE,rep(FALSE,length(controlvar_vec_funct)))  
    }else{
      if(!all(controls_funct %in% controlvar_vec_funct)) stop("Some listed control variables are not in the dataset.")
      vec_logic_inclusion_all=c(TRUE,TRUE,controlvar_vec_funct %in% controls_funct)
    }
  }else{
    if(!all(controls_funct %in% controlvar_vec_funct)) stop("Some listed control variables are not in the dataset.")
    vec_logic_inclusion_all=c(TRUE,TRUE,controlvar_vec_funct %in% controls_funct)
  }
  return(vec_logic_inclusion_all)
}


########

genGstarVec=function(x_funct,vec_logic_inclusion_all){
  x_incl_funct=as.matrix(x_funct[,vec_logic_inclusion_all]) 
  x_excl_funct=as.matrix(x_funct[,!vec_logic_inclusion_all])
  
  #(S'S)^{-1}S'Z[1,]
  g_funct=(solve(crossprod(x_incl_funct))%*%crossprod(x_incl_funct,x_excl_funct))[2,]
  g_star_funct=numeric(length(vec_logic_inclusion_all))
  g_star_funct[vec_logic_inclusion_all]=0              
  g_star_funct[!vec_logic_inclusion_all]=g_funct       
  return(g_star_funct)
}

genNuVec=function(x_funct,vec_logic_inclusion_all){
  x_incl_funct=as.matrix(x_funct[,vec_logic_inclusion_all]) 
  solve(crossprod(x_incl_funct))[2,2]
}

##########

pvals_test=function(fm,model_j,model_k,df,tv, method_test=c("LR","LM","Wald","Bootstrap"), methodMLE_funct="Nelder-Mead", pval_only=TRUE, report_comment=FALSE, B_funct=10){
  
  comment_funct=0
  
  #Prepare data
  y_var=get_lhs(fm)
  x_var=vars_in_formula(fm)
  
  if(has_rhs(gsub(tv, "", model_j))){
    model_j_controls_funct=vars_in_formula(gsub(tv, "", model_j))
  }else{
    model_j_controls_funct=NULL
  }
  if(has_rhs(gsub(tv, "", model_k))){
    model_k_controls_funct=vars_in_formula(gsub(tv, "", model_k))
  }else{
    model_k_controls_funct=NULL
  }
  
  #Matrix data
  x_funct=as.matrix(cbind(1,df[,x_var])) #add intercept
  y_funct=as.matrix(df[,y_var]) 
  nb_explanatoryVar_funct=dim(x_funct)[2]
  N_funct=nrow(df)
  
  #Get parameters of interest
  gj_star_funct=genGstarVec(x_funct,genVecInclusionControls(colnames(x_funct)[-c(1,2)],model_j_controls_funct))
  gk_star_funct=genGstarVec(x_funct,genVecInclusionControls(colnames(x_funct)[-c(1,2)],model_k_controls_funct))
  nu_j_funct=genNuVec(x_funct,genVecInclusionControls(colnames(x_funct)[-c(1,2)],model_j_controls_funct))
  nu_k_funct=genNuVec(x_funct,genVecInclusionControls(colnames(x_funct)[-c(1,2)],model_k_controls_funct))
  gamma_var_funct=!(genVecInclusionControls(colnames(x_funct)[-c(1,2)],model_j_controls_funct)==TRUE & genVecInclusionControls(colnames(x_funct)[-c(1,2)],model_k_controls_funct)==TRUE)
  
  #Get OLS estimates
  XtXm1_funct=solve(crossprod(x_funct))
  beta_OLS=XtXm1_funct %*% crossprod(x_funct, y_funct)
  res_funct=y_funct - x_funct %*% beta_OLS
  SSR_funct=as.numeric(crossprod(res_funct))
  sigma2_hat=SSR_funct/(N_funct - nb_explanatoryVar_funct)
  Var_sigma2_funct=2*sigma2_hat^2/(N_funct - nb_explanatoryVar_funct)
  Var_beta_funct=sigma2_hat*XtXm1_funct
  
  if(abs(nu_k_funct-nu_j_funct)<0.00000000001){
    comment_funct=1
    #Both have the same variance
    if(pval_only==FALSE) res=list(fm=fm, model_j=model_j, model_k=model_k,df=df,tv=tv)
    if(pval_only==TRUE) res=list()
    if("LR" %in% method_test) res=append(res,list(pval_LR=NA_real_))
    if("LM" %in% method_test) res=append(res,list(pval_LM=NA_real_))
    if("Wald" %in% method_test) res=append(res,list(pval_Wald=NA_real_))
    if(report_comment) res=append(res,list(comment=comment_funct))
    return(res)
  }
  
  bias_j_funct=as.numeric(t(gj_star_funct) %*% beta_OLS)
  bias_k_funct=as.numeric(t(gk_star_funct) %*% beta_OLS)
  sigma2_constr=(bias_j_funct^2-bias_k_funct^2)/(nu_k_funct-nu_j_funct)
  
  #For Wald test
  if("Wald" %in% method_test){
    #First order derivatives - Jacobian
    h_funct=bias_j_funct^2 - bias_k_funct^2 - (nu_k_funct-nu_j_funct) * sigma2_hat
    Gbeta_funct=2*bias_j_funct*gj_star_funct - 2*bias_k_funct*gk_star_funct
    est_var_funct=as.numeric(t(Gbeta_funct) %*% Var_beta_funct %*% Gbeta_funct) + (nu_k_funct-nu_j_funct)^2 * Var_sigma2_funct
    W_funct=as.numeric(h_funct^2 / est_var_funct)
    pvalWald_funct=1-pchisq(W_funct, df = 1)
  }
  
  #Unrestricted likelihood estimates: for LR test only
  if("LR" %in% method_test){
    #Unrestricted model
    init_unres_funct=c(beta_OLS, sigma2_hat)
    fit_unres_funct=optim(init_unres_funct, fn = function(par) -loglik_unrestricted(par, x_funct, y_funct),
                          method = methodMLE_funct, control=list(maxit=10000, reltol = 1e-12))
    ll_unres_funct=-fit_unres_funct$value
  }
  
  #Restricted likelihood estimates: for LR and LM tests
  if(("LR" %in% method_test) | ("LM" %in% method_test)){
    #Restricted model
    resModelFail=FALSE
    
    #Options
    opts_list <- list(
      "algorithm" = "NLOPT_LN_COBYLA",  # Derivative-free, handles constraints
      "xtol_rel" = 1.0e-8,
      "maxeval" = 10000
    )
    
    
    result_nlopt = nloptr(
      x0 = beta_OLS,
      eval_f = neg_loglik_restricted,
      lb = rep(-Inf,length(beta_OLS)), ub = rep(Inf,length(beta_OLS)),  # No box constraints on parameters
      eval_g_ineq = con_func,
      x_funct = x_funct,
      y_funct = y_funct,
      model_j_controls_funct=model_j_controls_funct, 
      model_k_controls_funct=model_k_controls_funct,
      gj_star_funct = gj_star_funct,
      gk_star_funct = gk_star_funct,
      nu_j_funct = nu_j_funct,
      nu_k_funct = nu_k_funct,
      opts = opts_list
    )
    
    ll_res_funct=-result_nlopt$objective
  }
  
  #LR test
  if(("LR" %in% method_test)){
    if(resModelFail){
      LR_stat_funct=NA
      pval_LR_funct=NA
    }else{
      LR_stat_funct=2 * (ll_unres_funct-ll_res_funct)
      pval_LR_funct=1 - pchisq(LR_stat_funct, df = 1)    
    }
  }
  
  #LM test
  if(("LM" %in% method_test)){
    if(resModelFail){
      LM_funct=NA
      pval_LM_funct=NA
    }else{
      #Get estimates from the restricted model
      beta_restrictLL=result_nlopt$solution
      num_funct=(as.double((gj_star_funct%*%beta_restrictLL)^2)-as.double((gk_star_funct%*%beta_restrictLL)^2))
      denom_funct=(nu_k_funct-nu_j_funct)
      sigma2_LLrestrict_funct=num_funct/denom_funct
      
      #Generate first derivatives
      residualsLLrestrict_funct=y_funct - as.vector(x_funct %*% beta_restrictLL) #Get the residuals
      s_beta_funct=as.vector( crossprod(x_funct, residualsLLrestrict_funct) / sigma2_LLrestrict_funct ) #First derivative of LL with coefs # X'(y - Xb)/sigma^2
      s_sig2_funct=-N_funct/(2*sigma2_LLrestrict_funct) + sum(residualsLLrestrict_funct^2)/(2*sigma2_LLrestrict_funct^2)  #First derivative of LL with sigma2
      s_funct=c(s_beta_funct, s_sig2_funct)
      
      #Fisher Info Matrix
      I_bb=crossprod(x_funct) / sigma2_LLrestrict_funct #X'X / σ^2
      I_bs2=matrix(0, nrow = nb_explanatoryVar_funct, ncol = 1) 
      I_s2s2=matrix(N_funct/(2*sigma2_LLrestrict_funct^2), nrow = 1, ncol = 1) # N/(2σ^4)
      I_funct=rbind(cbind(I_bb, I_bs2),cbind(t(I_bs2), I_s2s2))
      
      #Inverse of Fisher Info Mat
      I_inv=rbind(
        cbind( sigma2_LLrestrict_funct * solve(crossprod(x_funct)), matrix(0, nb_explanatoryVar_funct, 1) ),          # (X'X)^{-1} * σ^2
        cbind( matrix(0, 1, nb_explanatoryVar_funct), 2*sigma2_LLrestrict_funct^2/N_funct )                           # 2σ^4/n
      )
      #Jacobian of r(θ) = (a'β)^2 - (b'β)^2 - σ^2*c
      aTb=as.numeric(crossprod(gj_star_funct, beta_restrictLL))
      bTb=as.numeric(crossprod(gk_star_funct, beta_restrictLL))
      dr_dbeta=2*aTb*gj_star_funct - 2*bTb*gk_star_funct
      dr_dsig2=-(nu_k_funct-nu_j_funct)
      R_funct=c(dr_dbeta, dr_dsig2)   # (p+1)-vector
      # LM statistic: (R' I^{-1} s)' [R' I^{-1} R]^{-1} (R' I^{-1} s) ~ Chi^2_1
      Iv_s=I_inv %*% s_funct
      Iv_R=I_inv %*% R_funct
      #Compute the LM stat & p-value
      LM_funct=as.numeric( (t(R_funct) %*% Iv_s)^2 / (t(R_funct) %*% Iv_R) )
      pval_LM_funct=1-pchisq(LM_funct, df = 1)
    }}
  
  if("Bootstrap" %in% method_test){
    reg_full_funct=lm(as.formula(fm), data=df)
    fittedval_funct=fitted(reg_full_funct)
    resid_funct=resid(reg_full_funct)
    
    #Get bootstrapped stat
    T_boot=t(sapply(1:B_funct, function(i){
      df_funct_boot=df
      df_funct_boot$y=simulateData(fittedval_funct,resid_funct)
      MSE_j_boot=MSEstatReport(df=df_funct_boot, fm=fm, cm=model_j, tv=tv)$MSE
      MSE_k_boot=MSEstatReport(df=df_funct_boot, fm=fm, cm=model_k, tv=tv)$MSE
      return(unname(MSE_j_boot-MSE_k_boot))
    }))
    T_boot=as.vector(T_boot)
    
    #From Horowitz
    T_bar=mean(T_boot)
    reject_boot_funct=ifelse(0<=quantile(T_boot, c(0.05)) | 0>=quantile(T_boot, c(0.95)), 0.10, 1)
    reject_boot_funct=ifelse(0<=quantile(T_boot, c(0.025)) | 0>=quantile(T_boot, c(0.975)), 0.05, reject_boot_funct)
    reject_boot_funct=ifelse(0<=quantile(T_boot, c(0.005)) | 0>=quantile(T_boot, c(0.995)), 0.01, reject_boot_funct)
    percent_funct=mean(0<=T_boot)
    
    # reject_boot_uni_funct=ifelse(0<=quantile(T_boot, c(0.10)), 0.10, 1)
    # reject_boot_uni_funct=ifelse(0<=quantile(T_boot, c(0.05)), 0.05, reject_boot_funct)
    # reject_boot_uni_funct=ifelse(0<=quantile(T_boot, c(0.01)), 0.01, reject_boot_funct)
    pval_boot_uni_funct=ecdf(T_boot)(0)
  }
  
  if(pval_only==FALSE) res=list(fm=fm, model_j=model_j, model_k=model_k,df=df,tv=tv)
  if(pval_only==TRUE) res=list()
  if("LR" %in% method_test) res=append(res,list(pval_LR=pval_LR_funct))
  if("LM" %in% method_test) res=append(res,list(pval_LM=pval_LM_funct))
  if("Wald" %in% method_test) res=append(res,list(pval_Wald=pvalWald_funct))
  if("Bootstrap" %in% method_test) res=append(res,list(reject_boot=reject_boot_funct, pval_bootstrap_uni=pval_boot_uni_funct))
  if(report_comment) res=append(res,list(comment=comment_funct))
  
  return(res)
}

########

#Compute all possible model combinations that include treat
all_formulas_with_treat = function(f, treat) {
  #f: character string, for example "y ~ x1 + x2 + x3"
  #treat: character string, for example "treat"
  
  #split formula into LHS and RHS
  parts = strsplit(f, "~", fixed = TRUE)[[1]]
  if (length(parts) != 2L) stop("Invalid formula string.")
  
  response = trimws(parts[1])
  rhs      = trimws(parts[2])
  
  #get RHS variables
  predictors = trimws(strsplit(rhs, "\\+")[[1]])
  predictors = predictors[predictors != ""]   # drop empties if any
  
  #separate 'treat' from other predictors
  others = setdiff(predictors, treat)
  
  #all subsets of "others"
  if (length(others) == 0L) {
    subsets = list(character(0))   # only 'treat' will be used
  } else {
    subsets = unlist(
      lapply(0:length(others), function(k)
        combn(others, k, simplify = FALSE)
      ),
      recursive = FALSE
    )
  }
  
  #build formulas: response ~ treat [+ subset of others]
  formula_strings = vapply(subsets, function(vars) {
    rhs_terms = c(treat, vars)
    paste(response, "~", paste(rhs_terms, collapse = " + "))
  }, FUN.VALUE = character(1L))
  
  unique(formula_strings)
}

#Check that a specification includes linearly treat as a regressor
has_linear_treat = function(formulas, treat) {
  # Accept a single string, a character vector, or a list
  formulas = unlist(formulas)
  
  res=sapply(formulas, function(fs) {
    f = as.formula(fs)
    tm = attr(terms(f), "term.labels")  # model terms on RHS
    treat %in% tm                        # TRUE if treat is a main effect term
  })
  
  return(all(res))
  
}


########

#Reorder independent variable names
reorder_treat_first = function(formulas, treatment_variable) {
  
  reorder_one = function(fs) {
    # Convert to formula object
    f = as.formula(fs)
    
    # Extract LHS and RHS
    lhs = trimws(deparse(f[[2]]))
    rhs = trimws(deparse(f[[3]]))
    
    # Extract variables on RHS
    vars = trimws(strsplit(rhs, "\\+")[[1]])
    
    # Remove empty strings
    vars = vars[vars != ""]
    
    # Reorder: treat first, then all others
    new_vars = c(
      treatment_variable,
      setdiff(vars, treatment_variable)
    )
    
    # Build reformatted formula string
    paste(lhs, "~", paste(new_vars, collapse = " + "))
  }
  
  # Apply to list or vector
  vapply(formulas, reorder_one, FUN.VALUE = character(1))
}


########

#Get all variables in a formula string
vars_in_formula = function(f) {
  f = as.formula(f)
  tm = attr(terms(f), "term.labels")   # RHS terms
  all.vars(parse(text = tm))            # break interactions, etc.
}


########

#Function that returns TRUE if all variables listed in cm are also included in fm
#FALSE otherwise
extra_variables = function(fm, cm) {
  # Variables in the full model
  fm_vars = vars_in_formula(fm)
  
  # Variables in each candidate model (union of all)
  cm_vars = unique(unlist(lapply(cm, vars_in_formula)))
  
  # Variables present in cm but not fm
  ifelse(length(setdiff(cm_vars, fm_vars))>0,FALSE,TRUE)
}

########

#Function to order the models
order_models = function(formulas, treatment_variable = "treat") {
  
  # Helper to standardize one formula: treat first, controls sorted
  normalize_formula = function(fs) {
    f   = as.formula(fs)
    lhs = trimws(deparse(f[[2]]))
    rhs = trimws(deparse(f[[3]]))
    
    vars = trimws(strsplit(rhs, "\\+")[[1]])
    vars = vars[vars != ""]
    
    controls = setdiff(vars, treatment_variable)
    controls_sorted = sort(controls)
    
    paste(lhs, "~", paste(c(treatment_variable, controls_sorted), collapse = " + "))
  }
  
  # Normalize all formulas
  norm = vapply(formulas, normalize_formula, FUN.VALUE = character(1))
  
  # Count controls in each normalized formula
  n_controls = vapply(norm, function(fs) {
    f   = as.formula(fs)
    rhs = trimws(deparse(f[[3]]))
    vars = trimws(strsplit(rhs, "\\+")[[1]])
    sum(vars != treatment_variable)
  }, FUN.VALUE = integer(1))
  
  # Order: first by number of controls, then alphabetically
  o = order(n_controls, norm)
  
  norm[o]
}


########

get_lhs = function(model_string) {
  f = as.formula(model_string)
  deparse(f[[2]])
}

########

count_rhs_terms = function(f) {
  length(attr(terms(as.formula(f)), "term.labels"))
}





########

#Check if equation has a right-hand side side
has_rhs = function(fstring) {
  parts = strsplit(fstring, "~", fixed = TRUE)[[1]]
  if (length(parts) != 2) return(FALSE)
  rhs = trimws(parts[2])
  rhs != ""
}


########

#Try to load library
load_or_install = function(pkg) {
  # Try to load the package
  if (!requireNamespace(pkg, quietly = TRUE)) {
    
    message(sprintf("Package '%s' is not installed.", pkg))
    ans <- readline(prompt = "Do you want to install it? [y/n]: ")
    
    if (tolower(ans) %in% c("y", "yes")) {
      message(sprintf("Installing package '%s'...", pkg))
      install.packages(pkg)
      
      # Try loading again after installation
      if (!requireNamespace(pkg, quietly = TRUE)) {
        stop(sprintf("Installation of '%s' failed. Cannot continue.", pkg),
             call. = FALSE)
      }
      
      library(pkg, character.only = TRUE)
      message(sprintf("Package '%s' successfully installed and loaded.", pkg))
      
    } else {
      stop(sprintf("Package '%s' is required. Execution stopped.", pkg),
           call. = FALSE)
    }
    
  } else {
    # Package was already installed
    library(pkg, character.only = TRUE)
  }
}

########

#Function for non-linear optimization
neg_loglik_restricted = function(par_funct, x_funct, y_funct, model_j_controls_funct, model_k_controls_funct, gj_star_funct, gk_star_funct, nu_j_funct, nu_k_funct) {
  -loglik_restricted(par_funct, x_funct, y_funct, model_j_controls_funct, model_k_controls_funct, gj_star_funct, gk_star_funct, nu_j_funct, nu_k_funct)
}

########

# Define the constraint: sigma2 > 0
# For nloptr, we need to express this as a nonlinear inequality constraint:
# sigma2(par) > 0 => sigma2(par) >= 1e-6 (for numerical stability)
# Note : for nloptr, the constraint must be expressed in the negative domain
sigma2_constraint_function = function(par_funct, x_funct, y_funct, model_j_controls_funct, model_k_controls_funct, gj_star_funct, gk_star_funct, nu_j_funct, nu_k_funct) {
  sigma2_funct = as.double(((gj_star_funct %*% par_funct)^2 - (gk_star_funct %*% par_funct)^2) / (nu_k_funct - nu_j_funct))
  1e-6 - sigma2_funct  
}

########

# Wrap the constraint for nloptr
con_func = function(x, x_funct, y_funct, model_j_controls_funct, model_k_controls_funct, gj_star_funct, gk_star_funct, nu_j_funct, nu_k_funct) {
  c(sigma2_constraint_function(x, x_funct, y_funct, model_j_controls_funct, model_k_controls_funct, gj_star_funct, gk_star_funct, nu_j_funct, nu_k_funct))
}

########

#Simulate data with wild bootstrap
simulateData=function(fittedval_funct,resid_funct){
  #Generate Rademacher weights for wild bootstrap
  weights_funct=sample(c(-1,1), size=length(fittedval_funct), replace=TRUE)
  #Residual boostrap
  return(fittedval_funct+weights_funct*resid_funct)
}