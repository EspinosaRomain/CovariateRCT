library(mvtnorm)

computeStatsSelectionRules=function(main_data_funct,
                                    selectionRule=c("Final Sample"), #selectionRule=c("Final Sample","Pilot","Point Prior")
                                    pilot_df,
                                    priors_mu=NULL, 
                                    priors_var=NULL,
                                    alpha_funct=0.10,
                                    R_funct=100){
  
  controlvar_vec_funct=colnames(main_data_funct)[-c(1,2,(dim(main_data_funct)[2]))]

  #Full model  
  formula_string_funct=paste0("y ~ treatment + ",paste0(controlvar_vec_funct,collapse=" + "))
  formula_funct=as.formula(paste0("y ~ treatment + ",paste0(controlvar_vec_funct,collapse=" + ")))
  
  #Full model with main data
  full_reg_funct=summary(lm(formula_funct, data=main_data_funct))
  
  #Full model with main data
  pilot_full_reg_funct=summary(lm(formula_funct, data=pilot_df))
  
  all_combinations=unlist(
    sapply(1:length(controlvar_vec_funct), function(k) {
      if(k>0) combn(controlvar_vec_funct, k, simplify = FALSE)
    }),
    recursive = FALSE
  )
  
  res_list=append(list(MSEstatReportSelectionRule(df=main_data_funct,
                                                  fm=formula_string_funct,
                                                  cm="y ~ treatment",
                                                  tv="treatment", 
                                                  selectionRule=selectionRule, 
                                                  pilot_df=pilot_df, 
                                                  priors_mu=priors_mu, 
                                                  priors_var=priors_var,
                                                  alpha_funct=alpha_funct,
                                                  R_funct=R_funct)),
                  lapply(all_combinations, function(k) {
                    MSEstatReportSelectionRule(df=main_data_funct,
                                               fm=formula_string_funct,
                                               cm=paste0("y ~ treatment + ", paste0(unlist(k), collapse=" + ")),
                                               tv="treatment", 
                                               selectionRule=selectionRule, 
                                               pilot_df=pilot_df, 
                                               priors_mu=priors_mu, 
                                               priors_var=priors_var,
                                               alpha_funct=alpha_funct,
                                               R_funct=R_funct)
                  }))
  print(res_list)
  names(res_list)=paste0("elem_", seq_along(res_list))
  return(res_list)
}


MSEstatReportSelectionRule=function(df, 
                                          fm, 
                                          cm, 
                                          tv,
                                          selectionRule = c("Final Sample"), #c("Final Sample", "Pilot", "Point Prior")
                                          pilot_df=NULL, 
                                          priors_mu=NULL, #Vector: length=number of covariates. Put 0 if variable included in candidate model
                                          priors_var=NULL,  #Vector 
                                          R_funct=100, 
                                          alpha_funct=0.10){
  
  fm=unname(fm)
  print(cm)
  
  #Generate dummies
  finalSample_based_selection=ifelse("Final Sample" %in% selectionRule & cm!=fm,1,0)
  pilot_based_selection=ifelse("Pilot" %in% selectionRule & cm!=fm,1,0)
  priorPoint_based_selection=ifelse("Point Prior" %in% selectionRule & cm!=fm,1,0)
  
  if(missing(pilot_df) & pilot_based_selection==1) stop("Pilot Data missing")
  if(missing(priors_mu) & priorPoint_based_selection==1) stop("Point priors missing")
  
  
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
  var_funct=nu_funct*full_reg_funct$sigma^2
  #Stat
  stat_funct=as.numeric((bias_funct/full_reg_funct$sigma)^2/(g_star_funct%*%xtxm1_funct%*%g_star_funct))
  #Pvalue
  pval_finalSampleBased = pf(stat_funct, df1=1, df2=nrow(df)-dim(x_funct)[2], ncp=1)
  
  #Full model with pilot
  if(pilot_based_selection){
    x_pilot_funct=as.matrix(cbind(1,pilot_df[,x_var])) #add intercept
    y_pilot_funct=as.matrix(pilot_df[,y_var]) 
    xtxm1_pilot_funct=solve(crossprod(x_pilot_funct))
    full_pilot_reg_funct=summary(lm(as.formula(fm), data=pilot_df))
    
    #Bias
    bias_pilot_funct=as.numeric(stats::coef(full_pilot_reg_funct)[,1]%*%g_star_funct)
    
    #Stat with pilot
    stat_withPilot_funct=as.numeric((bias_pilot_funct/full_pilot_reg_funct$sigma)^2/(g_star_funct%*%xtxm1_pilot_funct%*%g_star_funct))
    #Non-centrality parameter corrected
    lambda_pilot_corrected=(g_star_funct%*%xtxm1_funct%*%g_star_funct)/(g_star_funct%*%xtxm1_pilot_funct%*%g_star_funct)
    #Pvalue
    pval_pilotBased = pf(stat_withPilot_funct, df1=1, df2=nrow(pilot_df)-dim(x_funct)[2], ncp=lambda_pilot_corrected)
  }

  #Selection with POINT prior beliefs
  if(priorPoint_based_selection){
    #Stat with point priors
    priors_mu_star=c(0,0,priors_mu)
    stat_pointPrior=as.numeric(full_reg_funct$sigma^2*(nu_3-nu_0)/((g_star_funct%*%priors_mu_star)^2)*(nrow(df)-ncol(x_funct)))
    pval_pointPriorBased=1-pchisq(stat_pointPrior, df=nrow(df)-dim(x_funct)[2])
  }

  listReturn=list(cm=cm,
                  estimate=estimate_funct, 
                  bias=bias_funct,
                  bias_sqrd_corr=bias_sqrd_corrected_funct,
                  var=var_funct, 
                  MSE=bias_sqrd_corrected_funct+var_funct, 
                  correction_unscaled=as.numeric(correction_funct), 
                  nu=nu_funct,
                  full_sigma2=full_reg_funct$sigma^2)
  if(finalSample_based_selection){
    listReturn$stat_finalSample=stat_funct
    listReturn$pval_finalSampleBased=pval_finalSampleBased
  }
  if(pilot_based_selection){
    listReturn$stat_withPilot=stat_withPilot_funct
    listReturn$lambda_ncp_pilot=as.numeric(lambda_pilot_corrected)
    listReturn$pval_pilotBased=pval_pilotBased
  }
  if(priorPoint_based_selection){
    listReturn$stat_pointPrior=stat_pointPrior
    listReturn$pval_pointPriorBased=round(pval_pointPriorBased,10)
  }
  
  return(listReturn)
}

