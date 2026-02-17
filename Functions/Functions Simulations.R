library(mvtnorm)

generateDatasetMulti=function(N_funct,EStreatment_funct,EScontrol_funct,CorrMatrix_funct,sigma2_funct,SD_vec_funct){
  k_funct=length(EScontrol_funct)
  VatMatrix_funct=diag(SD_vec_funct)%*%CorrMatrix_funct%*%diag(SD_vec_funct)
  RV=cbind(c(rep(0,round(N_funct/2)),rep(1,round(N_funct/2))),sapply(1:k_funct, function(i) rnorm(N_funct)))
  RV=sapply(1:(k_funct+1), function(i) RV[,i]-mean(RV[,i]))
  RVmod=RV%*%solve(chol(crossprod(RV)/N_funct))%*%chol(VatMatrix_funct)
  df_funct=data.frame(id=seq(1:N_funct))
  df_funct$treatment=RVmod[,1]
  total_controls=k_funct
  for(k_loop in 1:total_controls){
    var_funct=paste0("X",k_loop)
    df_funct[[paste0(var_funct,"")]]=RVmod[,(k_loop+1)]
  }
  noise_funct=rnorm(N_funct,mean=0,sd=sqrt(sigma2_funct))
  df_funct$y=(df_funct$treatment*EStreatment_funct+as.matrix(df_funct[,-c(1,2)])%*%EScontrol_funct)+noise_funct
  df_funct$y=df_funct$y-mean(df_funct$y)
  return(df_funct)
}
