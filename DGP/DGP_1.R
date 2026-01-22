#DGP #1: Univariate case

#Set Parameters
N=200
beta=1
sigma2=1
rho=0.3
CorrMat=matrix(data=c(1,rho,rho,1),nrow=2,ncol=2)
gamma=sqrt(1/(N*(1-rho^2)))
ESvec=c(gamma)
