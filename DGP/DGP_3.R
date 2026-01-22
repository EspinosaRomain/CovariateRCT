#DGP #3: Multivariate test: exact test
beta=1
N=200
rho_vec=c(0.5,0.5,0.5)
corr_other=0.2
CorrMat=matrix(data=corr_other, nrow = 4, ncol = 4)
CorrMat[2:4,1]=rho_vec
CorrMat[1,2:4]=rho_vec
diag(CorrMat)=1
CorrMat
SD_vec=c(1,rep(0.8,3))
VarMat=diag(SD_vec)%*%CorrMat%*%diag(SD_vec)
VarMat

#Computes nus
nu_0=solve(VarMat[1:1,1:1])[1,1]/N
nu_0
nu_1=solve(VarMat[1:2,1:2])[1,1]/N
nu_1
nu_2=solve(VarMat[1:3,1:3])[1,1]/N
nu_2
nu_3=solve(VarMat[1:4,1:4])[1,1]/N
nu_3

#Computes g's
g0=c(solve(VarMat[1:1,1:1])%*%VarMat[(1:1),(2:4)])
g1=c(0,(solve(VarMat[1:2,1:2])%*%VarMat[(1:2),(3:4)])[1,])
g2=c(0,0,(solve(VarMat[1:3,1:3])%*%VarMat[(1:3),(4:4)])[1,])
g3=rep(0,3)

ESvec=c(0.4,0.2,0.05)
#I set effect sizes such that MSE2 and MSE3 are equal in theory
sigma2=uniroot(function(s2_funct){
  MSE2_funct=as.double((g2%*%ESvec)^2+nu_2*s2_funct)
  MSE3_funct=as.double(nu_3*s2_funct)
  MSE2_funct-MSE3_funct
}, lower = 0, upper = 1000, tol=.Machine$double.eps)$root
sigma2

theoreticalMSE_0=as.double((g0%*%ESvec)^2+nu_0*sigma2)
theoreticalMSE_0
theoreticalMSE_1=as.double((g1%*%ESvec)^2+nu_1*sigma2)
theoreticalMSE_1
theoreticalMSE_2=as.double((g2%*%ESvec)^2+nu_2*sigma2)
theoreticalMSE_2
theoreticalMSE_3=as.double(nu_3*sigma2)
theoreticalMSE_3



