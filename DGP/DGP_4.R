#DGP #4: Selection rules
beta=1
N=200
if(!exists("EffectSizeMSE")) EffectSizeMSE=1
rho_vec=rep(0.5,3)
corr_other=0
CorrMat=matrix(data=corr_other, nrow = 4, ncol = 4)
CorrMat[2:4,1]=rho_vec
CorrMat[1,2:4]=rho_vec
diag(CorrMat)=1
CorrMat
SD_vec=c(1,rep(1,3))
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

sigma2=1.5
ES_scalar=uniroot(function(x_funct){
  ESvec_loop=rep(x_funct,3)
  MSE0_funct=as.double((g0%*%ESvec_loop)^2+nu_0*sigma2)
  MSE3_funct=as.double(nu_3*sigma2)
  MSE0_funct-EffectSizeMSE*MSE3_funct
}, lower = 0, upper = 1, tol=.Machine$double.eps)$root
ES_scalar
ESvec=rep(ES_scalar,3)

theoreticalMSE_0=as.double((g0%*%ESvec)^2+nu_0*sigma2)
theoreticalMSE_0
theoreticalMSE_1=as.double((g1%*%ESvec)^2+nu_1*sigma2)
theoreticalMSE_1
theoreticalMSE_2=as.double((g2%*%ESvec)^2+nu_2*sigma2)
theoreticalMSE_2
theoreticalMSE_3=as.double(nu_3*sigma2)
theoreticalMSE_3
