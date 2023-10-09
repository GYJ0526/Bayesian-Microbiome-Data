###########################################
# Bayesian modelingon Microbiome Project  #
# Simulation Study code				#
###########################################
library(BayesLogit)   # For rpg function -- install from ZIP file
library(mvtnorm)	
library(MCMCpack)     # For Iwish update of Sigmab
library(msm)          # For tnorm function

setwd('/Users/yeongjin.gwon/Dropbox/Mac/Desktop/Simulation/ZINBR')

### Basic simulation setting
SIM<-100
Sim1.eval<-matrix(0,SIM,4)

for (i in 1:SIM){
  cat("Running ith =",i,"\n")
  ### Data Generation
  source("DataGen.r")
  OUT<-Data.Gen(a0=0.5,a1=0,b0=2,b1=0,tphi=0.8,tsig2=1.5,iter=i,Dist=1)
  
  ### Display of the generated data
  y<-OUT$y
  N<-length(y)
  LogT<-OUT$LogT
  X<-OUT$X
  tmp<-table(y)/N*100
  barplot(tmp)

### Gibbs sampling
  source("SIM_ZINB_gibbsM3-HC.r")
  RES<-MCMC.M3(Data1=y,Data2=X,T=LogT,NREP=6000,NTHIN=5,NBURN=1000)

### Estimate overall effect
  source("SIM_ZINB_PostsummaryM3.r")
  OUT<-Post.summary(RES)

### Evaluation of performance
  Sim1.eval[i,]<-round(c(OUT$Est[7],OUT$SD[7],OUT$LHPD[7],OUT$UHPD[7]),3)

}

write.csv(Sim1.eval,file='1_C1_SimRes_M3_Error_N500.csv')