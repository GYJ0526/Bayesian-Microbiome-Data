#####################
# Posterior summary #
#####################

Post.summary<-function(Data){

  alpha<-Data$alpha
  beta<-Data$beta
  phi<-Data$phi
  sigma2<-Data$sigma2
  X<-Data$X
  y<-Data$y
  N<-dim(X)[1]
  B<-dim(alpha)[1]
  p<-dim(alpha)[2]

  m.alpha<-m.beta<-numeric(p)
  s.alpha<-s.beta<-numeric(p)
  m.R2<-s.R2<-0
  m.sigma2<-s.sigma2<-0
  L.a<-U.a<-L.b<-U.b<-numeric(p)
  L.r<-U.r<-0
  L.s<-U.s<-0
  L.g<-U.g<-0

## Compute marginalization for the group effect
  theta<-numeric(B)
  theta.0<-num<-den<-0
  z.a<-z.b<-c.a<-c.b<-0
  id.diab<-X[,2]
  
  OR.ft<-function(arg1,arg2,arg3){

    z.a<-exp(arg1%*%arg2)/(1+exp(arg1%*%arg2))
    z.b<-exp(arg1%*%arg3)
    theta.0<-z.a*z.b 
    return(theta.0)   
    
  }  

  id.0<-which(id.diab==0)
  L0<-length(id.0)
  id.1<-which(id.diab==1)
  L1<-length(id.1)

  theta.i<-numeric(L0)
  num<-numeric(B)
  for (b in 1:B){
    for (i in 1:length(id.0))
      theta.i[i]<-OR.ft(arg1=X[id.0[i],],arg2=alpha[b,],arg3=beta[b,])
    num[b]<-mean(theta.i)
#    num[b]<-sum(theta.i)
  }

  theta.i<-numeric(L1)
  den<-numeric(B)
  for (b in 1:B){
    for (i in 1:length(id.1))
      theta.i[i]<-OR.ft(arg1=X[id.1[i],],arg2=alpha[b,],arg3=beta[b,])
    den[b]<-mean(theta.i)
#    den[b]<-sum(theta.i)
  }

  OR<-den/num

### 95% HPD Interval
  library(boa)
  m.alpha<-colMeans(alpha)
  m.beta <-colMeans(beta)
  m.R2<-colMeans(phi)
  s.R2<-sd(phi)
  m.sigma2<-colMeans(sigma2)
  s.sigma2<-sd(sigma2)
  m.gamma<-mean(OR)
  s.gamma<-sd(OR)
  L.r<-as.numeric(boa.hpd(phi,alpha=0.05)[1])
  U.r<-as.numeric(boa.hpd(phi,alpha=0.05)[2])
  L.s<-as.numeric(boa.hpd(sigma2,alpha=0.05)[1])
  U.s<-as.numeric(boa.hpd(sigma2,alpha=0.05)[2])
  L.g<-as.numeric(boa.hpd(OR,alpha=0.05)[1])
  U.g<-as.numeric(boa.hpd(OR,alpha=0.05)[2])
  
  for (j in 1:p){
    s.alpha[j]<-sd(alpha[,j])
    s.beta[j]<-sd(beta[,j])
    L.a[j]<-as.numeric(boa.hpd(alpha[,j],alpha=0.05)[1])
    U.a[j]<-as.numeric(boa.hpd(alpha[,j],alpha=0.05)[2])
    L.b[j]<-as.numeric(boa.hpd(beta[,j],alpha=0.05)[1])
    U.b[j]<-as.numeric(boa.hpd(beta[,j],alpha=0.05)[2])
  }

### Result table
  Est<-round(c(m.alpha,m.beta,m.R2,m.sigma2,m.gamma),3)
  SD<-round(c(s.alpha,s.beta,s.R2,s.sigma2,s.gamma),3)
  L.HPD<-round(c(L.a,L.b,L.r,L.s,L.g),3)
  U.HPD<-round(c(U.a,U.b,U.r,U.s,U.g),3)

  OUT<-list(Est,SD,L.HPD,U.HPD)
  names(OUT)<-c("Est","SD","LHPD","UHPD")
  return(OUT)

}
