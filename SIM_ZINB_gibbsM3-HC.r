MCMC.M3<-function(Data1,Data2,T,NREP,NTHIN,NBURN){

  y<-Data1
  X<-Data2
  LogT<-T
  N<-dim(X)[1]
  p<-dim(X)[2]

  nrep<-NREP					
  thin<-NTHIN		           
  burn<-NBURN    		     	
  lastit<-(nrep-burn)/thin
  alpha<-numeric(p)
  beta<-numeric(p)

## Initial values
  y1<-rbinom(N,1,0.5)
  y1[y>0]<-1
  N0<-length(y[y==0])
  q<-rep(0.5,N)
  id<-1:N
  Acc<-0
  s<-0.003

  N1<-sum(y1)
  r2<-1
  alpha<-rnorm(p,0,1)
  beta <-rnorm(p,0,1)
  zeta<-rnorm(N,0,1)
  sigma2<-1

## MCMC samples
  Beta<-Alpha<-matrix(0,lastit,p)
  R2<-matrix(0,lastit,1)
  Sigma2<-matrix(0,lastit,1)
  Zeta<-matrix(0,lastit,N)

## Priors -- note diffuse uniform for r
  T0a<-diag(0.01,p)          	# Prior precision for alpha -- tightening can help in bordeline cases
  T0b<-diag(0.01,p)          	# Prior precision for beta

## Running MCMC
  for (i in 1:nrep){
  
## Update alpha
    Xalpha<-X%*%alpha
    w1<-rpg(N,1,Xalpha)                           	
    z1<-(y1-1/2)/w1
    va<-solve(crossprod(X*sqrt(w1))+T0a)
    ma<-va%*%(t(sqrt(w1)*X)%*%(sqrt(w1)*(z1)))
    alpha<-c(rmvnorm(1,ma,va))

## Update latent class indicator y1
    Xa<-X%*%alpha
    phi<-1/(1+exp(-Xa))
    Xb<-zeta
    q<-1/(1+exp(Xb))
    num<-phi*(q^r2)
    den<-phi*(q^r2)+1-phi
    den[den==0]<-0.001
    theta<-num/den
    theta[theta>1]<-1
    y1[y==0]<-rbinom(N0,1,theta[y==0]) 
    N1<-sum(y1)
    nis<-tapply(y1,id,sum) 

## Update beta
    vb<-solve(1/sigma2*t(X[y1==1,])%*%X[y1==1,]+T0b)
    temp<-as.vector(t(X[y1==1,])%*%as.matrix(zeta[y1==1]))
    mb<-vb%*%(temp/sigma2)
    beta<-c(rmvnorm(1,mb,vb))

## Update zeta
    N2<-length(nis[nis>0])
#   priormean<-X[y1==0,]%*%beta
    priormean<-0
    zeta[nis==0]<-rnorm(N-N2,priormean,sd=sqrt(sigma2))
    w2<-rpg(N2,y[y1==1]+r2,zeta[y1==1])                 	
    w2[w2<0.001]<-0.001
    z2<-(y[y1==1]-r2)/(2*w2)  
    vz<-1/(1/sigma2+w2)
    mz<-vz*(1/sigma2*as.vector(X[y1==1,]%*%beta)+z2*w2)
    zeta[nis>0]<-rnorm(N2,mz,sqrt(vz))  

# Update r2
#  rnew<-rtnorm(1,r2,sqrt(s),lower=0)
#  tmp1<-sum(dnbinom(y[y1==1],rnew,q[y1==1],log=T))
#  tmp2<-sum(dnbinom(y[y1==1],r2,q[y1==1],log=T))
#  tmp3<-dtnorm(r2,rnew,sqrt(s),0,log=T)
#  tmp4<-dtnorm(rnew,r2,sqrt(s),0,log=T)
#  ratio<-tmp1-tmp2+tmp3-tmp4
#  U<-runif(1,0,1)
#  if (log(U)<ratio) {
#    r2<-rnew
#    Acc<-Acc+1
#  }

# Update r2 using Gibbs as in Dadaneh et al and Zhou and Carin #
# Update latent counts, l        
    l<-rep(0,N1)
    ytmp<-y[y1==1]
    for(k in 1:N1)
      l[k]<-sum(rbinom(ytmp[k],1,round(r2/(r2+1:ytmp[k]-1),6)))
  
# Update r from conjugate gamma distribution given l and psi
    eta<-zeta[y1==1]
    psi<-1/(1+exp(-eta))
    psi[psi<0.0001]<-0.0001
    psi[psi>0.9999]<-0.9999
    r2<-rgamma(1,0.01+sum(l),0.01-sum(log(1-psi)))

# Update the variance of the random effect sigma2
    inv.a<-0.01
    inv.b<-0.01
    temp<-zeta[y1==1]-X[y1==1,]%*%beta
#   temp<-zeta-X%*%beta
    para1<-inv.a+N2/2
    para2<-inv.b+t(temp)%*%temp/2
    Sigmab<-rinvgamma(1,shape=para1,scale=para2)
    sigma2<-Sigmab

# Store
    if (i> burn & i%%thin==0) {
      k<-(i-burn)/thin
      Alpha[k,]<-alpha
      Beta[k,]<-beta
      R2[k,]<-r2
      Sigma2[k,]<-sigma2
      Zeta[k,]<-zeta
    }
  
    if (i%%1000==0) print(i)
  
  }

  OUT<-list(y,X,Alpha,Beta,Zeta,R2,Sigma2)
  names(OUT)<-c("y","X","alpha","beta","zeta","phi","sigma2")
  return(OUT)

}