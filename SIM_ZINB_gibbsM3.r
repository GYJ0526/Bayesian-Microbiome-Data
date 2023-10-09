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
  alpha0<-beta0<-numeric(p)

## MCMC samples
  Beta<-Alpha<-matrix(0,lastit,p)
  R2<-matrix(0,lastit,1)
  Sigma2<-matrix(0,lastit,1)
  XI<-matrix(0,lastit,N)

## Priors
  T0a<-diag(0.01,p)
  T0b<-diag(0.01,p)

## Initial values
  y1<-rbinom(N,1,0.5)
  y1[y>0]<-1
  N0<-length(y[y==0])
  q<-rep(0.5,N)
  id<-1:N
  Acc<-0
  s<-0.003

  r2<-1
  alpha<-rnorm(p,0,1)
  beta <-rnorm(p,0,1)
  zeta<-rnorm(sum(y1),0,1)
  sigma2<-1
  xi<-rnorm(N,0,sd=sqrt(sigma2))

## Running MCMC
  for (i in 1:nrep){
  
## Update alpha
    Xalpha<-X%*%alpha
    w1<-rpg(N,1,Xalpha)                           	
    z1<-(y1-1/2)/w1
    va<-solve(crossprod(X*sqrt(w1))+T0a)
    ma<-va%*%(T0a%*%alpha0+t(sqrt(w1)*X)%*%(sqrt(w1)*(z1)))
    alpha<-c(rmvnorm(1,ma,va))

## Update latent class indicator y1
    Xa<-X%*%alpha
    phi<-1/(1+exp(-Xa))
    Xb<-X%*%beta+xi
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
    Xbeta<-X[y1==1,]%*%beta+xi[y1==1]
    w2<-rpg(N1,y[y1==1]+r2,Xbeta)                 	# Polya weights
    w2[w2<0.001]<-0.001
    z2<-(y[y1==1]-r2)/(2*w2)                      	# Latent "response"
    vb<-solve(crossprod(X[y1==1,]*sqrt(w2))+T0b)
    mb<-vb%*%(t(sqrt(w2)*X[y1==1,])%*%(sqrt(w2)*(z2-xi[y1==1])))
    beta<-c(rmvnorm(1,mb,vb))

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
    eta<-X[y1==1,]%*%beta+xi[y1==1]
    psi<-1/(1+exp(-eta))
    psi[psi<0.001]<-0.001
    psi[psi>0.9999]<-0.999
    r2<-rgamma(1,0.01+sum(l),0.01-sum(log(1-psi)))

# Update the xi for random effects
    N2<-length(nis[nis>0])
    priorprec<-sigma2
    priormean<-X[y1==0,]%*%beta
    xi[nis==0]<-rnorm(N-N2,priormean,sqrt(priorprec))
    vb<-1/(sigma2+tapply(w2,id[y1==1],sum))
    mb<-vb*(tapply(w2*(z2-X[y1==1,]%*%beta),id[y1==1],sum))
    xi[nis>0]<-rnorm(N2,mb,sqrt(vb))
#   xi<-xi-mean(xi)

# Update the variance of the random effect sigma2
    inv.a<-0.01
    inv.b<-0.01
    para1<-inv.a+N2/2
    para2<-inv.b+(t(xi[nis>0])%*%(xi[nis>0]))/2
    Sigmab<-rinvgamma(1,shape=para1,scale=para2)
    sigma2<-Sigmab

# Store
    if (i> burn & i%%thin==0) {
      k<-(i-burn)/thin
      Alpha[k,]<-alpha
      Beta[k,]<-beta
      R2[k,]<-r2
      Sigma2[k,]<-sigma2
      XI[k,]<-xi
    }
  
    if (i%%100==0) print(i)
  
  }

  OUT<-list(y,X,Alpha,Beta,R2,XI,Sigma2)
  names(OUT)<-c("y","X","alpha","beta","phi","xi","sigma2")
  return(OUT)

}