### Data Generation Code ###
Data.Gen<-function(a0,a1,b0,b1,tphi,tsig2,iter,Dist){

  t.alpha<-c(a0,a1)
  t.beta<-c(b0,b1)
  #  t.alpha<-c(-0.5,0.2)
  #  t.beta<-c(2,0.5)
  t.phi<-tphi
  t.sigma2<-tsig2

### Data Generation 
  set.seed(997750375+iter)
  N<-500
  x1<-rbinom(N,1,0.5)
  X<-cbind(rep(1,N),x1)
  p<-ncol(X)

### Zero part
  eta1<-X%*%t.alpha
  phi<-1/(1+exp(-eta1))
  uu<-rbinom(N,1,phi)
  N1<-sum(uu)
  p0<-(N-N1)/N

### Count part
  LogT<-runif(N,4,6)
  eta2<-X[uu==1,]%*%t.beta+rnorm(N1,0,sd=sqrt(t.sigma2))+LogT[uu=1]
  if (Dist==1) {
    mu<-exp(eta2)
    y<-rep(0,N)
    y[uu==1]<-rpois(N1,mu)
  } else {
    pi<-1/(1+exp(-eta2))
    mu<-t.phi*pi/(1-pi)
    y<-rep(0,N)
    y[uu==1]<-rnbinom(N1,t.phi,mu=mu)   
  }
  pzero<-length(y[y==0])/N
  OUT<-list(y=y,X=X,LogT=LogT)
  return(OUT)
  
}

#OUT<-Data.Gen(a0=-0.5,a1=0.1,b0=0.5,b1=0.2,tphi=0.8,tsig2=1.5,iter=1,Dist=1)
