### Deviance Information Criteria ###
DIC.M3<-function(Data,ID,Model){
	
  PostSample<-Data
  y<-PostSample$y
  X<-PostSample$X
  Alpha<-PostSample$alpha
  Beta<-PostSample$beta
  Phi<-PostSample$phi
  Sigma2<-PostSample$sigma2
  XI<-PostSample$xi
  N<-dim(X)[1]
  B<-dim(Alpha)[1]
  p<-dim(Alpha)[2]

## Zero and count components are coded seperately
  Dev1<-0
  OUT<-numeric(N)
  L<-numeric(B)

  for (i in 1:N){
    for (b in 1:B){
      pi<-mu<-cnt1.0<-cnt1.c<-cnt2.c<-0
      pi<-as.numeric(exp(X[i,]%*%Alpha[b,])/(1+exp(X[i,]%*%Alpha[b,])))
      mu<-as.numeric(exp(X[i,]%*%Beta[b,]+XI[b,i]))
      if (y[i]==0) {
        cnt1.0<-log(1-pi+pi*(Phi[b]/(Phi[b]+mu))^Phi[b])    
      } else {
        cnt1.c<-log(pi)+lgamma(y[i]+Phi[b])-lgamma(Phi[b])-lgamma(y[i]+1)
        cnt2.c<-Phi[b]*(log(Phi[b])-log(Phi[b]+mu))+y[i]*(log(mu)-log(Phi[b]+mu))
      }
      L[b]<--2*(cnt1.0+cnt1.c+cnt2.c)
    }
    OUT[i]<-mean(L)
  }
  Dev1<-sum(OUT)

## Compute a goodness-of-fit measure DICbar
  Dev2<-0
  OUT<-numeric(N)

  for (i in 1:N){
    pi<-mu<-cnt1.0<-cnt1.c<-cnt2.c<-0
    pi<-as.numeric(exp(X[i,]%*%colMeans(Alpha))/(1+exp(X[i,]%*%colMeans(Alpha))))
    mu<-as.numeric(exp(X[i,]%*%colMeans(Beta)+mean(XI[,i])))
    if (y[i]==0) {
      cnt1.0<-log(1-pi+pi*(mean(Phi)/(mean(Phi)+mu))^mean(Phi))
    } else {
      cnt1.c<-log(pi)+lgamma(y[i]+mean(Phi))-lgamma(mean(Phi))-lgamma(y[i]+1)
      cnt2.c<-mean(Phi)*(log(mean(Phi))-log(mean(Phi)+mu))+y[i]*(log(mu)-log(mean(Phi)+mu))
    }
    OUT[i]<--2*(cnt1.0+cnt1.c+cnt2.c)
  }
  Dev2<-sum(OUT)

## Compute dimension penalty pD
  pD<-Dev1-Dev2
  DIC<-Dev2+2*pD
  list_out<-list(Dev1,Dev2,pD,DIC)
  names(list_out)<-c("Deviance","Goodness-of-fit","pD","DIC")

  setwd('/Users/yeongjin.gwon/OneDrive - University of Nebraska Medical Center/01 Academics/01 Research/05 Microbiome/4. Output/Selected taxa')
  write.csv(list_out,file=sprintf("DIC_y%d_M%d.csv",ID,Model))
  
}
