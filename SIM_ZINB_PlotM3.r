###############
# Trace plots #
###############
# Binary part alpha0-alpha2 for each taxa

Trace.plot<-function(Data,NREP,ID,Model){
	
  Postsample<-Data
  Alpha<-Postsample$alpha
  Beta<-Postsample$beta
  Phi<-Postsample$phi
  Sigma2<-Postsample$sigma2
  Xi<-Postsample$xi
  lastit<-NREP
   
  setwd('/Users/yeongjin.gwon/OneDrive - University of Nebraska Medical Center/01 Academics/01 Research/05 Microbiome/4. Output/Selected taxa')

  pdf(sprintf("Traceplot_y%d_M%d.pdf",ID,Model))	
  par(mfrow=c(3,3))
  plot(1:lastit,Alpha[,1],type="l",col="lightgreen",xlab="Iteration",ylab=expression(alpha[10]))
  abline(h=mean(Alpha[,1]),col="blue4")
  
  plot(1:lastit,Alpha[,2],type="l",col="lightgreen",xlab="Iteration",ylab=expression(alpha[11]))
  abline(h=mean(Alpha[,2]),col="blue4")

  plot(1:lastit,Alpha[,3],type="l",col="lightgreen",xlab="Iteration",ylab=expression(alpha[12]))
  abline(h=mean(Alpha[,3]),col="blue4")

  plot(1:lastit,Alpha[,4],type="l",col="lightgreen",xlab="Iteration",ylab=expression(alpha[13]))
  abline(h=mean(Alpha[,4]),col="blue4")
  
  plot(1:lastit,Alpha[,5],type="l",col="lightgreen",xlab="Iteration",ylab=expression(alpha[14]))
  abline(h=mean(Alpha[,5]),col="blue4")

  plot(1:lastit,Alpha[,6],type="l",col="lightgreen",xlab="Iteration",ylab=expression(alpha[15]))
  abline(h=mean(Alpha[,6]),col="blue4")

  plot(1:lastit,Alpha[,7],type="l",col="lightgreen",xlab="Iteration",ylab=expression(alpha[16]))
  abline(h=mean(Alpha[,7]),col="blue4")
  
  plot(1:lastit,Alpha[,8],type="l",col="lightgreen",xlab="Iteration",ylab=expression(alpha[17]))
  abline(h=mean(Alpha[,8]),col="blue4")
  
  plot(1:lastit,Alpha[,9],type="l",col="lightgreen",xlab="Iteration",ylab=expression(alpha[18]))
  abline(h=mean(Alpha[,9]),col="blue4")  

  par(mfrow=c(3,3))
  plot(1:lastit,Beta[,1],type="l",col="lightgreen",xlab="Iteration",ylab=expression(beta[10]))
  abline(h=mean(Beta[,1]),col="blue4")
  
  plot(1:lastit,Beta[,2],type="l",col="lightgreen",xlab="Iteration",ylab=expression(beta[11]))
  abline(h=mean(Beta[,2]),col="blue4")

  plot(1:lastit,Beta[,3],type="l",col="lightgreen",xlab="Iteration",ylab=expression(beta[12]))
  abline(h=mean(Beta[,3]),col="blue4")

  plot(1:lastit,Beta[,4],type="l",col="lightgreen",xlab="Iteration",ylab=expression(beta[13]))
  abline(h=mean(Beta[,4]),col="blue4")
  
  plot(1:lastit,Beta[,5],type="l",col="lightgreen",xlab="Iteration",ylab=expression(beta[14]))
  abline(h=mean(Beta[,5]),col="blue4")

  plot(1:lastit,Beta[,6],type="l",col="lightgreen",xlab="Iteration",ylab=expression(beta[15]))
  abline(h=mean(Beta[,6]),col="blue4")

  plot(1:lastit,Beta[,7],type="l",col="lightgreen",xlab="Iteration",ylab=expression(beta[16]))
  abline(h=mean(Beta[,7]),col="blue4")
  
  plot(1:lastit,Beta[,8],type="l",col="lightgreen",xlab="Iteration",ylab=expression(beta[17]))
  abline(h=mean(Beta[,8]),col="blue4")
  
  plot(1:lastit,Beta[,9],type="l",col="lightgreen",xlab="Iteration",ylab=expression(beta[18]))
  abline(h=mean(Beta[,9]),col="blue4")  
  
  plot(1:lastit,Phi[,1],type="l",col="lightgreen",xlab="Stored Iteration",ylab=expression(phi))
  abline(h=mean(Phi[,1]),col="blue4")

  plot(1:lastit,Sigma2[,1],type="l",col="lightgreen",xlab="Stored Iteration",ylab=expression(sigma[2]))
  abline(h=mean(Sigma2[,1]),col="blue4")
  dev.off()

}
