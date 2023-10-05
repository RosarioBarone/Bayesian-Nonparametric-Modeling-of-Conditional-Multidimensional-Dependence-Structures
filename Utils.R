##############################################################
#################### ROSARIO BARONE ##########################
#### Auxiliary functions DPM of conditional vine copulas #####
##############################################################

#### Inverse Fisher Transform ####
fish.inv <- function(beta){
  out<-(exp(2*beta) - 1) / (exp(2*beta) + 1)
  ### Fix numerical problems
  out[which(round(abs(out), digits = 3)==1)]=out[which(round(abs(out), digits = 3)==1)]*0.999 
  out[which(is.nan(out))]=0.999*sign(beta[which(is.nan(out))])
  return(out)
}

#### Conditional Gaussian bivariate copula density with linear link ####
dNormCondCopDens=function(u1,u2,b,x,log){
  theta<- x%*%b
  rho<-fish.inv(theta)
  x<-qnorm(u1)
  y<-qnorm(u2)
  cov = 1-rho^2
  cc = cov^(-0.5)*exp(-0.5*(rho^(2)*(x^2+y^2)-2*rho*x*y)*cov^(-1))
  if(isTRUE(log)) cc=log(cc)
  return(cc)
}

###  Conditional Gaussian vine copula density ####
dVineCopCond<-function(u,b,x,d,log){
  if(is.null(dim(u))) u<-matrix(u,nrow = 1,ncol = length(u))
  levels<-d-1
  density<-NULL
  #First tree-level
  j<-1
  par<-1
  while(j<d){
  if(dim(x)[1]>1) density<-cbind(density,dNormCondCopDens(u[,j],u[,j+1],b[,par],cbind(1,x[,,j]), log=TRUE) ) 
  if(dim(x)[1]==1) density<-cbind(density,dNormCondCopDens(u[,j],u[,j+1],b[,par],c(1,x[,,j]), log=TRUE) ) 
   j<-j+1
    par<-par+1
  }
  #Other tree-levels
  v<-c()  
  u.v<-NULL
  for (l in 2:levels){
    if(l>2)   u<-u.v
    j<-1
    while(j<=(d-l)){
      if(dim(x)[1]>1){
        v<-BiCopHfunc2(u[,j],u[,j+1],1,par=fish.inv(cbind(1,x[,,j])%*%b[,j]))
        v<-cbind(v,BiCopHfunc2(u[,j+2],u[,j+1],1,par=fish.inv(cbind(1,x[,,j+1])%*%b[,j+1])))
        density<-cbind(density,dNormCondCopDens(v[,1],v[,2],b[,par],cbind(1,x[,,par]), log=TRUE)) 
        par<-par+1
        u.v<-cbind(u.v,v)
        j=j+1
        }else{
        v<-BiCopHfunc2(u[,j],u[,j+1],1,par=fish.inv(c(1,x[,,j])%*%b[,j]))
        v<-cbind(v,BiCopHfunc2(u[,j+2],u[,j+1],1,par=fish.inv(c(1,x[,,j+1])%*%b[,j+1])))
        density<-cbind(density,dNormCondCopDens(v[,1],v[,2],b[,par],c(1,x[,,par]), log=TRUE)) 
        par<-par+1
        u.v<-cbind(u.v,v)
        j=j+1
      }
    }
  }
  out<-rowSums(density)
  if(isFALSE(log)) out<-exp(out)
  return(out)
}


####Gaussian bivariate copula density####
dNormCopDens=function(u1,u2,beta,log){
  rho<-fish.inv(beta)
  x<-qnorm(u1)
  y<-qnorm(u2)
  cov = 1-rho^2
  cc = cov^(-0.5)*exp(-0.5*(rho^(2)*(x^2+y^2)-2*rho*x*y)*cov^(-1))
  if(isTRUE(log)) cc=log(cc)
  return(cc)
}

dVineCop<-function(u,b,d,log){
  if(is.null(dim(u))) u<-matrix(u,nrow = 1,ncol = length(u))
  levels<-d-1
  density<-NULL
  #First tree-level
  j<-1
  par<-1
  while(j<d){
    density<-cbind(density,dNormCopDens(u[,j],u[,j+1],b[par], log=TRUE) ) 
    j<-j+1
    par<-par+1
  }
  #Other tree-levels
  v<-c()  
  u.v<-NULL
  for (l in 2:levels){
    if(l>2)   u<-u.v
    j<-1
    while(j<=(d-l)){
        v<-BiCopHfunc2(u[,j],u[,j+1],1,par=fish.inv(cbind(b[j])))
        v<-cbind(v,BiCopHfunc2(u[,j+2],u[,j+1],1,par=fish.inv(b[j+1])))
        density<-cbind(density,dNormCopDens(v[,1],v[,2],b[par], log=TRUE)) 
        par<-par+1
        u.v<-cbind(u.v,v)
        j=j+1
    }
  }
  out<-rowSums(density)
  if(isFALSE(log)) out<-exp(out)
  return(out)
}


# Relabel function useful in case of label switching

relabel.condVineDPM<-function(U,X,ncov,family_x,sim){
  take<-sim$take
  Psi<-sim$Psi
  Psi.mode<-sim$Psi.mode
  beta<-sim$beta
  XBin<-sum(family_x=="Binomial")
  XNorm<-sum(family_x=="Normal")
  d<-sim$d
  dd<-d*(d-1)/2
  if(XBin>0) theta<-sim$theta
  if(XNorm>0) mu<-sim$mu
  if(XNorm>0) sigma2<-sim$sigma2
  k<-Psi.mode
  library(extraDistr)
  MCMCiter<-length(take)
  w<-NULL
  for (i in take){
    psi<-Psi[,i]
    if(Psi.mode==2){
      w_i<-rbeta(1,shape1 = 1+length(which(psi==1)),shape2=length(psi)-length(which(psi==1)+1))
      w<-rbind(w,c(w_i,1-w_i))
    }else{
      w<-rbind(w,rdirichlet(1,1+as.vector(table(psi))))  
    }
  }
  apply(w, 2, mean)
  beta<-beta[,,,take]
  Psi<-Psi[,take]
  mu<-mu[take,,]
  sigma2<-sigma2[take,,]
  theta<-theta[take,,]
  if(is.na(dim(mu)[3])) dim(mu)[3]<-1
  if(is.na(dim(sigma2)[3])) dim(sigma2)[3]<-1
  if(is.na(dim(theta)[3])) dim(theta)[3]<-1
  
  simu<-list()
  simu$beta<-beta
  simu$mu<-mu
  simu$sigma2<-sigma2
  simu$theta<-theta
  simu$k<-Psi.mode
  simu$w<-w
  
  
  log.post=function(simu,U){
    n=dim(U)[1]
    beta<-simu$beta;
    k=simu$k; w=simu$w
    logpost=rep(0,dim(beta)[4])
    for (it in 1:dim(beta)[4]){
      lik=matrix(0,n,k)
      for(comp in 1:k) lik[,comp]=log(w[it,comp])+dVineCopCond(U,beta[,,k,it],X,d,log=TRUE)+dCovariates(X=X,family_x=family_x,mu=mu[it,k,],sigma2=sigma2[it,k,],theta=theta[it,k,])
      logpost[it]=sum(rowSums(lik))+sum(dnorm(beta[,,k,it],0,1,log = TRUE))
    }
    logpost}
   lp=log.post(simu,U)
  
  ind.map=order(lp,decreasing=TRUE)[1]
  map=list(beta=simu$beta[,,,ind.map],w=simu$w[ind.map,],mu=simu$mu[ind.map,,],sigma2=simu$sigma2[ind.map,,],theta=simu$theta[ind.map,,])
  pz.x.star=matrix(0,dim(U)[1],k)
  for (i in 1:dim(U)[1]){
    for (j in 1:k){
      if(XBin>1&XNorm>1) pz.x.star[i,j]=map$w[j]*dVineCopCond(U[i,],map$beta[,,j],array(X[i,,],dim = c(1,ncov,dd)),d,log=FALSE)+dCovariates(X=array(X[i,,],dim = c(1,ncov,dd)),family_x=family_x,mu=map$mu[j,],sigma2=map$sigma2[j,],theta=map$theta[j,])
      if(XBin<2&XNorm<2) pz.x.star[i,j]=map$w[j]*dVineCopCond(U[i,],map$beta[,,j],array(X[i,,],dim = c(1,ncov,dd)),d,log=FALSE)+dCovariates(X=array(X[i,,],dim = c(1,ncov,dd)),family_x=family_x,mu=map$mu[j],sigma2=map$sigma2[j],theta=map$theta[j])
      if(XBin<2&XNorm>1) pz.x.star[i,j]=map$w[j]*dVineCopCond(U[i,],map$beta[,,j],array(X[i,,],dim = c(1,ncov,dd)),d,log=FALSE)+dCovariates(X=array(X[i,,],dim = c(1,ncov,dd)),family_x=family_x,mu=map$mu[j,],sigma2=map$sigma2[j,],theta=map$theta[j])
      if(XBin<1&XNorm>2) pz.x.star[i,j]=map$w[j]*dVineCopCond(U[i,],map$beta[,,j],array(X[i,,],dim = c(1,ncov,dd)),d,log=FALSE)+dCovariates(X=array(X[i,,],dim = c(1,ncov,dd)),family_x=family_x,mu=map$mu[j],sigma2=map$sigma2[j],theta=map$theta[j,])
     }
    pz.x.star[i,]=pz.x.star[i,]/sum(pz.x.star[i,])
  }
  or.beta=array(0,dim=c(ncov+1,dd,k,MCMCiter))
  or.w<-matrix(0,MCMCiter,k); or.Psi<-matrix(0,dim(U)[1],MCMCiter)
  or.mu =array(0,dim = c(MCMCiter,k,XNorm))
  or.sigma2=array(0,dim = c(MCMCiter,k,XNorm))
  or.theta=array(0,dim = c(MCMCiter,k,XBin))
  pz.x=matrix(0,dim(U)[1],k)
  library(combinat)
  perma=permn(k)
  for ( it in 1:MCMCiter){
    if(it%%100==0){
      print(paste("Relabelling: iteration n",it))
    }
    kl.dist=rep(0,factorial(k))
    for (i in 1:dim(U)[1]){
      for(comp in 1:k){
        pz.x[i,comp]=w[it,comp]*dVineCopCond(U[i,],beta[,,comp,it],array(X[i,,],dim = c(1,ncov,dd)),d,log=FALSE)*exp(dCovariates(X=array(X[i,,],dim = c(1,ncov,dd)),family_x=family_x,mu=mu[it,comp,],sigma2=sigma2[it,comp,],theta=theta[it,comp,]))
      }
      pz.x[i,]=pz.x[i,]/sum(pz.x[i,])
      for (tau in 1:factorial(k)) kl.dist[tau]=kl.dist[tau]+sum(pz.x.star[i,]*log(pz.x[i,perma[[tau]]]))}
    best=order(kl.dist,decreasing=TRUE)[1]
    or.beta[,,,it]=beta[,,perma[[best]],it]
    or.mu[it,,]=mu[it,perma[[best]],]
    or.sigma2[it,,]=or.sigma2[it,perma[[best]],]
    or.theta[it,,]=or.theta[it,perma[[best]],]
    or.w[it,]=w[it,perma[[best]]]
  }
  results<-list()
  results$weights<-or.w
  results$beta<-or.beta
  results$mu<-or.mu
  results$sigma2<-or.sigma2
  results$theta<-or.theta
  return(results)
}

# Prefictive function for DPM of Conditional Vines

DPMofCondVineCopulaPredictive<-function(sim,pred.sample){
d<-sim$d
dd <- d*(d-1)/2
theta<-sim$theta
mu<-sim$mu
sigma2<-sim$sigma2
beta<-sim$beta
k<-sim$Psi.mode
family_x<-sim$family_x
BinX<-which(family_x=="Binomial")
NormX<-which(family_x=="Normal")
#Estimating Weight in post-process
if(k>1){ 
library(extraDistr)
w<-NULL
Psi<-sim$Psi
for (i in sim$take){
  psi<-Psi[,i]
  if(k==2){
    w_i<-rbeta(1,shape1 = 1+length(which(psi==1)),shape2=length(psi)-length(which(psi==1)+1))
    w<-rbind(w,c(w_i,1-w_i))
  }
  if(k>2) w<-rbind(w,rdirichlet(1,1+as.vector(table(psi)))) 
}
w.mean<-apply(w, 2, mean)
}
family <-rep(1,dd)
order <- 1:d
U.sim<-matrix(NA, nrow = pred.sample, ncol=d)
for (i in 1:pred.sample){
 it<-sample(x=sim$take,size=1)
 if(k==1) component=1
 if(k>1) component<-sample(x=1:k,size=1,prob=w.mean)
 par1<-NULL
  for(p in 1:dd){
    x<-c()
    x[1]<-1
    if(length(NormX)>0) x[NormX+1]<-rnorm(length(NormX),mu[it,component,],sigma2[it,component,])
    if(length(BinX)>0)  x[BinX+1]<-rbinom(n=length(BinX),size=1,theta[it,component,]) 
    if(k==1) par1 =c(par1,beta[,p,component,it]%*%x)
    if(k>1)  par1 =c(par1,beta[,p,component,it]%*%x)
  }
  par2 <- rep(0,dd)
  RVM <- D2RVine(order, family, fish.inv(par1), par2)
  U.sim[i,] <-RVineSim(1, RVM, U = NULL)
}
return(U.sim)
}

# Prefictive function for DPM of Vines

DPMofVineCopulaPredictive<-function(sim,pred.sample){
  d<-sim$d
  dd <- d*(d-1)/2
  beta<-sim$beta
  k<-sim$Psi.mode
  #Estimating Weight in post-process
  if(k>1){ 
    library(extraDistr)
    w<-NULL
    Psi<-sim$Psi
    for (i in sim$take){
      psi<-Psi[,i]
      if(k==2){
        w_i<-rbeta(1,shape1 = 1+length(which(psi==1)),shape2=length(psi)-length(which(psi==1)+1))
        w<-rbind(w,c(w_i,1-w_i))
      }
      if(k>2) w<-rbind(w,rdirichlet(1,1+as.vector(table(psi)))) 
    }
    w.mean<-apply(w, 2, mean)
  }
  family <-rep(1,dd)
  order <- 1:d
  U.sim<-matrix(NA, nrow = pred.sample, ncol=d)
  for (i in 1:pred.sample){
    it<-sample(x=1:dim(beta)[3],size=1) # sample(x=sim$take,size=1) #
    if(k==1) component=1
    if(k>1) component<-sample(x=1:k,size=1,prob=w.mean)
    par1<-NULL
    for(p in 1:dd){
      if(k==1) par1 =c(par1,beta[p,component,it])
      if(k>1)  par1 =c(par1,beta[p,component,it])
    }
    par2 <- rep(0,dd)
    RVM <- D2RVine(order, family, fish.inv(par1), par2)
    U.sim[i,] <-RVineSim(1, RVM, U = NULL)
  }
  return(U.sim)
}





# Function for evaluating the densities of the covariates (assumed to be independent)
dCovariates<-function(X,family_x,mu,sigma2,theta){
  BinX<-which(family_x=="Binomial")
  NormX<-which(family_x=="Normal")
  dCov<-NULL
  if(length(BinX)>0){
    if(dim(X)[1]>1) dCov<-cbind(dCov,rowSums(dbinom(x=X[,BinX,],size=1,prob=theta,log = TRUE)) )
    if(dim(X)[1]==1) dCov<-cbind(dCov,sum(dbinom(x=X[,BinX,],size=1,prob=theta,log = TRUE)) )
  }
  if(length(NormX)>0){
    if(dim(X)[1]>1) dCov<-cbind(dCov,rowSums(dnorm(X[,NormX,],mu,sd=sqrt(sigma2),log = TRUE)) )
    if(dim(X)[1]==1) dCov<-cbind(dCov,sum(dnorm(X[,NormX,],mu,sd=sqrt(sigma2),log = TRUE)) )
  }
  if(length(family_x)>1) dCov<-rowSums(dCov) 
  return(dCov)
}






