####################################################################################
########################### Code by: Rosario Barone   ##############################
#########    Dirichlet Process Mixture of Conditional Vine Copulas  ################
####################################################################################



DPMofConditionalVines<-function(U,X,family_x,ncov,M,m0,s0,HyperCov,JointSampling,s.prop,MCMCiter,burnin){
  require(VineCopula)
  require(mvtnorm)
  require(latex2exp)
  require(combinat)
  n<-dim(U)[1]
  d<-dim(U)[2]
  dd=d*(d-1)/2
  ncov=length(family_x)
  XBin<-sum(family_x=="Binomial")
  XNorm<-sum(family_x=="Normal")
  X.obs<-X
  beta<-array(0,dim=c(1+ncov,d*(d-1)/2,n,MCMCiter))
  beta[,,,1]=runif(length(beta[,,,1]),-0.005,0.005)
  psi<-rep(1,n)
  Psi<-matrix(0,n,MCMCiter)
  Psi[,1]<-psi
  if(isTRUE(JointSampling)) s.draw<-diag(s.prop,nrow=d*(d-1), ncol = d*(d-1)) else s.draw<-diag(s.prop,nrow=1+ncov, ncol =1+ncov)
  if(XNorm>0){
    mu<-array(0,dim=c(MCMCiter,n,XNorm)) ; sigma2<-array(0,dim=c(MCMCiter,n,length(which(family_x=="Normal")))); mu[1,,]=0.5;sigma2[1,,]=0.5;
  }
  if(XBin>0){ 
    theta<-array(0,dim=c(MCMCiter,n,XBin)); theta[1,,]=0.5;
  }
  
  
  ## Function fo updating the posterior distribution of the Normal covariates
  GibbsUpdateNormLik<-function(y){
    mu0_x<-HyperCov$mu0_x 
    s0_x<-HyperCov$s0_x
    n<-length(y)
    mu<-c()
    sigma2<-c()
    k0= mu0_x
    kn=k0+n
    mun<-n*mean(y)/kn
    nun<-n+1
    s2n<- (s0_x + n*var(y) + k0*n*(mean(y))^2/kn)/nun
    sigma2<-1/rgamma(1,shape=nun/2,rate=s2n*nun/2)
    mu<-rnorm(1,mean=mun,sd=sqrt(sigma2/n))
    out<-list()
    out$mu<-mu
    out$sigma2<-sigma2
    return(out)
  }
  
  ## Function fo updating the posterior distribution of the binary covariates
  
  GibbsUpdateBinomLik<-function(y){
    a<-HyperCov$a
    b<-HyperCov$b  
    theta<-rbeta(1,a+sum(y==1),b+sum(y==0))
    out<-list()
    out$theta<-theta
    return(out)
  }
  
  for (iter in 2:MCMCiter){
    
    if(iter%%100==0){
      print(paste("Iteration n",iter))
      par(mar=c(1,1,1,1))
      par(mfrow=c(3,1))
      Psi.max<-apply(Psi[,1:(iter-1)], 2, max)
      barplot(table(Psi.max), main="Observed Clusters")
      barplot(table(Psi[,(iter-99):(iter-1)]), main="Observations Across Clusters")
      plot(Psi.max)
    }
    
    psi.star = max(psi) + 1
    if(psi.star>n) psi.star=n
    for (i in 2:n){
      p_i<-NULL
      psi.minus<-psi[-i]
      psi.star<-length(unique(psi.minus)) + 1
      for(j in unique(psi.minus)){
        p_i<-c(p_i,log(length(which(psi.minus==j)))+dVineCopCond(U[i,],beta[,,j,iter-1],x=array(X[i,,],dim = c(1,ncov,dd)),d=d,log = TRUE) + dCovariates(X=array(X[i,,],dim = c(1,ncov,dd)),family_x=family_x,mu=mu[iter-1,j,],sigma2=sigma2[iter-1,j,],theta=theta[iter-1,j,]))
      }
      p_i<-c(p_i,log(M) -log(psi.star) + dCovariates(X=array(X[i,,],dim = c(1,ncov,dd)),family_x=family_x,mu=mu[iter-1,psi.star,],sigma2=sigma2[iter-1,psi.star,],theta=theta[iter-1,psi.star,]) +
               dVineCopCond(U[i,],beta[,,psi.star,iter-1],x=array(X[i,,],dim = c(1,ncov,dd)),d=d,log = TRUE) )
      psi[i]<-sample(1:length(p_i),size=1,prob = exp(p_i)/sum(exp(p_i)))
    }
    psi.new<-psi
    psi<-c()
    un<-unique(psi.new)
    for(i in 1:length(un)){
      psi[which(psi.new==un[i])]=i
    }
    psi.new<-psi
    #Relabel posterior values in order to fill the gap of empty clusters
    beta[,,1:length(un),iter-1]<-beta[,,un,iter-1]
    if(XNorm>0){
      mu[iter-1,1:length(un),]<-mu[iter-1,un,] ; sigma2[iter-1,1:length(un),]<-sigma2[iter-1,un,]
    }
    if(XBin>0) theta[iter-1,1:length(un),]<-theta[iter-1,un,] 
    ##Check if there is a cluster of size 1
    cluster.size<-NULL
    phi<-sort(unique(psi.new))
    for(tau in phi){
      cluster.size<-c(cluster.size,length(which(psi.new==phi[tau])))
    }
    single.cluster<-sum(which(cluster.size==1))
    
    if(single.cluster!=0){
      unchanged<-sample(x=c(1,NA), size = length(single.cluster), prob=c((max(psi)-1)/max(psi),1-(max(psi)-1)/max(psi) ))
      single.cluster<-single.cluster*unchanged
      single.cluster<-na.omit(single.cluster)
      if(sum(single.cluster)==0) single.cluster<-0
    }
    if(single.cluster!=0){
      reallocate<-phi[which(cluster.size==1)]
      single.obs<-which(psi.new%in%reallocate)
      clusters_set<-sort(phi[-reallocate])
      star<-cbind(rep(1:n),psi.new)
      
      if(length(clusters_set)==1){
        psi.new[single.obs]=clusters_set
      }else{
        
        for(i in single.obs){
          psi.minus<-psi.new[-i]
          p_i<-NULL
        for(j in clusters_set) p_i<-cbind(p_i,log(length(which(psi.minus==j)))+dVineCopCond(U[i,],beta[,,j,iter-1],array(X[i,,],dim = c(1,ncov,dd)),d=d,log = TRUE) + dCovariates(X=array(X[i,,],dim = c(1,ncov,dd)),family_x=family_x,mu=mu[iter-1,j,],sigma2=sigma2[iter-1,j,],theta=theta[iter-1,j,]))
        psi.new[i]=sample(clusters_set,size = 1 ,prob = exp(p_i)/sum(exp(p_i)))
        }
        }
      
    }
    
    #### Set clusters in order of appearence ###
    psi<-c()
    un<-unique(psi.new)
    for(i in 1:length(un)){
      psi[which(psi.new==un[i])]=i
    }
    psi.new<-psi
    #Relabel posterior values in order to fill the gap of empty clusters
    beta[,,1:length(un),iter-1]<-beta[,,un,iter-1]
    if(XNorm>0){
      mu[iter-1,1:length(un),]<-mu[iter-1,un,] ; sigma2[iter-1,1:length(un),]<-sigma2[iter-1,un,]
    }
    if(XBin>0) theta[iter-1,1:length(un),]<-theta[iter-1,un,] 
    
    Psi[,iter] = psi
    
    ##### Draw from the posterior distribution ####
    
    for (k in 1:n){
      if (length(which(psi==k))==0){
        beta[,,k,iter]<-t(rmvnorm(d*(d-1)/2,mean=rep(m0,ncov+1),sigma = diag(s0,ncov+1)))
        if(XNorm>0){
          mu[iter,k,]<-rnorm(XNorm,0,s0) ;sigma2[iter,k,]<-1/rgamma(XNorm,0.1,0.1)
        }
        if(XBin>0) theta[iter,k,]<-rbeta(XBin,1,1)
      }else{
        n.star<-which(psi==k)
        if(XNorm>0){
          if(XNorm>1){
            NormCov<-apply(X[,which(family_x=="Normal"),],2,GibbsUpdateNormLik)
            norm.par<-matrix(unlist(NormCov),nrow=XNorm,ncol=2)
            mu[iter,k,]<-norm.par[,1]
            sigma2[iter,k,]<-norm.par[,2]
          }else{
            NormCov<-GibbsUpdateNormLik(X[,which(family_x=="Normal"),])
            mu[iter,k,]<-NormCov$mu
            sigma2[iter,k,]<-NormCov$sigma2  
          }
        }
        if(XBin>0){
          if(XBin>1){
            theta[iter,k,]<-unlist(apply(X[,which(family_x=="Binomial"),], 2, GibbsUpdateBinomLik))
          }else{
            theta[iter,k,]<-unlist(GibbsUpdateBinomLik(X[,which(family_x=="Binomial"),]))
          }
        }
        if(isTRUE(JointSampling)){ 
          last.b<-beta[,,k,iter-1]
          prop.b<-matrix(rmvnorm(1,mean=c(beta[,,k,iter-1]),sigma = s.draw),nrow=1+ncov,ncol =d*(d-1)/2)
          l.lik.num<-sum(dVineCopCond(U[n.star,],prop.b,array(X[n.star,,],dim = c(length(n.star),ncov,dd)),d=d,log=T))
          l.lil.den<-sum(dVineCopCond(U[n.star,],last.b,array(X[n.star,,],dim = c(length(n.star),ncov,dd)),d=d,log=T))
          ACCEPT<-exp(l.lik.num+sum(dnorm(prop.b,m0,s0, log = TRUE))-l.lil.den-sum(dnorm(last.b,m0,s0, log = TRUE))); if(is.na(ACCEPT)) ACCEPT=0 #Fix numerical problems
          if(ACCEPT>runif(1)){beta[,,k,iter]<-prop.b}else{beta[,,k,iter]<-last.b} 
        }else{
          last.b<-beta[,,k,iter-1]  
          prop.b<-beta[,,k,iter-1]
          for (p in 1:(d*(d-1)/2)){  
            prop.b[,p]<-rmvnorm(1,mean = last.b[,p],sigma = s.draw)
            l.lik.num<-sum(dVineCopCond(U[n.star,],prop.b,array(X[n.star,,],dim = c(length(n.star),ncov,dd)),d=d,log=T))
            l.lil.den<-sum(dVineCopCond(U[n.star,],last.b,array(X[n.star,,],dim = c(length(n.star),ncov,dd)),d=d,log=T))
             ACCEPT<-exp(l.lik.num+sum(dnorm(prop.b[,p],m0,s0, log = TRUE))-l.lil.den-sum(dnorm(last.b[,p],m0,s0, log = TRUE))); if(is.na(ACCEPT)) ACCEPT=0 #Fix numerical problems
            if(ACCEPT>runif(1)){ beta[,p,k,iter]<-prop.b[,p] ; last.b[,p]<-prop.b[,p] }else{ beta[,,k,iter]<-last.b[,p]; prop.b[,p]<- last.b[,p]   }
          }
        }
      }
    }
  }
  Psi.max<-apply(Psi[,1:iter], 2, max)
  Psi.mode<-which(table(Psi.max[1:iter])==max(table(Psi.max[1:iter])))
  take<-which(Psi.max==Psi.mode)[-(1:burnin)]
  results<-list()
  results$d<-d
  results$take<-take
  results$Psi<-Psi
  results$beta<-beta
  results$Psi.mode<-Psi.mode
  results$Psi.max<-Psi.max
  if(XNorm>0){ 
    results$mu<-mu; results$sigma2<-sigma2 
    results$post.mean.mu<-apply(mu[take,,],1:2,mean)[1:Psi.mode]
    results$post.sd.mu<-apply(mu[take,,],1:2,sd)[1:Psi.mode]
    results$post.mean.sigma2<-apply(sigma2[take,,],1:2,mean)[1:Psi.mode]
    results$post.sd.sigma2<-apply(sigma2[take,,],1:2,sd)[1:Psi.mode]
  }
  if(XBin>0){
    results$theta<-theta
    results$post.mean.theta<-apply(theta[take,,],1:2,mean)[1:Psi.mode]
    results$post.sd.theta<-apply(theta[take,,],1:2,sd)[1:Psi.mode]
  }
  results$post.mean.beta<-apply(beta[,,,take], 1:3, mean)[,,1:Psi.mode]
  results$post.sd.beta<-apply(beta[,,,take], 1:3, sd)[,,1:Psi.mode]
  results$family_x<-family_x
  return(results)
}

