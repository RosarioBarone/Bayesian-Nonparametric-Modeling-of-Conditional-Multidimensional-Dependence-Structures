####################################################################################
########################### Code by: Rosario Barone   ##############################
####################################################################################

#setwd("")
library(VineCopula)
source("Utils.R")
source("DPMM of Cond Vines.R")



sim.cond.vine.mixture<-function(n,d,true.beta,true.mu.x,true.sigma.x,true.theta.x,true.w,famCop){
  require(CDVineCopulaConditional)
  set.seed(8)
  dd<-d*(d-1)/2
  true.k<-length(true.w)
  U<-matrix(NA, nrow = n, ncol=d) # Matrix of pseudo-observations
  X<-array(data=NA,dim = c(n,ncov,dd)) # Array of covariates
  n.k<-round(n*true.w)
  k_set<-cumsum(c(0,n.k))
  order <- 1:d
  for (j in 1:true.k){
    for(i in (k_set[j]+1):k_set[j+1]){
      par1<-NULL
      for(p in 1:dd){
        x<-NULL
        if(length(which(family_x=="Normal"))>0){
          for (nx in 1:length(which(family_x=="Normal"))) {
            x<-c(x,rnorm(1,true.mu.x[,j],true.sigma.x[,j]))
          }  
        } 
        if(length(which(family_x=="Binomial"))>0){
          for (nx in 1:length(which(family_x=="Binomial"))) {
            x<-c(x,rbinom(1,size=1,prob=true.theta.x[,j]))
          }  
        }
        x<-c(1,x)
        par1=c(par1,true.beta[j,p,]%*%x)
        X[i,,p]<-x[-1]
      }
      par2 <- rep(0,dd)
      RVM <- D2RVine(order, family=famCop, fish.inv(par1), par2)
      U[i,] <-RVineSim(1, RVM, U = NULL)
    }
  }
  sim<-c()
  sim$U<-U
  sim$X<-X
  return(sim)
}



d<-4
dd<-d*(d-1)/2
ncov<-1
true.w<-c(0.5,0.5)  #Weight of each mixture component
true.k<-length(true.w) # Number of mixture components
family_x<-c("Normal")
#True values of the conditional copula parameters (Intercept + Coeff) for the two mixture comonents
true.beta<-array(NA,dim = c(true.k,dd,1+ncov))
true.beta[1,,]<-matrix(c(-.3,-0.5,
                         -.3,-0.5,
                         -.3,-0.5,
                         -.3,-0.5,
                         -.3,-0.5,
                         -.3,-0.5),ncol=1+ncov, byrow=TRUE)
true.beta[2,,]<-matrix(c(0.3,0.5,
                         0.3,0.5,
                         0.3,0.5,
                         0.3,0.5,
                         0.3,0.5,
                         0.3,0.5),ncol=1+ncov, byrow=TRUE)
true.mu.x<-matrix(c(3,3),nrow=length(which(family_x=="Normal"))) #true value of the mean for the Normal covariates for each mixture component
true.sigma.x<-matrix(c(1,1), nrow = length(which(family_x=="Normal"))) #true value of the sd for the Normal covariates for each mixture component
true.theta<-NULL #true value of the parameter of the binary covariates for each mixture component
famCop<-rep(1,dd) #Define the copula family from "CDVineCopulaConditional"



data<-sim.cond.vine.mixture(n=50,d,true.beta,true.mu.x,true.sigma.x,true.theta.x,true.w,famCop)
U<-data$U
X<-data$X
plot(U[,1:2])

##Define List with Prior parameters for the covariates posterior parameters
#In this setting we only have 1 covariate assumed to be Normal 
HyperCov<-list()
HyperCov$mu0_x<-0   
HyperCov$s0_x<-1000


sim<-DPMofConditionalVines(#dimension of the vine
                           U, #Observations in [0,1]
                           X, #Covariates
                           family_x, #Distribution assumption on the covariates
                           ncov, #Number of covariates
                           M=1/20, #Concentration parameter
                           m0=0, # {base measure mean}
                           s0=1, # {base measure sd}
                           HyperCov, #list of the prior parameters for the covariates
                           JointSampling=FALSE, #Sampling strategy: if true, vine parameters are simualted
                           # jointly for each mixture component.
                           s.prop=0.02, #variance of the proposal distributions for the MH step
                           MCMCiter=1000, #Number of MCMC iterations
                           burnin=100 #Burnin
                           )


predictive<-DPMofCondVineCopulaPredictive(sim,pred.sample = 500)
par(mfrow=c(1,1))
plot(predictive[,1:2])
plot(predictive[,2:3])
plot(predictive[,3:4])

#If there is Label switching, we also provide a relabelling function:
#relabelling<-relabel.condVineDPM(U,X,ncov,family_x=c("Normal"),sim)

