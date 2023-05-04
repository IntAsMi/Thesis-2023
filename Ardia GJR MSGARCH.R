rm(list = ls())
load("AllForecasts.Rda")
load("data.Rda")
load("pricesdata.Rda")
load("tau.Rda")
# load("Forecasts of Value of the Portfolio k=4 k=11 k=36/Vhat36.Rda")
require('bayesGARCH'); require('invgamma'); require('EnvStats');
library(MSGARCH);library(freqdom);library(lubridate);library(truncnorm);library(truncdist)
nstep<-135;tsdata<-ts(data[,-1]);pricesdata <- as.data.frame(pricesdata)
datarows<-nrow(tsdata);Z<-ts(pricesdata[,-1]) ;Z_t<-Z[nrow(Z),]

for (ndpc in c(36)) {
  
  if (!file.exists(paste0("dpca_",ndpc,".RDS"))) {
    print(noquote('dpca file not found. Calculating DPCAs...'))
    dpca<-dpca(tsdata, q=10, Ndpc=ndpc) # of principal components, and a kernel estimator of 10 lags
    saveRDS(dpca,paste0('dpca_',ndpc,'.RDS'))
  }else {
    dpca <- readRDS(paste0("dpca_",ndpc,".RDS"))
  }
  chainsize<-100000;rem<-5000;nchains <- 2
  burnin<-chainsize-rem;k <- 1;
  
  
  #creating empty arrays
  ScoresN<-array(NA,dim=c(chainsize, nstep, ndpc))
  h<-array(NA,dim=c(chainsize, ndpc+nstep, ndpc))
  ht<-array(NA,dim=c(nstep, 1, ndpc))
  smpl.full<-array(NA,
                   dim=c(chainsize, 5, ndpc))
  
  spec <- CreateSpec(
    variance.spec = list(model = c("gjrGARCH"),garchOrder = c(1,1)),
    distribution.spec = list(distribution = c('std')),
    switch.spec = list(do.mix = F, K = 1),
    prior = list(mean = list(alpha0_1 = 0.001, alpha1_1 = 0.05,
                             alpha2_1 = 0.1, beta_1 = 0.8, nu_1 = 2.1)))
  
  # Extract the parameter values and bounds
  proposal <- function(theta, tuning) {
    alpha0 <- theta$alpha0
    alpha1 <- theta$alpha1
    alpha2 <- theta$alpha2
    beta <- theta$beta
    nu <- theta$nu
    
    alpha0.lower <- spec$prior$variance$lower[["alpha0"]]
    alpha0.upper <- spec$prior$variance$upper[["alpha0"]]
    
    alpha1.lower <- spec$prior$variance$lower[["alpha1"]]
    alpha1.upper <- spec$prior$variance$upper[["alpha1"]]
    
    alpha2.lower <- spec$prior$variance$lower[["alpha2"]]
    alpha2.upper <- spec$prior$variance$upper[["alpha2"]]
    
    beta.lower <- spec$prior$variance$lower[["beta"]]
    beta.upper <- spec$prior$variance$upper[["beta"]]
    
    nu.lower <- spec$prior$distribution$lower[["nu"]]
    nu.upper <- spec$prior$distribution$upper[["nu"]]
    
    # Generate a proposal value for each parameter using rtruncnorm from the truncnorm package
    alpha0.prop <- rtruncnorm(1, a = alpha0.lower, b = alpha0.upper, mean = alpha0, sd = tuning[["alpha0"]])
    alpha1.prop <- rtruncnorm(1, a = alpha1.lower, b = alpha1.upper, mean = alpha1, sd = tuning[["alpha1"]])
    alpha2.prop <- rtruncnorm(1, a = alpha2.lower, b = alpha2.upper, mean = alpha2, sd = tuning[["alpha2"]])
    beta.prop <- rtruncnorm(1, a = beta.lower, b = beta.upper, mean = beta, sd = tuning[["beta"]])
    nu.prop <- rtruncexp(1, a = nu.lower, b = nu.upper, rate = tuning[["nu"]], translate = 2)
    
    # Return a list of proposal values
    list(alpha0 = alpha0.prop,
         alpha1 = alpha1.prop,
         alpha2 = alpha2.prop,
         beta = beta.prop,
         nu = nu.prop)
  }
  
  for(k in 1:ndpc){
    #kth principal component
    print(Sys.time())
    y<-dpca$scores[,k]-median(dpca$scores[,k])
    #RUN THE SAMPLER
    addPriorConditions<-function(psi){psi[2]+psi[3] +psi[4]<1}#imposing the covariance-stationary conditions to hold, i.e. a +b < 1 
    MCMC<- MSGARCH::FitMCMC(spec,y,
                            ctr = list(nmcmc = chainsize,nburn = burnin,
                                       nthin = 1,do.se = FALSE,
                                       addPriorConditions=addPriorConditions),
                            proposal.fn=proposal,
                            tuning=list(alpha0=0.01,
                                        alpha1=0.01,
                                        alpha2=0.01,
                                        beta=0.01,
                                        nu=2.01))
    print(Sys.time())
    #form the posterior sample allowing a burn in of its size
    smpl.full[,,k] <- as.matrix(MCMC$par)
    pred <- predict(MCMC, newdata = NULL, nahead = nstep, do.return.draw = T,
                    do.cumulative = FALSE, ctr = list())
    ht[,,k]  <-pred$vol
    ScoresN[,,k] <- pred$draw + median(dpca$scores[,k])
    print(paste0('fitMCMC completed for PC ',k))
  }
  
  dimnames(smpl.full) <- list(NULL,colnames(MCMC$par),NULL)
  saveRDS(smpl.full, paste0('smplfull_GJR',ndpc,'_',nstep,'.RDS'))
  saveRDS(ScoresN,paste0('ScoresNGJR',ndpc,'_',nstep,'.RDS'))
  saveRDS(ht,paste0('htGJR',ndpc,'_',nstep,'.RDS'))
  
  if(file.exists(paste0('YYGJR',ndpc,'_',nstep,'.RDS'))) {
    YY <- readRDS(paste0('YYGJR',ndpc,'_',nstep,'.RDS'))
  }else{
    forecasts<-array(rep(NA, nstep*ndpc*rem),
                     dim=c(nstep,ndpc,rem))
    #Y<-array(rep(NA, (nrow(tsdata)+nstep)*(ndpc)*(rem)),
    # dim=c((nrow(tsdata)+nstep),(ndpc),(rem)))
    xhat<-array(NA, dim=c(nrow(tsdata)+nstep,ncol(tsdata),rem)) 
    xhatn<-array(NA, dim=c(rem, ncol(tsdata), nstep))
    Y<-array(NA, dim=c((nrow(tsdata)+nstep),(ndpc),(rem)))
    colnames(xhat)<-colnames(tsdata); colnames(xhatn)<-colnames(tsdata)
    for(g in 1:rem){ 
      for(i in 1:nstep){
        for(k in 1:ndpc){
          forecasts[i,k,g]<-ScoresN[g,i,k] }
      } }
    
    YY <- array(unlist(lapply(1:rem, function(x) rbind(dpca$scores[,1:ndpc], forecasts[,,x]))),
                dim = c((nrow(tsdata)+nstep),(ndpc),(rem)))
    saveRDS(YY,paste0('YYGJR',ndpc,'_',nstep,'.RDS'))
  }
  
  
  
  if (file.exists(paste0('xhat_GJR',ndpc,'_',nstep,'.RDS'))) {
    xhat <- readRDS(paste0('xhat_GJR',ndpc,'_',nstep,'.RDS'))	
  }else{
    for(g in 1:rem) {
      xhat[,,g] <- YY[,,g]%c%freqdom.transpose(rev(dpca$filters))
    }
    saveRDS(xhat,file = paste0('xhat_GJR',ndpc,'_',nstep,'.RDS'))
  }
  
  if (file.exists(paste0('xhatn_GJR',ndpc,'_',nstep,'.RDS'))) {
    xhatn <- readRDS(paste0('xhatn_GJR',ndpc,'_',nstep,'.RDS'))
  }else{
    for (i in 1:nstep) {
      for (j in 1:rem) {
        xhatn[j,,i] <- t(xhat[datarows+i,,j])
      }
      # print(i)
    }
    saveRDS(xhatn,file = paste0('xhatn_GJR',ndpc,'_',nstep,'.RDS'))
  }
  
  if (file.exists(paste0('zhatn_GJR',ndpc,'_',nstep,'.RDS'))) {
    zhatn <- readRDS(paste0('zhatn_GJR',ndpc,'_',nstep,'.RDS'))
  }else{
    #transforming into prices
    zhatn<-array(rep(NA, nstep*rem*ncol(tsdata)), dim=c(rem, ncol(tsdata), nstep))
    for(i in 1:rem){
      for(j in 1:ncol(tsdata)){
        for(k in 1:nstep){ if(k==1){
          zhatn[i,j,k]<-exp(xhatn[i,j,k])*(Z_t[j])
        }else{
          zhatn[i,j,k]<-exp(xhatn[i,j,k])*zhatn[i,j,k-1] }
        }
      }
      # print(noquote(paste0(round(i/rem*100,3),'%')))
    }
    saveRDS(zhatn,file = paste0('zhatn_GJR',ndpc,'_',nstep,'.RDS'))
  }
}

