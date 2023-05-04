rm(list = ls())
load("Data Set/yahfin/AllForecasts.Rda")
load("Data Set/yahfin/data.Rda")
load("Data Set/yahfin/pricesdata.Rda")
load("Data Set/yahfin/tau.Rda")
load("Forecasts of Value of the Portfolio k=4 k=11 k=36/Vhat36.Rda")
setwd()
require('EnvStats');library(doParallel);library(stats);
library(invgamma);library(lubridate);library(coda);library(R2jags)
library(freqdom);library(lubridate);library(rjags);
ndpc<-4;nstep<-135;tsdata<-ts(data[,-1]);pricesdata <- as.data.frame(pricesdata)
datarows<-nrow(tsdata);Z<-ts(pricesdata[,-1]) ;Z_t<-Z[nrow(Z),]

asyfn <- 'exp'

for (ndpc in c(4,11,36)) {
	#applying dynamic pca
	# ndpc<-36;
	if (!file.exists(paste0("dpca_",ndpc,".RDS"))) {
		dpca<-dpca(tsdata, q=10,Ndpc = ndpc) 
		saveRDS(dpca,paste0('dpca_',ndpc,'.RDS'))
	}else {
		dpca <- readRDS(paste0("dpca_",ndpc,".RDS"))
	}
	
	#predicting the one-step-ahead score for each PC
	chainsize<-1e+4;chains = 2
	burnin<-chainsize*1/2;k <- 1;
	rem<-5e+3#-burnin
	eps<-0.0001;K<-1
	
	#creating empty arrays
	k <- 1

	ScoresN<-array(NA,dim=c(rem, nstep, ndpc), dimnames = list(NULL,paste0('y.h',1:nstep),1:ndpc))
	h<-array(NA,dim=c(rem, datarows+nstep, ndpc))
	smplfull<-array(NA,dim=c(rem, 8, ndpc),dimnames = list(NULL,c("mu0", "sigma0", "arch", "garch", "v", "gamma1","gamma2", "loglik"),NULL))
	SMPL<-array(NA,dim=c(chainsize + rem, 8, ndpc),	dimnames = list(NULL,c("mu0", "sigma0", "arch", "garch", "v", "gamma1","gamma2", "loglik"),NULL))
	
	f_exp = function(u_t_1, gamma) {return(1 - exp(-gamma * u_t_1^2))	}
	f_logit = function(u_t_1, gamma) {return(1/(1 + exp(-gamma * u_t_1)))}
	library(parallel)
	options(mc.cores=chains)
	
	# generate a random value from a uniform distribution
	init <- function(chain_id) {c(mu0 = runif(1), sigma0 = runif(1),arch = runif(1),garch = runif(1),gamma = runif(1),df = 2 + runif(1,0,10))}
	
	for(k in 1:(ndpc)){
		print(paste0('start: ',k))
		#kth principal component
		y<-dpca$scores[,k]-median(dpca$scores[,k]);n<-length(y)
		
		# Specify the priors
		prior_mu0 <- uniform(eps,K) 
		prior_sigma0 <- uniform(eps,K) 
		prior_arch <- uniform(eps,K) # prior for the ARCH term in GARCH eq
		prior_garch <- uniform(eps,K) # prior for the GARCH term in GARCH eq
		prior_gamma <- uniform(eps,K) # prior for the asy parameter in GARCH eq
		prior_df <- jeffrey() # use an independent Jeffrey's prior for dof par
		
		print(a <- Sys.time())
		# Fit the model using Stan
		MCMC <- stan_garch(y, order = c(1,1,0), arma = c(0,0), genT = TRUE,
											 asym = asyfn,
											 iter = chainsize + rem,
											 warmup = floor((chainsize + rem)/5),
											 chains = chains,stepwise = T,
											prior_mu0 = prior_mu0,
											prior_sigma0 = prior_sigma0,
											prior_arch = prior_arch,
											prior_garch = prior_garch,
											prior_gamma = prior_gamma,
											prior_df = prior_df,
											verbose = TRUE, sum_constraint = TRUE,
											compute_criteria = F, init = init)
		print(Sys.time()-a)
		
		print(k)
		#POSTERIOR STATISTICS 
		sample <- extract_stan(MCMC,
													 pars = c(get_parameters(MCMC)), 
													 permuted = F,#extracted from each chain
													 inc_warmup = F)#not including burn in
		smplfull[,,k] <- sapply(1:length(dimnames(smplfull)[[2]]), function(y)
			sapply((dim(sample)[1]-rem+1):(dim(sample)[1]),
						 function(x) mean(sample[x,,y])))
		
		sample <- extract_stan(MCMC,
													 pars = c(get_parameters(MCMC)), 
													 permuted = F,#extracted from each chain
													 inc_warmup = T)#not including burn in
		SMPL[,,k] <- sapply(1:length(dimnames(smplfull)[[2]]), function(y)
			#taking the mean of each chain and par
			sapply(1:(dim(sample)[1]), function(x) mean(sample[x,,y])))
		
		ScoresN[,,k]<- as.matrix(posterior_predict(MCMC,draws = rem,h=nstep,
																							 seed = 123)) + 
			median(dpca$scores[,k])
		
		MCMC <- {}
		
		saveRDS(smplfull,paste0('smplfull_FSC',asyfn,ndpc,'_',nstep,'.RDS'))
		saveRDS(SMPL,paste0('SMPL_FSC',asyfn,ndpc,'_',nstep,'.RDS'))
		saveRDS(ScoresN,paste0('ScoresNFSC',asyfn,ndpc,'_',nstep,'.RDS'))
		
		for(i in 1:rem){ h0<-(smplfull[i,'sigma0',k])/(
			1-smplfull[i,'arch',k]-smplfull[i,'garch',k]) 
		for(j in 1:(n+1)){
			#if j==1 we use h0 else we use previous h
			if(j==1){ h[i,j,k]<-smplfull[i,'sigma0',k] + (smplfull[i,'garch',k]*h0)
			} else{ 
				h[i,j,k]<-smplfull[i,'sigma0',k] + smplfull[i,'arch',k] * y[j-1]^2*  + smplfull[i,'gamma1',k] * y[j-1]^2 * f_exp(y[j-1],smplfull[i,'gamma2',k]) +smplfull[i,'garch',k] * h[i,j-1,k]
			} }
		}
		
		for(i in 1:rem){
			for(j in 1:nstep){ 
				if(j!=1){
					h[i,(datarows+j),k]<-smplfull[i,'sigma0',k] + smplfull[i,'arch',k] *ScoresN[i,(j- 1),k]^2*  + smplfull[i,'gamma1',k] *ScoresN[i,(j- 1),k]^2 * f_exp(y[j-1],	smplfull[i,'gamma2',k]) +smplfull[i,'garch',k] * h[i,(datarows+(j-1)),k]
					if(is.na(h[i,(datarows+j),k])) stop()
				}
			}
		}
		
		saveRDS(h,paste0('h_FSC',asyfn,ndpc,'_',nstep,'.RDS'))
		
		gc(full = T)
		print(k)
	}
	
	forecasts<-array(NA, dim=c(nstep,ndpc,rem))
	xhat<-array(NA, dim=c(nrow(tsdata)+nstep,ncol(tsdata),rem),
							list(NULL,colnames(tsdata),NULL)) 
	xhatn<-array(NA, dim=c(rem, ncol(tsdata), nstep),
							 list(NULL,colnames(tsdata),NULL))
	Y<-array(NA, dim=c((nrow(tsdata)+nstep),(ndpc),(rem)))
	zhatn<-array(NA, dim=c(rem, ncol(tsdata), nstep))
	for(g in 1:rem){ 
		for(i in 1:nstep){
			for(k in 1:ndpc){
				forecasts[i,k,g]<-ScoresN[g,i,k] }
		} }
	
	for(g in 1:rem){
		Y[,,g]<-rbind(dpca$scores, forecasts[,,g]);
		xhat[,,g]<-((Y[,,g])%c%freqdom.transpose(rev(dpca$filters)))
		for(i in 1:nstep){ xhatn[g,,i]<-xhat[datarows+i,,g]
		}
	}
	# saveRDS(xhatn,file = paste0('xhatn_FSC',asyfn,ndpc,'_',nstep,'.rds'))
	# saveRDS(Y,paste0('YGJR_FSC',asyfn,ndpc,'_',nstep,'.RDS'))
	# saveRDS(xhat,file = paste0('xhat_FSC',asyfn,ndpc,'_',nstep,'.rds'))
	
	#transforming into prices
	for(i in 1:nstep){
		for(g in 1:rem){
			for(p in 1:ncol(tsdata)){ if(i==1){
				zhatn[g,p,i]<-exp(xhatn[g,p,i])*(Z[1187,p]) 
			}else{
				zhatn[g,p,i]<-(exp(xhatn[g,p,i]))*(zhatn[g,p,(i-1)]) }
			} }
	}
	saveRDS(zhatn,file = paste0('zhatn_FSC',asyfn,ndpc,'_',nstep,'.rds'))
}


