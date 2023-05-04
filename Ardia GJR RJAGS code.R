rm(list = ls())
load("AllForecasts.Rda")
load("data.Rda")
load("pricesdata.Rda")
load("tau.Rda")
setwd()
# install.packages('EnvStats');install.packages('doParallel');install.packages('stats'); install.packages('invgamma')
# install.packages('freqdom');install.packages('lubridate');install.packages('rjags');install.packages('coda')
require('EnvStats');library(doParallel);library(stats); library(invgamma)
library(freqdom);library(lubridate);library(rjags);library(coda);library(R2jags)
nstep<-135;tsdata<-ts(data[,-1]);pricesdata <- as.data.frame(pricesdata)
datarows<-nrow(tsdata);Z<-ts(pricesdata[,-1]) ;Z_t<-Z[nrow(Z),]

for (ndpc in c(4,11,36)) {
	setwd()
	if (!file.exists(paste0("dpca_",ndpc,".RDS"))) {
		print(noquote('dpca file not found. Calculating DPCAs...'))
		dpca<-dpca(tsdata, q=10, Ndpc=ndpc) # of principal components, and a kernel estimator of 10 lags
		saveRDS(dpca,paste0('dpca_',ndpc,'.RDS'))
	}else {
		dpca <- readRDS(paste0("dpca_",ndpc,".RDS"))
	}
	
	# specify the model
	# the way that the invgamma distr is implemented here is from
	# https://cran.r-project.org/web/packages/bridgesampling/vignettes/bridgesampling_example_jags.html
	model_string <- "model {
			h[1] <- 0.01
	    for (i in 2:n) {
					epsilon[i] ~ dnorm(0, 1)
	  			invvarpi[i] ~ dgamma((nu + 2) / 2, (nu + 2) / 2)
					varpi[i] <- 1/invvarpi[i]
					h[i] <- alpha0 + alpha1 * pow(y[i-1],2)* (y[i-1]>=0)  + alpha2 * pow(y[i-1],2) * (y[i-1]<0) + beta * h[i-1]
	        y[i] ~ dnorm(sqrt(((nu + 2) - 2) / (nu + 2) * varpi[i] * h[i]) * epsilon[i], 1)
	    }
	
	    alpha0 ~ dnorm(0.01, 1/10000) T(0,)
	    alpha1 ~ dnorm(0.01, 1/10000) T(0,1)
	    alpha2 ~ dnorm(0.01, 1/10000) T(0,1 - alpha1)
	    beta ~ dnorm(0.01, 1/10000) T(0, 1 - alpha1 - alpha2)
	    nu ~ dexp(2.1)T(2,)
	
			# Predictions
			y_pred[1] <- y[n]
			invvarpi_pred[1] <- invvarpi[n]
			varpi_pred[1] <- varpi[n]
			h_pred[1] <- h[n]
			epsilon_pred[1] <- epsilon[n]
		  for (i in (n+2):(n+136)) {
		    epsilon_pred[i-n] ~ dnorm(0, 1)
		    invvarpi_pred[i-n] ~ dgamma((nu + 2) / 2, (nu + 2) / 2)
				varpi_pred[i-n] <- 1/invvarpi_pred[i-n]
		    h_pred[i-n] <- alpha0 + alpha1 * pow(y_pred[i-n-1],2)* (y_pred[i-n-1]>=0) + alpha2 * pow(y_pred[i-n-1],2) * (y_pred[i-n-1]<0) + beta * h_pred[i-n-1]
				y_pred[i-n] ~ dnorm(sqrt(((nu + 2) - 2) / (nu + 2) * varpi_pred[i-n] * h_pred[i-n]) * epsilon_pred[i-n], 1)
		  }
	}"
	
	
	# write the model to a file
	writeLines(model_string, con = "model.txt")
	# Set up initial values
	inits <- function() {
		alpha1 = runif(1, 0, 1);alpha2 = runif(1, 0, 1-alpha1);
		beta = runif(1, 0, 1-alpha1-alpha2); 
		alpha0 = runif(1, 0, 1); nu = runif(1, 2.01, 10)
		return(list(alpha0 = alpha0,alpha1 = alpha1,alpha2 = alpha2,
								beta = beta,nu = nu))
	}
	
	# Define the parameters to monitor
	params <- c("alpha0", "alpha1","alpha2", "beta","nu",
							"y", "h", "y_pred","h_pred")
	
	chainsize<-10000;rem<-5000;nchains <- 2
	burnin<-chainsize-rem;k <- 1;
	
	#creating empty arrays
	ScoresN<-array(NA,dim=c(rem, nstep, ndpc),
								 dimnames = list(NULL,paste0('y.h',1:nstep),1:ndpc))
	h<-array(NA,dim=c(rem, datarows+nstep, ndpc))
	smplfull<-array(NA,dim=c(rem, 5, ndpc),
									dimnames = list(NULL,params[1:5],1:ndpc))
	SMPL<-array(NA,dim=c(chainsize, 5, ndpc),
							dimnames = list(NULL,params[1:5],1:ndpc)

	for(k in 1:ndpc){
		# #kth principal component
		y <- dpca$scores[,k]-median(dpca$scores[,k]);	n <- length(y)
		#RUN THE SAMPLER
		init_values <- lapply(1:nchains, function(x) inits())
	
		# Set up JAGS model & Burn-in
		model <- jags.model(textConnection(model_string),	data = list(y = y, n = n), 
												n.chains = nchains, inits = init_values,
												n.adapt=0, quiet=FALSE) 
		
		# Draw Samples
		print(a <- Sys.time())
		jags_samples_parallel <-
			jags.samples(model, variable.names = params,	n.iter = chainsize)
		print(b <- Sys.time());	print(b-a)
		
		# extract chains
		smplfull[,,k] <- chains <- sapply(params[1:5], function(x) 
			apply(jags_samples_parallel[[x]][1, (burnin+1):chainsize, 1:nchains], 1, mean))
		SMPL[,,k] <- chains <- sapply(params[1:5], function(x) 
			apply(jags_samples_parallel[[x]][1, 1:chainsize, 1:nchains], 1, mean))
		
		saveRDS(smplfull, paste0('smplfull_GJR',ndpc,'_',nstep,'.RDS'))
		saveRDS(SMPL, paste0('SMPL_GJR',ndpc,'_',nstep,'.RDS'))
		
		h[,,k] <- cbind(sapply(1:datarows, function(x) 
			apply(jags_samples_parallel$h[x, (burnin+1):chainsize, 1:nchains],1,mean)),
			sapply(2:(nstep+1), function(x) 
				apply(jags_samples_parallel$h_pred[x, (burnin+1):chainsize, 1:nchains],1,mean)))
		
		saveRDS(h,paste0('hGJR',ndpc,'_',nstep,'.RDS'))
		
		# plot(h[1,,1], type = 'l',ylim = range(h[,,1]))
		# for (i in (burnin+1):chainsize) lines(h[i,,1])
		
		ScoresN[,,k] <- (sapply(1:nstep, function(x) 
			apply(jags_samples_parallel$y_pred[x, (burnin+1):chainsize, 1:nchains],1,mean)) +
			median(dpca$scores[,k]))/10
		
		saveRDS(ScoresN,paste0('ScoresNGJR',ndpc,'_',nstep,'.RDS'))
		
		# plot(ScoresN[1,,1], type = 'l',ylim = range(ScoresN[,,1]))
		# for (i in 1:rem) lines(ScoresN[i,,1])
		
		print(k)
	}
	

	forecasts<-array(NA, dim=c(nstep,ndpc,rem))
	xhat<-array(NA, dim=c(nrow(tsdata)+nstep,ncol(tsdata),rem),list(NULL,colnames(tsdata),NULL)) 
	xhatn<-array(NA, dim=c(rem, ncol(tsdata), nstep),list(NULL,colnames(tsdata),NULL))
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
	saveRDS(xhatn,file = paste0('xhatn_GJR',ndpc,'_',nstep,'.RDS'))
	saveRDS(Y,paste0('Y_GJR',ndpc,'_',nstep,'.RDS'))
	# saveRDS(xhat,file = paste0('xhat_GJR',ndpc,'_',nstep,'.RDS'))
	#transforming into prices
	for(i in 1:nstep){
		for(g in 1:rem){
			for(p in 1:ncol(tsdata)){ if(i==1){
				zhatn[g,p,i]<-exp(xhatn[g,p,i])*(Z[1187,p]) 
			}else{
				zhatn[g,p,i]<-(exp(xhatn[g,p,i]))*(zhatn[g,p,(i-1)]) }
			} }
	}
	saveRDS(zhatn,file = paste0('zhatn_GJR',ndpc,'_',nstep,'.RDS'))
	
}
}





fdat <- as.data.frame(readxl::read_xlsx('FORECASTS.xlsx','ALL'))
ftime <- as.Date(as.data.frame(readxl::read_xlsx('FORECASTS.xlsx','Sheet1'))[,2])
#using the weightings of the portfolio which are saved as "tau" we can calculate the value of our portfolio
#V represents true values and V_N represents future sample
V<-rep(0, nrow(Z))
for(t in 1:nrow(Z)){
	for(p in 1:ncol(tsdata)){
		V[t]<-V[t]+(tau[p]*Z[t,p]) }
}
V <- ts(data = V,
				start = c(year(as.Date(pricesdata[1,1])),yday(pricesdata[1,1])),
				frequency = 365*(datarows/1820))
png('plot/vop.png',
		width = 12*580, height = 3*580, units = "px", pointsize = 85,
		bg = "white")
ts.plot(V, ylab = 'Value of Portfolio (VoP)',xlab = '', lwd = 3)
dev.off()

V_N<-array(rep(0, rem*nstep), dim=c(rem,nstep)) 
for(g in 1:rem){
	for(i in 1:nstep){
		for(p in 1:ncol(tsdata)){
			V_N[g,i]<-V_N[g,i]+(tau[p]*zhatn[g,p,i]) }
	}}
vn <- ts(t(V_N),
				 start = c(year(as.Date(ftime[1])),yday(ftime[1])),
				 frequency = 365*(datarows/1820))
vn.med <- ts(apply(V_N,2,median),
						 start = c(year(as.Date(ftime[1])),yday(ftime[1])),
						 frequency = 365*(datarows/1820))
ts.plot(as.numeric(vn.med), ylab = 'Prediction Median VoP',xlab = '')


summary(V_N[,1]);summary(V_N[,5])

vn_sum_gjr <- data.frame(t(matrix(summary(V_N[,1]),1,6)),
												 row.names = attr(summary(V_N[,1]),'names'))

xtable::xtable()


PercentageGains<-array(rep(0, rem*nstep), dim=c(rem,nstep))

for(g in 1:rem){
	for(i in 1:nstep){ 
		PercentageGains[g,i]<-(V_N[g,i]-V[1186])/V[1186]
	} }
ts.plot(ts(apply(PercentageGains,2,median),
					 start = c(year(as.Date(pricesdata[datarows,1])),
					 					yday(pricesdata[datarows,1])),
					 frequency = 365*datarows/1820))

#plotting prob density of percentage gains of the value of the portfolio (PGVP)
# for the progressing steps ahead.
for (i in seq(1,nstep,3)) {
	plot(density(PercentageGains[,i]))
	Sys.sleep(.3)
}
#plotting the median of the PGVP gains ahead
plot(apply(PercentageGains,2,median),type = 'l')

days.ah <- c(1,5,10,50,100)
DirectVaR<-sapply(days.ah, function(x) quantile(PercentageGains[,x],0.05))
xtable::xtable(DirectVaR<-setNames(data.frame((matrix(DirectVaR,1,5))),
																	 paste0('1186 + ',days.ah)))

#this is for the one-step ahead... can be adjust for how many steps ahead
#estimating ES for one-step ahead but can easily be adjusted 
#by changing PercentageGains[g,#]
EstimatedES<-array(NA,list(rem,length(days.ah)))
for (i in 1:length(days.ah)) {
	count <- 1
	for(g in 1:rem){
		if(PercentageGains[g,days.ah[i]]<=DirectVaR[i]){
			EstimatedES[count,i]<-PercentageGains[g,days.ah[i]]
			count<-count+1 }
	}
}
EES<-apply(na.omit(EstimatedES),2,mean)
xtable::xtable(EES<-setNames(data.frame((matrix(EES,1,5))),
														 paste0('1186 + ',days.ah)))



#to cross validate our model for the time period 10/09/2020 to 04/05/2021 
#we ran the code to create a 236 step ahead forecast. 
#This forecast was then adjusted to match up with the trading days 
#experienced in real life to create a rem*nsteps matrix which represents 
#the realisations for each time point
#we call this matrix lastV_N11
#NewV is a vector representing the actual 135 values observed
#Cross validating our model

##

fy_star<-c()
for(i in 1:135){
	y<-vn[i,]
	eplot<-epdfPlot(y, plot.it=FALSE)
	y_star<-lastV_N[i]
	fy_star[i]<-0
	if((y_star>=min(eplot$x)) && (y_star<max(eplot$x))){
		minabs<-abs(y_star-eplot$x) 
		fy_star[i]<-log(eplot$f.x[which(minabs==min(minabs))])
	} 
}

xtable::xtable(setNames(data.frame((matrix(summary(fy_star),1,6))),
												attr(summary(fy_star),'names')))
