if (FALSE){

# load package
library(tweedie)
library(rbenchmark)
library(ggplot2)
library(coda)



options(error = recover)
#setwd("~/2011/cplm")
setwd("C:\\Documents and Settings\\CAB2007\\My Documents\\2011\\cplm")
load("./data/fineroot.RData")
source("./R/classMethods.R")
source("./R/cpglm.R")

#dyn.load("./src/cplm.so")
link="log"
control=list()
trace=T


dyn.load("./src/cpglm.dll")  
    
mf <- match.call(bcpglm,call("bcpglm",RLD~ factor(Zone)*factor(Stock),
	data=fineroot))
                    

# default fit, 3 chains, no initial values 
set.seed(10)
fit1 <- bcpglm(RLD~ factor(Zone)*factor(Stock), data=fineroot)
gelman.diag(fit1$sims.list)
# have serious autocorrelations                             
acfplot(fit1$sims.list,lag.max=20)                             
summary(fit1)                             
                             
# supply initial values
# run the profiled likelihood approach to
# get rough estimates                             
M <- cpglm(RLD~ factor(Zone)*factor(Stock),
  data=fineroot,method="profile", 
	control=list(decimal=1,trace=FALSE))
n.chains <- 3
# inits must be a list of length n.chains           
init.vals <- vector("list",n.chains)
# each element is a vector (beta, phi, p)           
init.vals <- lapply(init.vals, function(x) 
                      c(coef(M)+rnorm(length(coef(M))),
                      M$phi+runif(1), 
                      runif(1,(M$p+1)/2,(M$p+2)/2)))
fit2 <- bcpglm(RLD~ factor(Zone)*factor(Stock), data=fineroot,
               inits=init.vals, n.chains=n.chains)
               
# run longer chains to get better mixing
fit3 <- bcpglm(RLD~ factor(Zone)*factor(Stock), data=fineroot,
               n.chains=3, n.iter=11000, n.burnin=1000,
               n.thin=10)
# since the slot sims.list is a "mcmc.list" object
# we can apply various methods defined from the package "coda"                 
gelman.diag(fit3$sims.list)
xyplot(fit3$sims.list)               
acfplot(fit3$sims.list,lag.max=20)                  
densityplot(fit3$sims.list)               
                 
# I also defined plot and summary directly for a "bcpglm" object                 
summary(fit3)
plot(fit3)                 
     


}
