if (FALSE){

# load package
library(tweedie)
library(rbenchmark)
library(ggplot2)
library(coda)
library(lme4)


options(error = recover)
#setwd("~/2011/cplm")
setwd("C:\\Documents and Settings\\CAB2007\\My Documents\\2011\\cplm")
load("./data/fineroot.RData")
load("./data/insLoss.RData")
source("./R/classMethods.R")
source("./R/cpglm.R")
source("./R/utilities.R")
source("./R/bcpglm.R")

#dyn.load("./src/cplm.so")
link="log"
control=list()
trace=T
n.chains=3
n.iter=20
n.burnin=10
n.sims=10
n.thin=max(1, floor(n.chains * (n.iter - n.burnin) / n.sims))
n.report=10000
inits=NULL
bound.p=c(1.01,1.99)
phi.shape=0.001
phi.scale=0.001
tune.iter=5000
n.tune=10
tune.weight=0.25
bound.phi =100
prior.beta.var=NULL
prior.beta.mean=NULL

dyn.load("src/cpglm_bayes.dll")  




bcpglm <- function(formula, link = "log", data, inits = NULL,
                   weights, offset, subset, na.action, contrasts = NULL, 
                   n.chains=3, n.iter=2000, n.burnin=floor(n.iter/2),
                   n.thin=max(1, floor(n.chains * (n.iter - n.burnin) / n.sims)),
                   n.sims=1000, n.report=1000, prior.beta.mean, prior.beta.var, 
                   phi.shape=0.001, phi.scale=0.001, bound.p=c(1.01,1.99),...) {}    

set.seed(10)

mf <- match.call(bcpglm,call("bcpglm",increLoss~ factor(year)+factor(lag),
                             data=insLoss))




dyn.load("src/cpglm_bayes.dll")             

dyn.unload("src/cpglm_bayes.dll")             

# fit the fineroot data with Bayesian models
# first use the latent variable approach
set.seed(10)
fit1 <- bcpglm(RLD~ factor(Zone)*factor(Stock), data=fineroot,
               n.iter=11000, n.burnin=1000,
               n.thin=10, n.report=5000, method="latent")
gelman.diag(fit1$sims.list)
# diagnostic plots                             
acfplot(fit1$sims.list,lag.max=20)
xyplot(fit1$sims.list)                              
densityplot(fit1$sims.list)               

# summary results
summary(fit1)                             

# now try direct density evaluation (This is much slower)
set.seed(10)
fit2 <- bcpglm(RLD~ factor(Zone)*factor(Stock), data=fineroot,
               n.iter=11000, n.burnin=1000,
               n.thin=10, n.report=5000)
gelman.diag(fit2$sims.list)
summary(fit2)                             


# now fit the Bayesian model to an insurance loss triangle 
# (see Peters et al. 2009)
fit3 <- bcpglm(increLoss~ factor(year)+factor(lag), data=insLoss,
               tune.iter=5000, n.tune = 10,
               n.iter=11000, n.burnin=1000,
               n.thin=10, n.report=5000, bound.p=c(1.1,1.95))
gelman.diag(fit3$sims.list)
summary(fit3)                             



}
