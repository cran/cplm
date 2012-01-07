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
source("R/cpglm.R")
source("R/cpglmm.R")
source("R/bcpglmm.R")
source("R/classMethods.R")
source("R/utilities.R")


#dyn.load("./src/cplm.so")
link="log"
control=list()
trace=T
n.chains=3
n.iter=5000
n.burnin=1000
n.sims=1000
n.thin=max(1, floor(n.chains * (n.iter - n.burnin) / n.sims))
n.report=10000
bound.p=c(1.01,1.99)
bound.phi=100
prior.Sigma = NULL
prior.beta.mean=NULL 
prior.beta.var=NULL
inits = NULL
tune.iter=5000
n.tune=20
tune.weight=0.25



bcpglmm <- function(formula, link = "log", data, inits = NULL,
                   weights, offset, subset, na.action, contrasts = NULL, 
                   n.chains=3, n.iter=2000, n.burnin=floor(n.iter/2),
                   n.thin=max(1, floor(n.chains * (n.iter - n.burnin) / n.sims)),
                   n.sims=1000, n.report=1000, prior.beta.mean, prior.beta.var, 
                   phi.shape=0.001, phi.scale=0.001, bound.p=c(1.01,1.99),...) {}    

set.seed(10)
fineroot$a <- rnorm(nrow(fineroot))
formula <- RLD~ Zone+(1+a|Plant)+(1|Stock) 
contrasts = NULL
call <- match.call(bcpglmm,call("bcpglmm", RLD~ Zone+(1+a|Plant)+(1|Stock),data=fineroot))


dyn.load("src/bcpglmm.dll") 
a= .Call("init")


a= .Call("finish")
dyn.unload("src/bcpglmm.dll")

set.seed(10)
fit1 <- bcpglmm(RLD~Zone*Stock+(1|Plant), data = fineroot,  
                   n.chains=3, n.iter=25000, n.burnin=5000,
                   n.sims=2000, n.report=5000)

fit2 <- cpglmm(RLD~Zone*Stock+(1|Plant), data = fineroot)


x=rnorm(10)
y=rnorm(10)
z=rnorm(10)
a=cov(cbind(x,y,z))
x2= rnorm(3)
m = rnorm(3)
solve(a)

xx=cbind(x,y,z)
a= diag(1,nrow=3)
b=.Call("inv",m=as.integer(3),n=as.integer(3),x=a)

t(x2-m)%*%solve(a)%*%(x2-m) *(-0.5)

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
