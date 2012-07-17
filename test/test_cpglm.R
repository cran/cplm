###
if (FALSE) {
library(tweedie)
library(lme4)
library(rbenchmark)

options(error = recover)
#setwd("~/2011/cplm")
setwd("C:\\Documents and Settings\\CAB2007\\My Documents\\2012\\cplm")
load("./data/fineroot.RData")
load("./data/InsReserve.RData")

source("./R/classMethods.R")
source("./R/utilities.R")

link="log"
control=list()
trace=T
contrasts = NULL



cpglm <- function(formula, link = "log", data, weights, offset, 
                  subset, na.action = NULL, contrasts = NULL, 
                  control = list(), chunksize = 0, ...) {}
fm <- RLD ~ factor(Zone)
call <- match.call(cpglm,call("cpglm",quote(fm), 
  data=fineroot, chunksize = 250))
                    
dyn.load("src/cpglm_em.dll")

dyn.unload("src/cpglm_em.dll")


# profile likelihood         
fit1 <- cpglm(RLD~ factor(Zone)*factor(Stock),
  data=fineroot)
     
# residual and qq plot
parold <- par(mfrow=c(2,2), mar=c(5,5,2,1))
# 1. regular plot
r1 <- resid(fit1)/sqrt(fit1$phi)
plot(r1~fitted(fit1), cex=0.5)
qqnorm(r1, cex=0.5)
# 2. quantile residual plot to avoid overlappling
u <- ptweedie(fit1$y, fit1$p, fitted(fit1), fit1$phi)
u[fit1$y == 0] <- runif(sum(fit1$y == 0), 0, u[fit1$y == 0])
r2 <- qnorm(u)
plot(r2~fitted(fit1), cex=0.5)
qqnorm(r2, cex=0.5)

par(parold)





# profile likelihood         
fit1 <- cpglm(RLD~ factor(Zone)*factor(Stock),
  data=fineroot)


# MCEM fit

set.seed(12)
fit2 <- cpglm(RLD~ factor(Zone)*factor(Stock),
  data=fineroot, method="MCEM",

set.seed(12)
fit2 <- cpglm(RLD~ factor(Zone)*factor(Stock),
	data=fineroot, method="MCEM",

  control=list(init.size=5,sample.iter=50,
              max.size=200,fixed.size=FALSE),
  pstart=1.6)



# show the iteration history of p
plot(fit2$theta.all[,ncol(fit2$theta.all)],
     type="o", cex=0.5)

# compare the two
summary(fit1)
summary(fit2)


# insurance loss data
fit3 <- cpglm(increLoss~factor(year)+factor(lag), data=insLoss)              

# get the standard error for phi and p              
sqrt(diag(attr(fit3$vcov,"phi_p")))              
              
}