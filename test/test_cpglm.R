

###
if (FALSE) {
library(tweedie)
library(rbenchmark)

options(error = recover)
#setwd("~/2011/cplm")
setwd("C:\\Documents and Settings\\CAB2007\\My Documents\\2011\\cplm")
load("./data/fineroot.RData")
source("./R/classMethods.R")
#dyn.load("./src/cplm.so")
link="log"
control=list()
trace=T


mf <- match.call(cpglm,call("cpglm",RLD~ factor(Zone)*factor(Stock),
  data=fineroot,control=list(maxit=150,sample.iter=20),
      pstart=1.4))
                    
dyn.load("./src/cpglm.dll")
dyn.unload("./src/cpglm.dll")

# MCEM fit
set.seed(10)
fit1 <- cpglm(RLD~ factor(Zone)*factor(Stock),
	data=fineroot,
  control=list(init.size=5,sample.iter=50,
              max.size=2000,fixed.size=FALSE),
  pstart=1.6)

# profile likelihood         
fit2 <- cpglm(RLD~ factor(Zone)*factor(Stock),
	data=fineroot,method="profile", 
	control=list(decimal=1))      

# compare the two 
summary(fit1)
summary(fit2)

}