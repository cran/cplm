if (FALSE){

cpglmm <- function(formula, link="log", data, weights, offset,
                  subset, na.action, betastart=NULL, phistart=NULL, 
                  pstart=NULL, contrasts = NULL, control = list()) {}
    
library(ggplot2)
library(lme4)
library(tweedie)

setwd("C:\\Documents and Settings\\cab2007\\My Documents\\2011\\cplm")

load("./data/fineroot.RData")
source("R/cpglm.R")
source("R/cpglmm.R")
source("R/bcpglm.R")
source("R/classMethods.R")
dyn.load("src/cplm.dll")
dyn.unload("src/cplm.dll")

call <- match.call(cpglmm,call("cpglmm",RLD~  Zone*Stock + (1|Plant), 
            link="log", data = fineroot))  
formula <-   RLD~   Zone*Stock + (1|Plant)
link <- "log"
control <- list()

dyn.load("src/cpglmm_lap.dll")  
a=.Call("init")

f0 <- glmmPQL(RLD~  Stock*Zone, random=~1|Plant, 
              family=tweedie(var.power=1.414, link.power=0),
              data=fineroot )
f1 <- cpglmm(RLD~  Stock*Zone+ (1|Plant) , 
            link="log", data = fineroot)

sigma(f1)
VarCorr(f1)

}

