if (FALSE){

cpglmm <- function(formula, link="log", data, weights, offset,
                  subset, na.action, inits=NULL,
                   contrasts = NULL, control = list(),
                   basisGenerators = c("tp","tpU","bsp","sp2d")) {}
    
library(ggplot2)
library(lme4)
library(amer)
library(tweedie)

setwd("C:\\Documents and Settings\\cab2007\\My Documents\\2011\\cplm")

load("./data/fineroot.RData")
source("R/cpglm.R")
#source("R/cpglmm.R")
source("R/bcpglm.R")
source("R/classMethods.R")
source("R/utilities.R")
source("R/spline.R")
dyn.load("src/cplm.dll")
dyn.unload("src/cplm.dll")


fineroot$x1 <- rnorm(nrow(fineroot))
fineroot$x2 <- rnorm(nrow(fineroot))
call <- amer:::expand.call(cpglmm,call("cpglmm",RLD~  Stock + Spacing +  (1|Plant) + (1|Zone), 
             data = fineroot))
formula <-   RLD~   Stock + Spacing +  (1|Plant) + (1|Zone)
link <- "log"
control <- list()
contrasts <- NULL

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

