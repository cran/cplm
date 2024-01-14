### R code from vignette source 'cplm.Rnw'

###################################################
### code chunk number 1: options
###################################################
options(prompt = "R> ", digits = 4, show.signif.stars = FALSE)


###################################################
### code chunk number 2: cplm.Rnw:48-70
###################################################
mypar <- function(){
  parold <- par(mar = c(3.5, 2.5, 2, 0.5), tck = -0.02, 
    mgp = c(1.3, 0.2, 0), cex.axis = 0.7, cex.lab = 0.7)
  parold
}

# no margin in lattice
theme.novpadding <-
   list(layout.heights =
        list(top.padding = 0,
       main.key.padding = 0,
 	    key.axis.padding = 0,
 	    axis.xlab.padding = 0,
 	    xlab.key.padding = 0,
 	    key.sub.padding = 0,
 	    bottom.padding = 0),
        layout.widths =
        list(left.padding = 0,
 	    key.ylab.padding = 0,
 	    ylab.axis.padding = 0,
 	    axis.key.padding = 0,
 	    right.padding = 1))


###################################################
### code chunk number 3: cplm.Rnw:116-121
###################################################
library(cplm)
obj <- new("cplm")  # create an object of class "cplm"
names(obj)          # or slotNames(obj)
obj$inits           # or obj@inits
obj[["inits"]]      # or slot(obj, "inits")


###################################################
### code chunk number 4: cplm.Rnw:125-126 (eval = FALSE)
###################################################
## showMethods(classes = "cpglm")


###################################################
### code chunk number 5: cplm.Rnw:129-132 (eval = FALSE)
###################################################
## class ? bcplm
## method ? resid("cpglm") # the "resid" method for class "cpglm" 
## method ? resid("cpglmm")


###################################################
### code chunk number 6: cplm.Rnw:174-177 (eval = FALSE)
###################################################
## da <- subset(AutoClaim, IN_YY == 1) # use data in the Yip and Yau paper
## da <- transform(da, CLM_AMT5 = CLM_AMT5/1000,
##                 INCOME = INCOME/10000)


###################################################
### code chunk number 7: cplm.Rnw:180-181 (eval = FALSE)
###################################################
## summary(da$CLM_AMT5); sum(da$CLM_AMT5 == 0)/nrow(da)


###################################################
### code chunk number 8: cplm.Rnw:190-192 (eval = FALSE)
###################################################
## P1 <- cpglm(CLM_AMT5 ~ CAR_USE + MARRIED + AREA + MVR_PTS, data = da)
## summary(P1)


###################################################
### code chunk number 9: cplm.Rnw:220-221 (eval = FALSE)
###################################################
## c(coef(P1), p = P1$p, phi = P1$phi)


###################################################
### code chunk number 10: cplm.Rnw:234-236 (eval = FALSE)
###################################################
## P1.big <- cpglm(CLM_AMT5 ~ CAR_USE + MARRIED + AREA + MVR_PTS, 
##     data = da, chunksize = 500)


###################################################
### code chunk number 11: cplm.Rnw:258-259
###################################################
head(FineRoot)


###################################################
### code chunk number 12: cplm.Rnw:263-264
###################################################
load("cpglmm.RData")


###################################################
### code chunk number 13: cplm.Rnw:266-268 (eval = FALSE)
###################################################
## f0 <- cpglmm(RLD ~ Stock * Zone +  (1|Plant), data = FineRoot)
## summary(f0)


###################################################
### code chunk number 14: cplm.Rnw:295-296
###################################################
inherits(f0, "mer")


###################################################
### code chunk number 15: cplm.Rnw:299-302
###################################################
fixef(f0)  #coefficients
ranef(f0)
VarCorr(f0)  #variance components


###################################################
### code chunk number 16: cplm.Rnw:306-308 (eval = FALSE)
###################################################
## FineRoot$Plant <- factor(FineRoot$Plant)
## f1 <- cpglmm(RLD ~ Stock + Spacing + (1|Plant), data = FineRoot)


###################################################
### code chunk number 17: cplm.Rnw:311-315 (eval = FALSE)
###################################################
## f2 <- update(f1, . ~ . + (1|Zone))
## f3 <- update(f1, . ~ . + (1|Plant:Zone))
## # test the additional random effect
## anova(f1, f2)


###################################################
### code chunk number 18: cplm.Rnw:328-329 (eval = FALSE)
###################################################
## anova(f1, f3)


###################################################
### code chunk number 19: cplm.Rnw:343-346 (eval = FALSE)
###################################################
## f4 <- cpglmm(RLD ~  Stock * Zone +  (1|Plant), 
##             data = FineRoot, optimizer = "L-BFGS-B", 
##             control = list(trace = 2, PQL.init = FALSE))


###################################################
### code chunk number 20: cplm.Rnw:366-368 (eval = FALSE)
###################################################
## f5 <- cpglmm(RLD ~  Stock * Zone +  (1|Plant), 
##             data = FineRoot, nAGQ = 10)


###################################################
### code chunk number 21: cplm.Rnw:376-377 (eval = FALSE)
###################################################
## P2 <- cpglmm(CLM_AMT5 ~ CAR_USE + MARRIED + AREA + tp(MVR_PTS), data = da)


###################################################
### code chunk number 22: cplm.Rnw:388-389 (eval = FALSE)
###################################################
## plotF(P2, rug = FALSE)


###################################################
### code chunk number 23: cplm.Rnw:428-432 (eval = FALSE)
###################################################
## set.seed(10)
## B1 <- bcplm(increLoss ~ factor(year) + factor(lag), data = ClaimTriangle, 
##             n.chains = 2, n.iter = 7000, tune.iter = 4000,
##             n.burnin = 2000, n.thin = 5, bound.p = c(1.1, 1.95))


###################################################
### code chunk number 24: cplm.Rnw:456-457 (eval = FALSE)
###################################################
## summary(gelman.diag(B1$sims.list)[[1]][, 1])


###################################################
### code chunk number 25: cplm.Rnw:464-466 (eval = FALSE)
###################################################
## xyplot(B1$sims.list[, c(1:2, 20, 21)], xlab = NULL)
## densityplot(B1$sims.list[, c(1:2, 20, 21)], ylab = NULL)


###################################################
### code chunk number 26: cplm.Rnw:497-501 (eval = FALSE)
###################################################
## set.seed(10)
## B2 <- bcplm(RLD ~ Stock * Zone + (1|Plant), 
##             data = FineRoot, n.iter = 11000, 
##             n.burnin = 1000, n.thin = 10)


###################################################
### code chunk number 27: cplm.Rnw:504-505 (eval = FALSE)
###################################################
## summary(B2)


###################################################
### code chunk number 28: cplm.Rnw:551-555 (eval = FALSE)
###################################################
## da <- transform(da, P1 = fitted(P1), P2 = fitted(P2), P3 = fitted(P3))
## gg <- gini(loss = "CLM_AMT5", score  = paste("P", 1:3, sep = ""), 
##            data = da)
## gg


###################################################
### code chunk number 29: cplm.Rnw:574-576 (eval = FALSE)
###################################################
## theme_set(theme_bw())
## plot(gg, overlay = FALSE)


