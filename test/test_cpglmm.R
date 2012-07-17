if (FALSE){

cpglmm <- function(formula, link="log", data, weights, offset,
                  subset, na.action, inits=NULL,
                   contrasts = NULL, control = list(), nAGQ = 1, 
                   basisGenerators = c("tp","tpU","bsp","sp2d")) {}
    
library(ggplot2)
library(lme4)
library(amer)
library(tweedie)

setwd("C:\\Documents and Settings\\cab2007\\My Documents\\2012\\cplm")

load("./data/fineroot.RData")
source("R/cpglm.R")
#source("R/cpglmm.R")
source("R/classMethods.R")
source("R/utilities.R")
source("R/spline.R")


da <- subset(AutoClaim, IN_YY == TRUE)
da$JOBCLASS <- factor(da$JOBCLASS)
da$CLM_AMT5 <- da$CLM_AMT5 / 1000

call <- amer:::expand.call(cpglmm,call("cpglmm",CLM_AMT5 ~ tp(MVR_PTS, k = 5) + AREA + CAR_USE + MARRIED + (1|JOBCLASS), 
             data = da))
formula <- CLM_AMT5 ~ tp(MVR_PTS, k = 5) + AREA + CAR_USE + MARRIED + (1|JOBCLASS)
link <- "log"
control <- list(trace = 1)
contrasts <- NULL
nAGQ = 1
inits = NULL
doFit = TRUE


call <- amer:::expand.call(cpglmm,call("cpglmm", RLD ~ Stock * Zone + (1|Plant), data = FineRoot))
formula <- RLD ~ Stock * Zone + (1|Plant)
link <- "log"
control <- list(trace = 1)
contrasts <- NULL
nAGQ = 1
inits = NULL
doFit = TRUE
optimizer = "nlminb"



dyn.load("src/cpglmm.dll")  
a=.Call("init")



dyn.unload("src/cpglmm.dll")
f0 <- glmmPQL(RLD~  Stock * Zone, random=~1|Plant, 
              family=tweedie(var.power=1.44, link.power=0),
              data=FineRoot, niter = 10 )
f1 <- cpglmm(RLD~  Stock*Zone+ (1|Plant) , 
            link="log", data = fineroot)

f1 <- function()
  cpglmm(LOSS ~ LOGTIV + factor(BLD_COV) + factor(STR_MULT) 
         + factor(PROTECT_CD) + POL_AGE
         + CREDIT_SCORE + COLLEGE_PCT + BUS_AGE
         + CNTY_BURGLARY + CNTY_CPI + CNTY_LLEGAL + UNEMP_RATE  
         + (1|COUNTY),
         data = dat2, control = list(PQL.init = FALSE))

f2 <- function()
  cpglmm(LOSS ~ LOGTIV + factor(BLD_COV) + factor(STR_MULT) 
          + factor(PROTECT_CD) + POL_AGE
          + CREDIT_SCORE + COLLEGE_PCT + BUS_AGE
          + CNTY_BURGLARY + CNTY_CPI + CNTY_LLEGAL + UNEMP_RATE  
          + (1|COUNTY),
          data = dat2, control = list(PQL.init = TRUE))
a = f1()
b= f2()

call <- amer:::expand.call(cpglmm,call("cpglmm",LOSS ~ LOGTIV + factor(BLD_COV) + factor(STR_MULT) 
                                       + factor(PROTECT_CD) + POL_AGE
                                       + CREDIT_SCORE + COLLEGE_PCT + BUS_AGE
                                       + CNTY_BURGLARY + CNTY_CPI + CNTY_LLEGAL + UNEMP_RATE  
                                       + (1|COUNTY),
                                       data = dat2))
formula <-   LOSS ~ LOGTIV + factor(BLD_COV) + factor(STR_MULT) +
  factor(PROTECT_CD) + POL_AGE+ 
  CREDIT_SCORE + COLLEGE_PCT + BUS_AGE+
  CNTY_BURGLARY + CNTY_CPI + CNTY_LLEGAL + UNEMP_RATE  +
  (1|COUNTY)



f1 <- function()
  cpglmm(RLD~  factor(Stock) * Zone +  (1|Plant) , 
         data = FineRoot)

f2 <- function()
  cpglmm2(RLD~  factor(Stock) * Zone +  (1|Plant) , 
          data = FineRoot, control = list(PQL.init = TRUE, trace = 1))

system.time(f1())
system.time(f2())
library(rbenchmark)
benchmark(f1(), f2(), replications =10)
sigma(f1)
VarCorr(f1)









#######################################################
##           Compound Poisson GLMM                   ##
## Author: Wayne Zhang, actuary_zhang@hotmail.com    ##
#######################################################

cpglmm2 <- function(formula, link = "log", data, weights, offset,
                    subset, na.action, inits = NULL, 
                    contrasts = NULL, control = list(),
                    basisGenerators = c("tp", "tpU", "bsp", "sp2d"),
                    optimizer = "nlminb", doFit = TRUE, nAGQ = 1) {
  
  call <- amer:::expand.call()  
  if (missing(data)) 
    data <- environment(formula)   
  link.power <- make.link.power(link)
  # identify smooth terms 
  formula <- eval(call$formula)
  tf <- terms.formula(formula, specials = eval(call$basisGenerators, 
                                               parent.frame(2)))
  n.f <- length(unlist(attr(tf, "specials")))
  # create model frame and get factor list  
  if (n.f) {
    call2 <- as.list(call)[-1]
    call2 <- call2[-match(c("link", "inits", "control", 
                            "optimizer", "doFit", "nAGQ"), names(call2))]
    setup <- do.call(frFL, as.list(call2))
    fr <- setup$m$fr 
    FL <- setup$m$FL
  } else {  
    fr <- lme4:::lmerFrames(call, formula, contrasts)
    FL <- lme4:::lmerFactorList(formula, fr, 0L, 0L)
  }
  
  # set control parameters
  ctr <- do.call(cplm.control, control)
  FL$dims["mxit"] <- ctr$max.iter
  FL$dims["mxfn"] <- ctr$max.fun
  # get dims 
  dm <- lme4:::mkZt(FL, NULL)
  dm$dd["verb"] <- ctr$trace
  # update nAGQ in dd
  if ((nAGQ <- as.integer(nAGQ)) < 1) nAGQ <- 1L
  if (nAGQ %% 2 == 0) nAGQ <- nAGQ + 1L # reset nAGQ to be an odd number
  dm$dd["nAGQ"] <- as.integer(nAGQ)
  AGQlist <- .Call(lme4:::lme4_ghq, nAGQ)
  M1 <- length(levels(dm$flist[[1]]))
  n <- ncol(dm$Zt)
  q <- dm$dd[["q"]]
  d <- dm$dd[["p"]]
  # default offset and prior wts
  if (is.null(fr$wts) || length(fr$wts) == 0)  
    fr$wts <- as.double(rep(1, n)) 
  if (is.null(fr$off) || length(fr$off) == 0)
    fr$off <- as.double(rep(0, n)) 
  if (M1 >= n) {
    msg1 <- "Number of levels of a grouping factor for the random effects\n"
    msg3 <- "n, the number of observations"
    if (dm$dd["useSc"]) 
      stop(msg1, "must be less than ", msg3)
    else if (M1 == n) 
      message(msg1, "is *equal* to ", msg3)
  }
  
  # initial values
  if (!is.null(inits)){    
    check.inits.cpglmm(inits, dm$dd['p'], dm$dd['nt'])
  } else {
    # generate starting values
    inits <- cpglm.init(fr, link.power)
  }
  names(inits$beta) <- names(fr$fixef)
  
  # input cpglmm class for optimization 
  # muEta and var are not initialized so that lmer is fitted when PQL.init is TRUE
  ans <- new(Class = "cpglmm", env = new.env( ), nlmodel = (~I(x))[[2]], 
             frame = fr$mf, call = call, flist = dm$flist, 
             Zt = dm$Zt, X = fr$X, y = as.numeric(fr$Y), pWt = fr$wts, 
             Gp = unname(dm$Gp), dims = dm$dd, 
             ST = dm$ST, A = dm$A, Cm = dm$Cm, Cx = (dm$A)@x, L = dm$L, 
             deviance = dm$dev, fixef = inits$beta, ranef = numeric(q), 
             u = numeric(q), eta = numeric(n), 
             mu = numeric(n), resid = numeric(n), 
             sqrtXWt = as.matrix(numeric(n)), sqrtrWt = numeric(n), 
             RZX = matrix(0, q, d), RX = matrix(0, d, d), 
             ghx = AGQlist[[1]], ghw = AGQlist[[2]], 
             p = inits$p, phi = inits$phi, link.power = as.double(link.power), 
             bound.p = ctr$bound.p, formula = formula, contrasts = contrasts,
             model.frame = fr$mf, inits = inits, vcov = matrix(0, d, d))
  
  # run the PQL method to improve initial values
  if (ctr$PQL.init) {
    ST <- ans@ST
    eta <- as.numeric(fr$X %*% inits$beta  + fr$off)    
    fam <- tweedie(1.5, link.power)        
    while (1) {
      etaold <- eta
      mu.eta.val <- fam$mu.eta(eta)
      mu <- fam$linkinv(eta)
      ans@y <- eta + (fr$Y - mu)/mu.eta.val - fr$off
      ans@pWt <- fr$wts * mu.eta.val^2/fam$variance(mu)
      .Call(lme4:::mer_optimize, ans)
      .Call(lme4:::mer_update_ranef, ans)
      .Call(lme4:::mer_update_mu, ans)
      eta <- ans@eta + fr$off
      if (sum((eta - etaold)^2) < 1e-06 * sum(eta^2)) 
        break    
    }
    # reset y and pWt 
    ans@y <- as.numeric(fr$Y)
    ans@pWt <- fr$wts        
    # reset ST if the estimate from PQL is degenerate
    if (any(sapply(ans@ST, det) < 1e-8)) ans@ST <- ST
  }    
  
  # initialize muEta and var to fit glmer  
  ans@muEta <-  numeric(n) 
  ans@var <- numeric(n)
  ans@offset <- fr$off
  
  # return cpglmm object if model fitting is turned off
  if (!doFit)  return(ans)
  
  # run optimization
  cpglmm_dev <- function(parm){
    # set parameters
    .Call("cpglmm_setPars", ans, parm)
    # update u
    .Call("cpglmm_update_u", ans)
    #invisible(.Call(lme4:::mer_update_ranef, ans))  
    u <- ans@u
    nlvl <- length(levels(ans@flist[[1]]))
    nre <- length(u)/nlvl
    flist.int <- as.integer(ans@flist[[1]])
    out <-  matrix(0, nlvl, nAGQ^nre)
    q <- unname(ans@dims["q"])
    pos <- 1
    for (i in 1:nAGQ){
      z <- rep(ans@ghx[i], nlvl)
      const <- log(ans@ghw[i])  + ans@ghx[i]^2  
      ux <- as.vector(u + sqrt(2) * sqrt(ans$phi) * solve(ans@L, z, "L"))
      ans@u <- as.vector(ux)
      .Call("cpglmm_update_mu", ans)
      ly <- log(dtweedie(ans@y, mu = ans@mu, phi = ans@phi, power = ans@p))
      ly2 <- tapply(ly, flist.int, sum) 
      lu <- tapply(ux^2, rep(1:nlvl, nre), sum)/(2 * ans@phi)
      out[, pos] <- const + ly2 - lu  
      pos <- pos + 1
    }
    maxx <- apply(out, 1, max)
    llik <- sum(log(rowSums(exp(sweep(out, 1, maxx)))) + maxx)
    -2 * llik + 2*log(det(ans@L)) 
  }
  
  # initial values for theta, beta, log(phi), p, 
  parm <- c(.Call(lme4:::mer_ST_getPars, ans), ans$fixef, log(ans$phi), ans$p) 
  parm <- unname(parm)
  # set bounds
  n.parm <- length(parm)
  lower <- rep(-Inf, n.parm)
  upper <- rep(Inf, n.parm)
  # reset bounds for diagonal elements of theta
  lower[1:sum(sapply(ans$ST, ncol))] <- 0
  # reset bounds for p
  lower[n.parm] <- ans$bound.p[1]
  upper[n.parm] <- ans$bound.p[2]
  
  # run optimizations
  rslt <- cplm:::cplm_optim(parm, cpglmm_dev, lower = lower, upper = upper, 
                            control = ctr, optimizer = optimizer)
  ans@dims[["cvg"]] <- as.integer(rslt$convergence)
  if (rslt$convergence)  warning(rslt$message)
  
  # update ans using the found optima 
  cpglmm_dev(rslt$par)    
  
  # update random effects 
  invisible(.Call(lme4:::mer_update_ranef, ans))  
  # update sigmaML to be used in postVar
  dev <- ans@deviance 
  dev['sigmaML'] <- sqrt(ans@phi)
  ans@deviance <- dev 
  # update vcov in ans
  invisible(.Call(lme4:::mer_update_RX, ans))
  ans@vcov <- vcov(ans)
  ans
  
}        






}

