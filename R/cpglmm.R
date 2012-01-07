#######################################################
##           Compound Poisson GLMM                   ##
## Author: Wayne Zhang, actuary_zhang@hotmail.com    ##
#######################################################

cpglmm <- function(formula, link="log", data, weights, offset,
                  subset, na.action, inits=NULL, 
                  contrasts = NULL, control = list(),
                  basisGenerators = c("tp","tpU","bsp","sp2d"),
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
    call2 <- call2[-match(c("link","inits","control", 
                            "optimizer", "doFit", "nAGQ"), names(call2))]
    #setup <- do.call(amer:::amerSetup, as.list(call2))     
    setup <- do.call(frFL, as.list(call2))
    fr <- setup$m$fr 
    FL <- setup$m$FL
  } else {  
    fr <- lme4:::lmerFrames(call, formula, contrasts)
    FL <- lme4:::lmerFactorList(formula, fr, 0L, 0L)
  }
  
  # set control parameters
  ctr <- do.call(cpglmm.control, control)
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
  # default offset and prior wts
  if (is.null(fr$wts) || length(fr$wts) == 0)  
    wts <- as.double(rep(1,n)) else 
    wts <- fr$wts 
  if (is.null(fr$off) || length(fr$off) == 0)
    off <- as.double(rep(0,n)) else 
    off <- fr$off
  if (M1 >= n) {
        msg1 <- "Number of levels of a grouping factor for the random effects\n"
        msg3 <- "n, the number of observations"
        if (dm$dd["useSc"]) 
            stop(msg1, "must be less than ", msg3)
        else if (M1 == n) 
            message(msg1, "is *equal* to ", msg3)
    }
  
  # check initial values
  if (!is.null(inits)){
    check.inits.cpglmm(inits, dm$dd['p'], dm$dd['nt'])
    betastart <- inits$beta 
    names(betastart) <- names(fr$fixef)
    phistart <- inits$phi 
    pstart <- inits$p
  } else {
    pstart <- 1.5
    # generate starting values
    glmFit <- glm.fit(fr$X, fr$Y, weights = wts, offset = off, 
        family = tweedie(var.power = pstart, link.power = link.power), 
        intercept = attr(attr(fr$mf, "terms"), "intercept") > 0)
    betastart <- glmFit$coefficients
    names(betastart) <- names(fr$fixef)
    phistart <- sum((fr$Y - glmFit$fitted.values)^2 / glmFit$fitted.values^pstart) /
                    glmFit$df.residual
    inits <- list(beta = betastart, phi = phistart, p = pstart)
  }
  # input cpglmm class  for optimization 
  ans <- new(Class = "cpglmm", env = new.env( ), nlmodel = (~I(x))[[2]], 
        frame = fr$mf, call = call, flist = dm$flist, 
        Zt = dm$Zt, X = fr$X, y = as.numeric(fr$Y), pWt = wts, 
        offset = off, Gp = unname(dm$Gp), dims = dm$dd, 
        ST = dm$ST, A = dm$A, Cm = dm$Cm, Cx = (dm$A)@x, L = dm$L, 
        deviance = dm$dev, fixef = betastart, ranef = numeric(dm$dd[["q"]]), 
        u = numeric(dm$dd[["q"]]), eta = numeric(dm$dd[["n"]]), 
        mu = numeric(dm$dd[["n"]]), muEta = numeric(dm$dd[["n"]]), 
        var = numeric(dm$dd[["n"]]), resid = numeric(dm$dd[["n"]]), 
        sqrtXWt = as.matrix(numeric(dm$dd[["n"]])), sqrtrWt = numeric(dm$dd[["n"]]), 
        RZX = matrix(0, dm$dd[["q"]], dm$dd["p"]), RX = matrix(0, dm$dd["p"], dm$dd["p"]), 
        ghx = AGQlist[[1]], ghw = AGQlist[[2]], 
        p = pstart, phi = phistart, link.power = as.double(link.power), 
        bound.p = ctr$bound.p, formula = formula, contrasts = contrasts,
        model.frame = fr$mf, inits = inits, vcov = matrix(0, dm$dd["p"], dm$dd["p"]))
  # return cpglmm object if model fitting is turned off
  if (!doFit)
    return(ans)
  
  # run optimization
  if (optimizer == "nlminb") {
    invisible(.Call("cpglmm_optimize",ans))
    if (ans@dims[["cvg"]] > 6) 
        warning(lme4:::convergenceMessage(ans@dims[["cvg"]]))
  } else {    
    # function to be fed to optimizers (return deviance)
    # parm: theta+beta+log(phi)+p
    cpglmm_dev <- function(parm){             
      .Call("cpglmm_update_dev",ans, parm) 
    }
    # initial values theta, beta, log(phi), p, 
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

    if (optimizer == "bobyqa"){ 
      rslt <- bobyqa(parm, cpglmm_dev, lower = lower, upper = upper, 
         control = list(iprint = ctr$trace,rhobe = 0.02, rhoend = 2e-7))
      ans@dims[["cvg"]] <- as.integer(rslt$ierr)
      if (rslt$ierr) warning(rslt$msg)
    } else 
    if (optimizer == "L-BFGS-B"){
      rslt <- optim(parm, cpglmm_dev, gr = NULL, 
          method = "L-BFGS-B", lower = lower, upper = upper, 
          control = list(trace = ctr$trace, maxit = ctr$max.iter,
                         maxfun = ctr$max.fun))
      ans@dims[["cvg"]] <- as.integer(rslt$convergence)
      if (rslt$convergence) warning(rslt$message)
    }

    # update ans using the found optima 
    invisible(.Call("cpglmm_update_dev", ans, rslt$par))    
  }
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




  
