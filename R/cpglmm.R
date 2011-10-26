#######################################################
##           Compound Poisson GLMM                   ##
## Author: Wayne Zhang, actuary_zhang@hotmail.com    ##
#######################################################

cpglmm <- function(formula, link="log", data, weights, offset,
                  subset, na.action, betastart=NULL, phistart=NULL, 
                  pstart=NULL, contrasts = NULL, control = list()) {
    
  call <- match.call()  
  if (missing(data)) 
    data <- environment(formula)   
  link.power <- make.link.power(link)
  # create model frame 
  fr <- lme4:::lmerFrames(call, formula, contrasts)
  offset <- wts <- NULL
  if (length(fr$wts)) 
        wts <- fr$wts
  if (length(fr$off)) 
        offset <- fr$off

  if (!is.null(betastart)){
    if (length(betastart) != ncol(fr$X))
      stop(gettextf("number of 'betastart' is %d should 
                    equal %d (number of mean parameters)", 
                length(betastart), ncol(fr$X)), domain = NA)
    }
  if (!is.null(phistart) && length(phistart)>1) 
    stop("multiple values specified for 'phistart'")
  if (!is.null(phistart) && phistart<=0)
    stop("value of 'phistart' should be greater than 0")
  if (!is.null(pstart) && length(pstart)>1) 
    stop("multiple values specified for 'pstart'")
  if (!is.null(pstart) && (pstart<=1 || pstart>=2))
    stop("value of 'pstart' should be between 1 and 2")   
  if (is.null(pstart))
    pstart <- 1.5
  
  # generate starting values
  glmFit <- glm.fit(fr$X, fr$Y, weights = wts, offset = offset, 
      family = tweedie(var.power=pstart,link.power=link.power), 
      intercept = attr(attr(fr$mf, "terms"), "intercept") > 0)

  if (is.null(betastart))
      betastart <- glmFit$coefficients
  names(betastart) <- names(fr$fixef)
  if (is.null(phistart))
      phistart <- sum((fr$Y-glmFit$fitted.values)^2/glmFit$fitted.values^pstart)/
                    glmFit$df.residual
                  
  # get factor list and dims                  
  FL <- lme4:::lmerFactorList(formula, fr, 0L, 0L)
  ctr <- do.call(cpglmm.control, control)
  FL$dims["mxit"] <- ctr$max.iter
  FL$dims["mxfn"] <- ctr$max.fun
  dm <- lme4:::mkZt(FL, NULL)
  dm$dd["verb"] <- ctr$trace
  AGQlist <- .Call(lme4:::lme4_ghq, 1)                  
  M1 <- length(levels(dm$flist[[1]]))
  n <- ncol(dm$Zt)
  if (M1 >= n) {
        msg1 <- "Number of levels of a grouping factor for the random effects\n"
        msg3 <- "n, the number of observations"
        if (dm$dd["useSc"]) 
            stop(msg1, "must be less than ", msg3)
        else if (M1 == n) 
            message(msg1, "is *equal* to ", msg3)
    }
        
  # environment to evaluate tweedie density
  tmp.env <- new.env( )
  environment(ldtweedie) <- tmp.env
                    
  # input cpglmm class  for optimization 
    ans <- new(Class = "cpglmm", env = tmp.env, nlmodel = (~I(x))[[2]], 
        frame = fr$mf, call = call, flist = dm$flist, 
        Zt = dm$Zt, X = fr$X, y = fr$Y, pWt = unname(glmFit$prior.weights), 
        offset = unname(fr$off), Gp = unname(dm$Gp), dims = dm$dd, 
        ST = dm$ST, A = dm$A, Cm = dm$Cm, Cx = (dm$A)@x, L = dm$L, 
        deviance = dm$dev, fixef = betastart, ranef = numeric(dm$dd[["q"]]), 
        u = numeric(dm$dd[["q"]]), eta = unname(glmFit$linear.predictors), 
        mu = unname(glmFit$fitted.values), muEta = numeric(dm$dd[["n"]]), 
        var = numeric(dm$dd[["n"]]), resid = unname(glmFit$residuals), 
        sqrtXWt = as.matrix(numeric(dm$dd[["n"]])), sqrtrWt = numeric(dm$dd[["n"]]), 
        RZX = matrix(0, dm$dd[["q"]], dm$dd["p"]), RX = matrix(0, dm$dd["p"], dm$dd["p"]), 
        ghx = AGQlist[[1]], ghw = AGQlist[[2]], 
        p=pstart, phi=phistart,link.power=as.double(link.power), 
        bound.p=ctr$bound.p)

  # run optimization
  invisible(.Call("cpglmm_optimize",ans,ldtweedie))
  if (ans@dims[["cvg"]] > 6) 
        warning(lme4:::convergenceMessage(ans@dims[["cvg"]]))
  invisible(.Call("mer_update_RX",ans))                  
  ans
}        

cpglmm.control <- function(max.iter=300L,
                       max.fun=900L,               
                       bound.p=c(1.01,1.99),
                       trace=0){         
  if (!is.numeric(max.iter) || max.iter <= 0) 
        stop("value of 'max.iter' must be > 0")
  if (!is.numeric(max.fun) || max.fun <= 0) 
        stop("value of 'max.fun' must be > 0")
  if (!is.numeric(bound.p) || length(bound.p)!=2)
        stop("'bound.p' must be of length 2")
  if (min(bound.p)<1)
        stop("invalid lower bound in 'bound.p'")          
  if (!is.numeric(trace) && !is.logical(trace))
        stop("'trace' must be logical or numeric")
  bound.p <- sort(bound.p)
  
  list(max.iter=as.integer(max.iter),
       max.fun=as.integer(max.fun),
       bound.p=as.numeric(bound.p),
       trace=as.integer(trace))  
}
    
# return -2*loglike given a list x
ldtweedie <- function(x){
  -2*sum(log(dtweedie(y=x[[1]],mu=x[[2]],phi=x[[3]],p=x[[4]])))	
}
  
