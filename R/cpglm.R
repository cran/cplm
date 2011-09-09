#######################################################
##             compound Poisson GLM                  ##
## Author: Wayne Zhang, actuary_zhang@hotmail.com    ##
#######################################################

cpglm <- function(formula, link = "log", data, weights, offset, 
                  subset, na.action, betastart=NULL, phistart=NULL, 
                  pstart=NULL, contrasts = NULL, control=list(),
                  method ="MCEM", ...) {
  
  call <- match.call()  
  if (missing(data)) 
    data <- environment(formula)   
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights",
               "na.action", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  Y <- model.response(mf, "any")
  X <- if (!is.empty.model(mt)) 
        model.matrix(mt, mf, contrasts)
  weights <- as.vector(model.weights(mf))
  offset <- as.vector(model.offset(mf))
  link.power <- make.link.power(link)

  if (!is.null(weights) && !is.numeric(weights)) 
        stop("'weights' must be a numeric vector")
    if (!is.null(weights) && any(weights <= 0)) 
        stop("negative or zero weights not allowed")
  if (!is.null(offset)) {
    if (length(offset) != NROW(Y)) 
      stop(gettextf("number of 'offset' is %d should 
                    equal %d (number of observations)", 
                length(offset), NROW(Y)), domain = NA)
    }
  if (!is.null(betastart)){
    if (length(betastart) != ncol(X))
      stop(gettextf("number of 'betastart' is %d should 
                    equal %d (number of mean parameters)", 
                length(betastart), ncol(X)), domain = NA)
    }
  if (!is.null(phistart) && length(phistart)>1) 
    stop("multiple values specified for 'phistart'")
  if (!is.null(phistart) && phistart<=0)
    stop("value of 'phistart' should be greater than 0")
  if (!is.null(pstart) && length(pstart)>1) 
    stop("multiple values specified for 'pstart'")
  if (!is.null(pstart) && (pstart<=1 || pstart>=2))
    stop("value of 'pstart' should be between 1 and 2")   

  if (method=="MCEM")
    cpfit <- cpglm_em(X,Y,weights=weights,offset=offset,
                     link.power=link.power,
                     betastart=betastart,phistart=phistart,pstart=pstart,
                     intercept=attr(mt, "intercept") > 0L,control=control)
  if (method=="profile")    
    cpfit <- cpglm_profile(X,Y,weights=weights,offset=offset,
                     link.power=link.power,contrasts=contrasts,control=control,
                      intercept=attr(mt, "intercept") > 0L)
  
  class(mt) <- "terms"
  ans <- new("cpglm", 
             coefficients=cpfit$coefficients, 
             residuals=cpfit$residuals,
             fitted.values=cpfit$fitted.values,
             linear.predictors=cpfit$linear.predictors,
             weights=cpfit$weights,
             df.residual=cpfit$df.residual,
             deviance=cpfit$deviance,
             aic=cpfit$aic,
             offset=cpfit$offset,
             prior.weights=cpfit$prior.weights,               
             call=call,
             formula=formula,
             data=data,             
             control=cpfit$control,
             contrasts=contrasts,
             p=cpfit$p,
             phi=cpfit$phi,             
             theta=cpfit$theta,
             theta.all=cpfit$theta.all,
             vcov=cpfit$vcov,
             iter=cpfit$iter,
             converged=cpfit$converged,
             method=method,
             y=Y,
             link.power=link.power,
             na.action=attr(mf, "na.action"),
             model.frame = mf
             )  
  return(ans)
}


# function to run the MCEM 
cpglm_em <- function(X,Y,weights=NULL,offset=NULL,
                      link.power=0,
                      betastart,phistart,pstart,
                      intercept = TRUE,
                      control=list()){
    # set control options                        
    control <- do.call("cpglm.control", control)                   
    if (!is.null(pstart)){
      if (pstart<control$bound.p[1] || pstart>control$bound.p[2])
        stop ("value of 'pstart' outside the 'control$bount.p'")
    } 
    X <- as.matrix(X)          
    # get names
    xnames <- dimnames(X)[[2L]]
    ynames <- if (is.matrix(Y)) 
        rownames(Y) else 
        names(Y)
    
    # default weights and offsets if NULL    
    nobs <- NROW(Y)
    if (is.null(weights))     
      weights <- rep.int(1, nobs)
    if (is.null(offset)) 
        offset <- rep.int(0, nobs)          
    # generating starting values if necessary
    if (is.null(pstart)) 
      pstart <- sum(control$bound.p)/2
    if (is.null(betastart) || is.null(phistart)) {
      fit.start <- glm(Y~-1+X,weights=weights,offset=offset,
                  family=tweedie(var.power=pstart,
                                 link.power=link.power))
      if (is.null(betastart))
        betastart <- as.numeric(fit.start$coefficients)
      if (is.null(phistart))
        phistart <- sum(residuals(fit.start,"pearson")^2)/
          df.residual(fit.start)
    }
    
    out <- .Call("cpglm_em",
                 X=as.double(X),
                 Y=as.double(Y),
                 ygt0= as.integer(which(Y>0L)),
                 offset=as.double(offset),
                 weights=as.double(weights),
                 beta=as.double(betastart),
                 phi=as.double(phistart),
                 p=as.double(pstart),
                 link.power=as.double(link.power),
                 bound=as.double(control$bound.p),
                 init.size=as.integer(control$init.size), 
                 sample.iter=as.integer(control$sample.iter),
                 max.iter=as.integer(control$max.iter),
                 epsilon1=as.double(control$epsilon1),
                 epsilon2=as.double(control$epsilon2),
                 alpha=as.double(control$alpha),                   
                 ck = as.double(control$k),                                    
                 fixed.size=as.integer(control$fixed.size),
                 trace=as.integer(control$trace),
                 max.size=as.integer(control$max.size),
                 beta.step=as.integer(control$beta.step))
    out$vcov <- svd.inv(out$hess)
    out <- out[!(names(out)=="hess")]
    out$df.residual <- nrow(X) - ncol(X)                         
    out$deviance <- sum(tweedie.dev(Y,out$fitted.values, out$p)) 
    out$aic <- -2 * sum(log(dtweedie(Y, mu = out$fitted.values, 
                phi = out$phi, power = out$p))) + 2*(ncol(X) +2)
    out$prior.weights <- weights
    out$offset <- offset 
    out$converged <- as.logical(out$converged)                       
    out$control <- control
    names(out$coefficients) <- xnames
    names(out$residuals) <- names(out$fitted.values) <-
      names(out$linear.predictors) <- names(out$weights) <- ynames
    return(out)        
}   

# function to implement the automatic profile likelihood approach 
cpglm_profile <- function(X,Y,weights=NULL,offset=NULL,
                      link.power=0, intercept=TRUE, 
                      contrasts, control=list()){
  control <- do.call("cpglm.control", control)
  if (control$trace) {
        cat("---\n This function is based on 'tweedie.profile' in the 'tweedie' package;\n")
        cat(" If it fails, try using  method=\"series\"\n")
        cat(" rather than the default  method=\"inversion\"\n")
        cat(" Another possible reason for failure is the range of p:\n")
        cat(" Try a different boundary for 'bound.p'\n---\n")
  }  
  pl <- ceiling(control$bound.p[1]*10)/10
  pu <- floor(control$bound.p[2]*10)/10    
  for (i in 1:control$decimal){  
    p.vec <- seq(pl,pu,by=10^(-i))
    fit <- tweedie_profile(X=X, Y=Y, weights=weights, offset=offset,
                           link.power=link.power,p.vec=p.vec,
                           verbose=control$trace,intercept=intercept)
    cc <- 10^i                            
    if (i==control$decimal)
      break else{
      pl <- max(ceiling(control$bound.p[1]*10*cc)/(10*cc), 
               fit$p.max-1/cc+1/(cc*10))                 
      pu <- min(floor(control$bound.p[2]*10*cc)/(10*cc), 
               fit$p.max+1/cc-1/(cc*10))     
     }                            
  }
  # fit glm using the maximized p
  fit2 <- glm.fit(X,Y,weights=weights,offset=offset,
                  family=tweedie(var.power=fit$p.max,
                                 link.power=link.power),
                  intercept=intercept)  
  class(fit2) <- "glm"
  out <- c(list(
             deviance=sum(tweedie.dev(Y, fit2$fitted.values,fit$p.max)),
             aic=-2*fit$L.max+2*(fit2$rank+2),
             control=control,
             p=fit$p.max,
             phi=fit$phi.max,             
             theta=c(fit2$cofficients,fit$p.max,fit$phi.max),
             theta.all=matrix(c(fit2$cofficients,fit$p.max,fit$phi.max),
                              nrow=1),
             vcov=vcov(fit2),
             offset=offset),
             fit2[c("coefficients","residuals","fitted.values",
                    "linear.predictors","iter","weights",
                    "prior.weights","df.residual","converged")])  
  return(out)  
}               


# function to compute log density 
dtweedie.nlogl <- function(phi, y, mu, power) {
    ans <- -2 * sum(log(dtweedie(y = y, mu = mu, phi = phi, 
        power = power)))
    if (is.infinite(ans)) {
        ans <- sum(tweedie.dev(y = y, mu = mu, power = power))/length(y)
    }    
    #attr(ans, "gradient") <- dtweedie.dldphi(y = y, mu = mu, 
    #    phi = phi, power = power)
    ans
}
    
tweedie_profile <- function (X,Y,weights=NULL,offset=NULL,
                       p.vec = NULL, link.power=0,  
                       method = "inversion",  verbose = FALSE,
                       intercept =TRUE) {
    if (is.logical(verbose)) 
        verbose <- as.numeric(verbose)    
    np <- length(p.vec)
    if (np<1)
      stop ("'p.vec' must have at least one element")
    nY <- length(Y)                                     
    L <- phi.vec <- rep(NA, np)
                                         
    for (i in (1:np)) {
        p <- p.vec[i]
        if (verbose) 
            cat(paste("p= ", p, "\n", sep = ""))
        catch.possible.error <- try(fit.model <- glm.fit(x = X, 
            y = Y, weights = weights, offset = offset, 
            family = tweedie(var.power = p, link.power = link.power),
            intercept = intercept), 
            silent = TRUE)
        skip.obs <- FALSE
        if (class(catch.possible.error) == "try-error") 
            skip.obs <- TRUE        
        if (skip.obs) {
            warning(paste("  Problem near p= ", 
                p, "; this error reported:\n     ", catch.possible.error, 
                " Examine the data and function inputs carefully."))
            mu <- rep(NA, nY)
        } else 
            mu <- fitted(fit.model)        
        if (verbose) 
            cat("* Phi estimation")
        if (skip.obs) {
            if (verbose) 
                cat("; but skipped for this obs\n")
            phi.vec[i] <- NA
        } else {
            if (verbose) 
                cat(" (using optimize): ")
            phi.est <- sum(tweedie.dev(y = Y, mu = mu, power = p))/nY
            low.limit <- min(0.001, phi.est/2)    
            ans <- optimize(f = dtweedie.nlogl, maximum = FALSE, 
                    interval = c(low.limit, 10 * phi.est), power = p, 
                    mu = mu, y = Y)
            phi.vec[i] <- ans$minimum
            if (verbose) 
                  cat(" Done (phi =", phi.vec[i], ")\n")
        }
        if (verbose) {
            cat("* Computing the log-likelihood ")
            cat("(method =", method, "):")
        }
        if (skip.obs) {
            if (verbose) 
                cat(" but skipped for this obs\n")
            L[i] <- NA
        } else {
            if (method == "saddlepoint") 
                L[i] <- dtweedie.logl.saddle(y = Y, mu = mu, 
                  power = p, phi = phi.vec[i], eps = 1/6) else 
                L[i] <- switch(pmatch(method, c("interpolation", 
                        "series", "inversion"), nomatch = 2), 
                        `1` = dtweedie.logl(mu = mu, power = p, phi = phi.vec[i], y = Y), 
                        `2` = sum(log(dtweedie.series(y = Y,  mu = mu, power = p, phi = phi.vec[i]))), 
                        `3` = sum(log(dtweedie.inversion(y = Y, mu = mu, power = p, phi = phi.vec[i]))))
        }
       if (verbose) 
            cat(" L =", L[i], "\n")
    }
    L.max <- max(L)
    p.max <- p.vec[L == L.max]
    phi.max <- phi.vec[L == L.max]
                      
    out <- list(p = p.vec, phi=phi.vec, L=L, 
                p.max = p.max, phi.max = phi.max, L.max = L.max, 
                method = method)
    return(out)  
}
  
# function to take inverse of a matrix using svd 
svd.inv <- function(x){
	sx <- svd(x)
	return(sx$v%*% diag(1/sx$d)%*%t(sx$u))	
}
    
# function to compute the link.power needed in tweedie
make.link.power <- function(link) {
  if (!is.character(link) && !is.numeric(link))
    stop("link.power must be either numeric or character.")
  if (is.character(link)){  
    okLinks <- c("log", "identity", "sqrt","inverse")
    if (link %in% okLinks) 
      switch(link,log=0, identity=1, sqrt=0.5, inverse=-1) else
      stop("invalid link function!")
  } else 
    link  
}

# control options intializer
cpglm.control <- function(init.size=100L,
                       sample.iter=50L,
                       max.size=10000L,
                       max.iter=200,
                       epsilon1=1e-03,
                       epsilon2=1e-04,
                       alpha =0.25,
                       k=5,                       
                       bound.p=c(1.01,1.99),
                       fixed.size=TRUE,   
                       beta.step=10,
                       trace=TRUE,
                       profile.method="inversion",
                       decimal=3){
  if (!is.numeric(init.size) || init.size <= 0)
        stop("value of sample.size should be an integer and >0")
  if (!is.numeric(sample.iter) || sample.iter <= 0)
        stop("value of sample.iter should be an integer and >0")   
  if (!is.numeric(epsilon1) || epsilon1 <= 0) 
        stop("value of 'epsilon1' must be > 0")
  if (!is.numeric(epsilon2) || epsilon2 <= 0) 
        stop("value of 'epsilon2' must be > 0") 
  if (!is.numeric(alpha) || alpha <= 0 || alpha>=1) 
        stop("value of 'alpha' must be between 0 and 1")               
  if (!is.numeric(k) || k <= 0) 
        stop("value of 'k' must be > 0")         
  if (!is.numeric(max.iter) || max.iter <= 0) 
        stop("value of 'maxit' must be > 0")
  if (min(bound.p)<1 || max(bound.p)>2)
        stop("value of 'bound.p' must be between 1 and 2")
  if (!is.numeric(fixed.size) && !is.logical(fixed.size))
        stop("'fixed.size' must be logical or numeric")
  if (!is.numeric(beta.step) || beta.step <= 0) 
        stop("value of 'beta.step' must be greater than 0")          
  if (!is.numeric(trace) && !is.logical(trace))
        stop("'trace' must be logical or numeric")
  if (!is.numeric(decimal) || decimal<=0 )
        stop("'decimal' must be a positive integer")
  if (!(profile.method %in% c("series","inversion",
                              "interpolation","saddlepoint")))
        stop("invalid 'profile.method'")
  bound.p <- sort(bound.p)
  fixed.size <- as.logical(fixed.size)
  trace <- as.logical(trace)
  
    list(init.size=init.size,
         sample.iter=sample.iter,
         max.iter=max.iter,
         epsilon1 = epsilon1,
         epsilon2=epsilon2,
         alpha=alpha,
         k=k,
         fixed.size=fixed.size,
         max.size=max.size,
         bound.p=bound.p,
         beta.step=beta.step,
         trace=trace,
         profile.method=profile.method,
         decimal=decimal)  
}






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
                    
dyn.load("./src/cpglm_em.dll")
dyn.unload("./src/cpglm_em.dll")

# MCEM fit
set.seed(11)
fit1 <- cpglm(RLD~ factor(Zone)*factor(Stock),
	data=fineroot,
  control=list(init.size=5,sample.iter=60,
              max.size=3000,fixed.size=FALSE),
  pstart=1.6)

# profile likelihood         
fit2 <- cpglm(RLD~ factor(Zone)*factor(Stock),
	data=fineroot,method="profile", 
	control=list(decimal=1))      

# compare the two 
summary(fit1)
summary(fit2)

}
