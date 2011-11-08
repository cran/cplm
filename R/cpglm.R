#######################################################
##             compound Poisson GLM                  ##
## Author: Wayne Zhang, actuary_zhang@hotmail.com    ##
#######################################################

cpglm <- function(formula, link = "log", data, weights, offset, 
                  subset, na.action, inits=NULL, contrasts = NULL, 
                  control=list(), method ="profile", ...) {

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
  n.obs <- NROW(X)
  n.beta <- NCOL(X)
  
  # check arguments 
  check.args.cplm(call,n.obs)
  
  if (method=="MCEM")
    cpfit <- cpglm_em(X,Y,weights=weights,offset=offset,
                     link.power=link.power,inits=inits,
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
             call=call,
             formula=formula,           
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
             link.power=link.power,
             model.frame = mf,
             na.action = attr(mf,"na.action"),
             offset = cpfit$offset,
             prior.weights =cpfit$prior.weights,
             y = Y,
             inits = inits)  
  return(ans)
}


# function to run the MCEM 
cpglm_em <- function(X,Y,weights=NULL,offset=NULL,
                      link.power=0, inits=NULL,
                      intercept = TRUE,
                      control=list()){
    # set control options                        
    control <- do.call("cpglm.control", control)                   

    X <- as.matrix(X)          
    # get names
    xnames <- dimnames(X)[[2L]]
    ynames <- if (is.matrix(Y)) 
        rownames(Y) else 
        names(Y)
    
    # default weights and offsets if NULL    
    n.obs <- NROW(Y)
    if (is.null(weights))     
      weights <- rep.int(1, n.obs)
    if (is.null(offset)) 
        offset <- rep.int(0, n.obs)          
    # generating starting values if necessary
    if (!is.null(inits)){
      check.inits.cpglm(inits, NCOL(X))
      betastart <- inits$beta
      phistart <- inits$phi
      pstart <- inits$p
    } else {
      pstart <- 1.5
      fit.start <- glm(Y~-1+X,weights=weights,offset=offset,
                  family=tweedie(var.power=pstart,
                                 link.power=link.power))
      betastart <- as.numeric(fit.start$coefficients)
      phistart <- sum(residuals(fit.start,"pearson")^2)/
          df.residual(fit.start)
    }
    
    out <- .Call("cpglm_em",
                 X=as.double(X),
                 Y=as.double(Y),
                 ygt0= as.integer(which(Y>0L)-1),
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
                phi = out$phi, power = out$p))) + 2*ncol(X)
    out$prior.weights <- weights
    out$offset <- offset 
    out$converged <- as.logical(out$converged)                       
    out$control <- control
    names(out$coefficients) <- xnames
    names(out$residuals) <- names(out$fitted.values) <-
      names(out$linear.predictors) <- names(out$weights) <- ynames
    return(out)        
}   

# function to implement the  profile likelihood approach 
cpglm_profile <- function(X,Y,weights=NULL,offset=NULL,
                      link.power=0, intercept=TRUE, 
                      contrasts, control=list()){
  control <- do.call("cpglm.control", control)
   
  # profiled likelihood 
  llik_profile <- function(parm){
    phi <- exp(parm[1])
    p <- parm[2]
    fit2 <- glm.fit(X,Y,weights=weights,offset=offset,
                  family=tweedie(var.power=p,
                                 link.power=link.power),
                  intercept=intercept) 
    -2*sum(log(dtweedie.series(Y,p,fit2$fitted.values,phi)))    
  }
  
  # generate starting values for phi
  pstart <- 1.5
  fit <- glm.fit(X,Y,weights=weights,offset=offset,
                  family=tweedie(var.power=pstart,
                                 link.power=link.power),
                  intercept=intercept)  
  mu <- fit$fitted.values
  phistart <- sum((Y-mu)^2/mu^pstart)/fit$df.residual
  parm <- c(log(phistart),pstart)
  
  # optimize the profiled loglikelihood
  opt_ans <- optim(parm,llik_profile,gr=NULL,method="L-BFGS-B",
                      lower=c(-Inf,control$bound.p[1]),
                      upper=c(Inf,control$bound.p[2]),
                      control=list(trace=control$trace))
  p.max <- opt_ans$par[2]
  phi.max <- exp(opt_ans$par[1])

  # fit glm using the optimized index parameter
  fit <- glm.fit(X,Y,weights=weights,offset=offset,
                  family=tweedie(var.power=p.max,
                                 link.power=link.power),
                  intercept=intercept)
  class(fit) <- "glm" 
  
  # compute vcov for p and phi  
  llik_profile2 <- function(parm){
    phi <- parm[1]
    p <- parm[2]
    fit2 <- glm.fit(X,Y,weights=weights,offset=offset,
                  family=tweedie(var.power=p,
                                 link.power=link.power),
                  intercept=intercept) 
    -sum(log(dtweedie.series(Y,p,fit2$fitted.values,phi)))    
  }
  pm <- c(phi.max,p.max) 
  hs <- hess(pm,llik_profile2)
  dimnames(hs) <- list(c("phi","p"),c("phi","p"))  
  vc <- vcov(fit)
  attr(vc,"phi_p") <- solve(hs)
    
  # return results
  out <- c(list(
             deviance=sum(tweedie.dev(Y, fitted(fit),p.max)),
             aic=dtweedie.nlogl(Y,fitted(fit),phi.max,p.max)+2*fit$rank,
             control=control,
             p=p.max,
             phi=phi.max,             
             theta=c(fit$cofficients,phi.max,p.max),
             theta.all=matrix(c(fit$cofficients,phi.max,p.max),
                              nrow=1),
             vcov=vc,
             offset=offset),
             fit[c("coefficients","residuals","fitted.values",
                    "linear.predictors","iter","weights",
                    "prior.weights","df.residual","converged")])  
  return(out)  
}               


# function to compute log density 
dtweedie.nlogl <- function(y, mu, phi,power) {
    ans <- -2 * sum(log(dtweedie(y = y, mu = mu, phi = phi, power = power)))
    if (is.infinite(ans)) {
        ans <- sum(tweedie.dev(y = y, mu = mu, power = power))/length(y)
    }    
    #attr(ans, "gradient") <- dtweedie.dldphi(y = y, mu = mu, 
    #    phi = phi, power = power)
    ans
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
                       k=5,                       
                       bound.p=c(1.01,1.99),
                       fixed.size=TRUE,   
                       beta.step=10,
                       trace=0){
  if (!is.numeric(init.size) || init.size <= 0)
        stop("value of sample.size should be an integer and >0")
  if (!is.numeric(sample.iter) || sample.iter <= 0)
        stop("value of sample.iter should be an integer and >0")   
  if (!is.numeric(epsilon1) || epsilon1 <= 0) 
        stop("value of 'epsilon1' must be > 0")
  if (!is.numeric(epsilon2) || epsilon2 <= 0) 
        stop("value of 'epsilon2' must be > 0") 
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
  bound.p <- sort(bound.p)
  fixed.size <- as.logical(fixed.size)
  trace <- as.integer(trace)
  
    list(init.size=init.size,
         sample.iter=sample.iter,
         max.iter=max.iter,
         epsilon1 = epsilon1,
         epsilon2=epsilon2,
         k=k,
         fixed.size=fixed.size,
         max.size=max.size,
         bound.p=bound.p,
         beta.step=beta.step,
         trace=trace)  
}

# function to compute gradient
grad <- function(parm, fun){
  n <- length(parm)
  eps <- 0.001
  gd <- rep(NA,n)
  for (i in 1:n){
    parm[i] <- parm[i]- eps
    g1 <- fun(parm)
    parm[i] <- parm[i]+2*eps
    g2 <- fun(parm)
    gd[i] <- (g2-g1)/(2*eps)
  }
  return(gd)
}

# function to compute hessian
hess <- function(parm, fun){
  n <- length(parm)
  eps <- 0.001
  hn <- matrix(0,n,n)
  for (i in 1:n){
    parm[i] <- parm[i]- eps
    g1 <- grad(parm,fun)
    parm[i] <- parm[i]+2*eps
    g2 <- grad(parm,fun)
    hn[i,] <- (g2-g1)/(2*eps)
  }
  return(hn)  
}



