###########################################################
# check arguments 
###########################################################

check.args.cplm <- function(call,n.obs){
  ## checking arguments  
  if (!is.null(call$weights)){
    if (!is.numeric(call$weights)) 
        stop("'weights' must be a numeric vector")
    if (any(call$weights <= 0)) 
        stop("negative or zero weights not allowed")
  }
  if (!is.null(call$offset)) {
    if (length(call$offset) != n.obs) 
      stop(gettextf("number of 'offset' is %d should 
                    equal %d (number of observations)", 
                length(call$offset), n.obs), domain = NA)
    }
}

check.args.bcplm <- function(call, n.beta, n.chains){
  # check counts related inputs
  if (!is.null(call$n.chains) && (!is.numeric(call$n.chains) 
      || call$n.chains < 1))
    stop("'n.chains' must be greater than 1" )
  if (!is.null(call$n.burnin) && !is.null(call$n.iter) && 
      call$n.burnin >= call$n.iter)
  	stop("'n.burnin' should be less than 'n.iter'" )
  if (!is.null(call$prior.beta.mean) && 
  	length(call$prior.beta.mean) != n.beta)
  	stop(gettextf("'prior.beta.mean' should be of length %d"), n.beta)  
  if (!is.null(call$prior.beta.mean) && 
  	length(call$prior.beta.mean) != n.beta)
  	stop(gettextf("'prior.beta.mean' should be of length %d"), n.beta)
}

###########################################################
# Check initial values
###########################################################

# check initial values in cpglm
check.inits.cpglm <- function(inits, n.beta){
  if (!("beta" %in% names(inits)))
    stop("the 'beta' component in 'inits' is missing") 
  if (!("phi" %in% names(inits)))
    stop("the 'phi' component in 'inits' is missing")
  if (!("p" %in% names(inits)))
    stop("the 'p' component in 'inits' is missing")
  if (length(inits$beta) != n.beta)
    stop(gettextf("number of 'beta' in 'inits' is %d, but should 
                    equal %d (number of mean parameters)", 
                length(inits$beta), n.beta, domain = NA))
    
  if (length(inits$phi) > 1 || inits$phi <= 0) 
    stop("'phi' in 'inits' should be of length 1 and greater than 0")
  if (length(inits$p) > 1 || inits$p <= 1 || inits$p >= 2) 
    stop("'p' in 'inits' should be of length 1 and between 1 and 2")
}

# check initial values in cpglmm
check.inits.cpglmm <- function(inits, n.beta, n.term){
  if (!("Sigma" %in% names(inits)))
    stop("the 'Sigma' component in 'inits' is missing") 
  if (length(inits$Sigma) != n.term) 
    stop(gettextf("'Sigma' in 'inits' should be of length %d", n.term))
}

# check initial values in bcpglm
check.inits.bcpglm <- function(inits, n.beta, n.chains){
  if (length(inits) != n.chains)
    stop(gettextf("'inits' should be of length %d", n.chains))
  invisible(lapply(inits, check.inits.cpglm, n.beta))
}

# check initial values in bcpglmm
check.inits.bcpglmm <- function(inits, n.beta, n.term, n.chains){
  if (length(inits) != n.chains)
    stop(gettextf("'inits' should be of length %d", n.chains))
  invisible(lapply(inits, function(x) check.inits.cpglmm(x, n.beta, n.term)))
}  


###########################################################
# default control options   
###########################################################
# control options intializer
cpglm.control <- function(bound.p = c(1.01,1.99),
                        trace = 0,
                        max.iter = 200L){ 
  if (min(bound.p)<1 || max(bound.p)>2)
        stop("value of 'bound.p' must be between 1 and 2")        
  if (!is.numeric(trace) && !is.logical(trace))
        stop("'trace' must be logical or numeric")
  if (!is.numeric(max.iter) || max.iter < 0)
        stop("'max.iter' must be greater than 0")
  list(bound.p = sort(bound.p),         
        trace = as.integer(trace),
        max.iter = as.integer(max.iter))  
}


# update control parameters in cpglmm 
cpglmm.control <- function(max.iter = 300L,
                       max.fun = 20000L,               
                       bound.p = c(1.01, 1.99),
                       trace = 0){         
  if (!is.numeric(max.iter) || max.iter <= 0) 
        stop("value of 'max.iter' must be > 0")
  if (!is.numeric(max.fun) || max.fun <= 0) 
        stop("value of 'max.fun' must be > 0")
  if (!is.numeric(bound.p) || length(bound.p) != 2)
        stop("'bound.p' must be of length 2")
  if (min(bound.p) < 1 || max(bound.p) > 2)
        stop("invalid bounds in 'bound.p'")          
  if (!is.numeric(trace) && !is.logical(trace))
        stop("'trace' must be logical or numeric")
  
  list(max.iter = as.integer(max.iter),
       max.fun = as.integer(max.fun),
       bound.p = as.numeric(sort(bound.p)),
       trace = as.integer(trace))  
}

###########################################################
# parse and default prior info for the variance component
###########################################################

# inverse gamma
igamma <- function(scale = 0.01, shape = 0.01){
  return(list(igamma.scale = scale, igamma.shape = shape))
}  

# inverse wishart
iwish <- function(df = 3, scale = diag(1, df)){
  return(list(iwish.df = df, iwish.scale = scale))
}  

# default prior info for the variance component
prior.Sigma.default <- function(Sigma){
  lapply(Sigma, function(x){
        nc <- ncol(x)
        if (nc == 1) 
          igamma(scale = 0.001, shape = 0.001) else
          iwish(df = as.double(nc))
  })
}

###########################################################
# numerical derivatives  
###########################################################

# function to compute gradient
grad <- function(parm, fun, ...){
  n <- length(parm)
  eps <- 0.001
  gd <- rep(NA, n)
  for (i in 1:n){
    parm[i] <- parm[i] - eps
    g1 <- fun(parm, ...)
    parm[i] <- parm[i] + 2 * eps
    g2 <- fun(parm, ...)
    gd[i] <- (g2 - g1) / (2 * eps)
  }
  return(gd)
}

# function to compute hessian
hess <- function(parm, fun, ...){
  n <- length(parm)
  eps <- 0.001
  hn <- matrix(0, n, n)
  for (i in 1:n){
    parm[i] <- parm[i] - eps
    g1 <- grad(parm, fun, ...)
    parm[i] <- parm[i] + 2 * eps
    g2 <- grad(parm, fun, ...)
    hn[i,] <- (g2 - g1) / ( 2 * eps)
  }
  return(hn)  
}


###########################################################
# glm related   
###########################################################

# construct model frame in cpglm   
cpglm.mf <- function(mf, contrast){  
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
    n.obs <- nrow(X)
    if (is.null(weights))
      weights <- rep(1, n.obs)    
    if (is.null(offset))
      offset <- rep(0, n.obs)
    return (list(mf = mf, wts = weights, off = offset,
                 Y = Y, X = X))
  }

# fit a Tweedie glm given a model frame 
cpglm.fit <- function(fr, p = 1.5, link.power = 0) {
  fm <- tweedie(var.power = p,link.power = link.power)
  int <- attr(attr(fr$mf,"terms"), "intercept") > 0L
  glm.fit(fr$X, fr$Y, weights = fr$wts, offset = fr$off,
          family = fm, intercept = int)
}
  
# compute fitted values of for bigglm  
fitted.bigglm <- function(object, data, ...){
  # get chunks of data
  tt <- terms(object)
  n <- object$n
  beta <- coef(object)
  cursor <- 0
  eta <- offset <- pwts <- c()
  datafun <- function(){
    if (cursor >= n)
        return(NULL)
    start <- cursor + 1
    cursor <<- cursor + min(object$call$chunksize, n - cursor)
    data[start:cursor, ]
  }
  # get stats for each chunk
  while(!is.null(chunk <- datafun())){
    mf <- model.frame(tt, chunk)
    mm <- model.matrix(tt, mf)
    if(is.null(off <- model.offset(mf))) 
      off <- rep(0, nrow(mm))  
    if (!is.null(object$weights))
      w <- model.frame(object$weights, chunk)[[1]] else 
      w <- rep(1, nrow(mm))
    eta <- c(eta, mm %*% beta + off)    
    offset <- c(offset, off)
    pwts <- c(pwts, w)
  }
  # compute stats to be returned
  mu <- object$family$linkinv(eta)
  dmu <- object$family$mu.eta(eta)
  wts <- pwts * dmu * dmu / (object$family$variance(mu))
  y <- eval(object$call$formula[[2]], data)
  res <- (y - mu) / dmu
  list(linear.predictors = eta, 
        fitted.values = mu,
        offset = offset,
        prior.weights = pwts,
        weights = wts,
        residuals = res )
}


###########################################################
# general utility functions  
###########################################################

# function to compute log density 
dtweedie.nlogl <- function(y, mu, phi,power) {
    ans <- -2 * sum(log(dtweedie(y = y, mu = mu, phi = phi, power = power)))
    if (is.infinite(ans)) {
        ans <- sum(tweedie.dev(y = y, mu = mu, power = power)) / length(y)
    }    
    #attr(ans, "gradient") <- dtweedie.dldphi(y = y, mu = mu, 
    #    phi = phi, power = power)
    ans
}
    
  
# function to take inverse of a matrix using svd 
svd.inv <- function(x){
  sx <- svd(x)
  return(sx$v %*% diag(1 / sx$d) %*% t(sx$u))	
}
    
# function to compute the link.power needed in tweedie
make.link.power <- function(link) {
  if (!is.character(link) && !is.numeric(link))
    stop("link.power must be either numeric or character.")
  if (is.character(link)){  
    okLinks <- c("log", "identity", "sqrt","inverse")
    if (link %in% okLinks) 
      switch(link, log = 0, identity = 1, sqrt = 0.5, inverse = -1) else
      stop("invalid link function!")
  } else 
    link  
}

