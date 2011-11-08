
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

check.args.bcplm <- function(call,n.beta, n.chains){
  # check counts related inputs
  if (!is.null(call$n.chains) && (!is.numeric(call$n.chains) 
      || call$n.chains <1))
    stop("'n.chains' must be greater than 1" )
  if (!is.null(call$n.burnin) && !is.null(call$n.iter) && 
      call$n.burnin>=call$n.iter)
  	stop("'n.burnin' should be less than 'n.iter'" )
  if (!is.null(call$prior.beta.mean) && 
  	length(call$prior.beta.mean)!=n.beta)
  	stop(gettextf("'prior.beta.mean' should be of length %d"),n.beta)  
  if (!is.null(call$prior.beta.mean) && 
  	length(call$prior.beta.mean)!=n.beta)
  	stop(gettextf("'prior.beta.mean' should be of length %d"),n.beta)
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
    
  if (length(inits$phi)>1 || inits$phi<=0) 
    stop("'phi' in 'inits' should be of length 1 and greater than 0")
  if (length(inits$p)>1 || inits$p<=1 || inits$p>=2) 
    stop("'p' in 'inits' should be of length 1 and between 1 and 2")
}

# check initial values in cpglmm
check.inits.cpglmm <- function(inits, n.beta, n.term){
  if (!("Sigma" %in% names(inits)))
    stop("the 'Sigma' component in 'inits' is missing") 
  if (length(inits$Sigma)!=n.term) 
    stop(gettextf("'Sigma' in 'inits' should be of length %d", n.term))
}

# check initial values in bcpglm
check.inits.bcpglm <- function(inits, n.beta, n.chains){
  if (length(inits)!=n.chains)
    stop(gettextf("'inits' should be of length %d", n.chains))
  invisible(lapply(inits, check.inits.cpglm))
}

# check initial values in bcpglmm
check.inits.bcpglmm <- function(inits, n.beta, n.term, n.chains){
  if (length(inits)!=n.chains)
    stop(gettextf("'inits' should be of length %d", n.chains))
  invisible(lapply(inits, check.inits.cpglmm))
}  


###########################################################
# parse and default prior info for the variance component
###########################################################

# inverse gamma
igamma <- function(scale=0.01,shape=0.01){
  return(list(igamma.scale=scale, igamma.shape=shape))
}  

# inverse wishart
iwish <- function(df=3, scale=diag(1,df)){
  return(list(iwish.df=df, iwish.scale=scale))
}  

# default prior info for the variance component
prior.Sigma.default <- function(Sigma){
  lapply(Sigma, function(x){
        nc <- ncol(x)
        if (nc==1) 
          igamma(scale=0.001,shape=0.001) else
          iwish(df=as.double(nc))
  })
}
