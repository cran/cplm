#######################################################
##           Bayesian compound Poisson GLM           ##
## Author: Wayne Zhang, actuary_zhang@hotmail.com    ##
#######################################################

bcpglm <- function(formula, link = "log", data, inits = NULL,
                   weights, offset, subset, na.action, contrasts = NULL, 
                   n.chains=3, n.iter=2000, n.burnin=floor(n.iter/2),
                   n.thin=max(1, floor(n.chains * (n.iter - n.burnin) / n.sims)),
                   n.sims=1000, n.report=1000, prior.beta.mean, prior.beta.var, 
                   bound.phi=100, bound.p=c(1.01,1.99), method="dtweedie",
                   tune.iter=4000, n.tune=10, tune.weight=0.25,...) {

  call <- match.call()  
  if (missing(data)) 
    data <- environment(formula)   
  mf <- match.call(expand.dots = FALSE)
  # construct model frame
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
  n.beta <- NCOL(X)
  
  ## checking arguments	
  if (!is.null(weights)){
    if (!is.numeric(weights)) 
        stop("'weights' must be a numeric vector")
    if (any(weights <= 0)) 
        stop("negative or zero weights not allowed")
  }
  if (!is.null(offset)) {
    if (length(offset) != NROW(Y)) 
      stop(gettextf("number of 'offset' is %d should 
                    equal %d (number of observations)", 
                length(offset), NROW(Y)), domain = NA)
    }
  if (length(bound.p)!=2)
    stop ("'bound.p' must be of length 2") 
  if (bound.p[1]<1 || bound.p[2]>2)
  	stop("'bound.p' should be between 1 and 2")  
  
  # check initial values
  check.inits(inits,n.chains,bound.p,ncol(X))

  # check counts related inputs
  if (!is.numeric(n.chains) || n.chains <1)
  	stop("'n.chains' must be greater than 1" )
  if (n.burnin>=n.iter || n.sims>=n.iter)
  	stop("'n.burnin' or 'n.sims' should be less than 'n.iter'" )
  if (!missing(prior.beta.mean) && 
  	length(prior.beta.mean)!=n.beta)
  	stop(gettextf("'prior.beta.mean' should be of length %d"),n.beta)  
  if (!missing(prior.beta.mean) && 
  	length(prior.beta.mean)!=n.beta)
  	stop(gettextf("'prior.beta.mean' should be of length %d"),n.beta)
  		  
  # default prior mean and var if missing           
  if (missing(prior.beta.mean))
    prior.beta.mean <- rep(0, n.beta)			
  if (missing(prior.beta.var))
  	prior.beta.var <- rep(10000, n.beta)
  
  # run MCMC 
	bfit <- bcpglm_gibbs(X=X,Y=Y,weights=weights,offset=offset,
                      link.power=link.power, inits=inits, 
                      n.chains=n.chains, n.iter=n.iter, 
                      n.burnin=n.burnin, n.sims=n.sims, n.thin=n.thin, 
                      n.report=n.report, prior.beta.mean=prior.beta.mean, 
                      prior.beta.var=prior.beta.var, intercept=attr(mt, "intercept") > 0L,
                      bound.phi=bound.phi, bound.p=bound.p, n.tune=n.tune,
                      tune.iter=tune.iter, tune.weight=tune.weight, method=method)
  
  # coerce sims.list to mcmc.list from coda
  ans <- new("bcpglm", 
             n.chains=as.integer(n.chains), 
             n.iter=as.integer(n.iter), 
             n.burnin=as.integer(n.burnin),
             n.thin=as.integer(n.thin), 
             n.sims=as.integer(bfit$dims["n.sims"]), 
             sims.list=bfit$sims.list,
             link.power=link.power,
             call=call,
             formula=formula,
             data=data,
             model.frame = mf,
             contrasts=contrasts,
             inits = bfit$inits
             )  
  return(ans)
}


# function to run the MCMC using latent variables
bcpglm_gibbs <- function(X,Y,weights=NULL,offset=NULL,
                      link.power=0, inits=NULL, n.chains, n.iter, 
                      n.burnin, n.sims, n.thin, n.report,
                      prior.beta.mean, prior.beta.var, 
                      bound.phi=100, bound.p=c(1.01,1.99),
                      tune.iter=2000, n.tune=10, tune.weight=0.25,
                      intercept=T, method="dtweedie"){

  X <- as.matrix(X)          
  # get names
  xnames <- dimnames(X)[[2L]]
  ynames <- if (is.matrix(Y)) 
        rownames(Y) else 
        names(Y)
  n.beta <- NCOL(X)
  ygt0 <- as.integer(which(Y>0L))
  # default weights and offsets if NULL    
  n.obs <- NROW(Y)
  if (is.null(weights))     
      weights <- rep.int(1, n.obs)
  if (is.null(offset)) 
        offset <- rep.int(0, n.obs)   
    
  n.keep <- floor((n.iter-n.burnin) / n.thin)
  n.sims <- n.chains * n.keep  
  # dimensions used in simulation
  dims <- list(n.obs= n.obs,
           n.beta=n.beta,
           n.pos= length(ygt0),
           n.chains=as.integer(n.chains), 
           n.iter=as.integer(n.iter), 
           n.burnin=as.integer(n.burnin),
           n.thin=as.integer(n.thin), 
           n.keep=as.integer(n.keep),
           n.sims=as.integer(n.sims),
           n.report=as.integer(n.report),
           tune.iter=as.integer(tune.iter),
           n.tune=as.integer(n.tune))
  dims <- unlist(dims)

  # generating scale matrix in metropolis update 
  fit.start <- cpglm_profile(X=X,Y=Y,weights=weights,offset=offset,
                      link.power=link.power, intercept=intercept)
                  
  ebeta.var <- fit.start$vcov*2.38^2/n.beta
  ephi.var <- 2.38^2 * attr(fit.start$vcov, "phi_p")[1,1]
  ep.var <- 2.38^2 * attr(fit.start$vcov, "phi_p")[2,2]
    
  # generate initial values if necessary
  if (is.null(inits)) {
      pstart <- fit.start$p
      betastart <- as.numeric(fit.start$coefficients)
  	  phistart <- fit.start$phi
		  inits.start <- c(betastart, phistart, pstart)
		  inits <- vector("list",n.chains)
		  cv <- chol(ebeta.var)
	    inits[[1]] <- c(betastart, phistart, pstart)
	    if (n.chains>1){
		    for (i in 2:n.chains)
			    inits[[i]] <- c(betastart + t(cv)%*%rnorm(n.beta),
						        runif(1,phistart/2,1.5*phistart),
						        runif(1,(bound.p[1]+pstart)/2,(bound.p[2]+pstart)/2))
	    }
  }

  # run MCMC   
    # input for the C function 	     
    input <- list(X=X,
               y=as.double(Y),
               ygt0= as.integer(which(Y>0L)-1),
               offset=as.double(offset),
               pWt=as.double(weights),
               mu = double(dims["n.obs"]),
               eta = double(dims["n.obs"]),
               inits = inits,
               beta=as.double(inits[[1]][1:n.beta]),
               phi=as.double(inits[[1]][n.beta+1]),
               p=as.double(inits[[1]][n.beta+2]),
               link.power=as.double(link.power),
               pbeta.mean=as.double(prior.beta.mean),
               pbeta.var=as.double(prior.beta.var),
               bound.phi=as.double(bound.phi),
               bound.p=as.double(bound.p),    
               ebeta.var=as.double(ebeta.var),
               ep.var= as.double(ep.var),
               ephi.var = as.double(ephi.var),               
               dims=dims,
               tune.weight=as.double(tune.weight))
  
  if (method=="dtweedie")
    sims.list<- .Call("bcpglm_gibbs_tw",input) 
  if (method=="latent")
    sims.list<- .Call("bcpglm_gibbs_lat",input)
                  
  sims.list <- lapply(sims.list, function(x){ 
                  dimnames(x) <- list(NULL, c(xnames,"phi","p"))
                  return(x)})  
  # coerce to mcmc object                  
  sims <- lapply(sims.list, as.mcmc)
  sims <- as.mcmc.list(sims)
                         
  out <- list(dims = dims,sims.list=sims, inits = inits)
  return(out)
}

# check initial values of bcpglm input inits
check.inits <- function(inits, n.chains, bound.p, n.beta){
  if (!is.null(inits)){
  	if (!is.list(inits))
  		stop("initial values must be in a list in 'inits'") else 
  	if (any(sapply(inits,length)!=n.beta+2))
  		stop(gettextf("each set of initial values must be of length %d"), 
  				  n.beta+2) else 
  	if (any(sapply(inits,"[", n.beta+1)<0))
  		stop(gettextf("dispersion parameter (element %d) must be positive", 
  			  	n.beta+1)) else      	
  	if (any(sapply(inits,"[", n.beta+2)<bound.p[1]) ||
  			any(sapply(inits,"[", n.beta+2)>bound.p[2]))
  		stop(gettextf("index parameter (element %d) out of specified bound", 
  				  n.beta+2))
    }  	
}

