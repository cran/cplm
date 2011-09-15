#######################################################
##           Bayesian compound Poisson GLM           ##
## Author: Wayne Zhang, actuary_zhang@hotmail.com    ##
#######################################################

bcpglm <- function(formula, link = "log", data, inits = NULL,
                   weights, offset, subset, na.action, contrasts = NULL, 
                   n.chains=3, n.iter=2000, n.burnin=floor(n.iter/2),
                   n.thin=max(1, floor(n.chains * (n.iter - n.burnin) / n.sims)),
                   n.sims=1000, n.report=1000, prior.beta.mean, prior.beta.var, 
                   phi.shape=0.001, phi.scale=0.001, bound.p=c(1.01,1.99),...) {
  
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
	bfit <- bcpglm_gibbs_lat(X=X,Y=Y,weights=weights,offset=offset,
                      link.power=link.power, inits=inits, 
                      n.chains=n.chains, n.iter=n.iter, 
                      n.burnin=n.burnin, n.sims=n.sims, n.thin=n.thin, 
                      n.report=n.report, prior.beta.mean=prior.beta.mean, 
                      prior.beta.var=prior.beta.var, phi.shape=phi.shape, 
                      phi.scale=phi.scale, bound.p=bound.p)
  
  # coerce sims.list to mcmc.list from coda
  sims <- lapply(bfit$sims.list, as.mcmc)
  sims <- as.mcmc.list(sims)
  ans <- new("bcpglm", 
             n.chains=as.integer(n.chains), 
             n.iter=as.integer(n.iter), 
             n.burnin=as.integer(n.burnin),
             n.thin=as.integer(n.thin), 
             n.sims=as.integer(bfit$n.sims), 
             sims.list=sims,
             summary = summary(sims),
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
bcpglm_gibbs_lat <- function(X,Y,weights=NULL,offset=NULL,
                      link.power=0, inits=NULL, n.chains, n.iter, 
                      n.burnin, n.sims, n.thin, n.report,
                      prior.beta.mean, prior.beta.var, phi.shape=0.001, 
                      phi.scale=0.001, bound.p=c(1.01,1.99)){

    X <- as.matrix(X)          
    # get names
    xnames <- dimnames(X)[[2L]]
    ynames <- if (is.matrix(Y)) 
        rownames(Y) else 
        names(Y)
    n.beta <- NCOL(X)
    # default weights and offsets if NULL    
    nobs <- NROW(Y)
    if (is.null(weights))     
      weights <- rep.int(1, nobs)
    if (is.null(offset)) 
        offset <- rep.int(0, nobs)   
    
    # generating scale matrix in metropolis update of beta
    pstart <- if (!is.null(inits))
    			inits[[1]][n.beta+2] else 
    			sum(bound.p)/2
    fit.start <- glm(Y~-1+X,weights=weights,offset=offset,
                  family=tweedie(var.power=pstart,
                                 link.power=link.power))
    ebeta.var <- vcov(fit.start)
    
    # generate initial values if necessary
    if (is.null(inits)) {                             
	      betastart <- as.numeric(fit.start$coefficients)
    	  phistart <- sum(residuals(fit.start,"pearson")^2)/
          			  df.residual(fit.start)
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
    
    # start MCMC 
    sims.list <- vector("list",n.chains)
    n.keep <- floor((n.iter-n.burnin) / n.thin)
    n.sims <- n.chains * n.keep
    	if (n.report>0){
    		cat(paste(rep("-",50), collapse=""))
    		cat("\n")
	    	cat("Markov Chain Monte Carlo starts...\n")
    		cat(paste(rep("-",50), collapse=""))	   
    		cat("\n") 	
    	}
    for (i in 1:n.chains){
    	if (n.report>1)
	    	cat("Start Markov chain", i, ":\n")    	
	    sims.list[[i]] <- .Call("bcpglm_gibbs_lat",
                 X=as.double(X),
                 Y=as.double(Y),
                 ygt0= as.integer(which(Y>0L)),
                 offset=as.double(offset),
                 weights=as.double(weights),
                 beta=as.double(inits[[i]][1:n.beta]),
                 phi=as.double(inits[[i]][n.beta+1]),
                 p=as.double(inits[[i]][n.beta+2]),
                 link.power=as.double(link.power),
                 pbeta.mean=as.double(prior.beta.mean),
                 pbeta.var=as.double(prior.beta.var),
                 pphi.shape =as.double(phi.shape),
                 pphi.scale =as.double(phi.scale),
                 bound.p=as.double(bound.p),                 
                 ebeta.var=as.double(ebeta.var),
                 n.iter=as.integer(n.iter),
                 n.burnin=as.integer(n.burnin),
                 n.thin=as.integer(n.thin),
                 n.keep=as.integer(n.keep),
                 n.report=as.integer(n.report))
     dimnames(sims.list[[i]]) <- list(NULL, c(xnames,"phi","p"))
     if (n.report>1){
        cat(paste(rep("-",50), collapse=""))
        cat("\n")
      }
    }
    if (n.report>0){
    	cat("End of Markov Chain Monte Carlo simulation.\n")
    	cat(paste(rep("-",50), collapse=""))	    	
    	cat("\n")
    }     
    out <- list(n.chains=n.chains, n.iter=n.iter, 
                n.burnin=n.burnin, n.thin=n.thin, 
                n.sims=n.sims, sims.list=sims.list,
                inits = inits)
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
    