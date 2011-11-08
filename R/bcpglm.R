#######################################################
##           Bayesian compound Poisson GLM           ##
## Author: Wayne Zhang, actuary_zhang@hotmail.com    ##
#######################################################

bcpglm <- function(formula, link = "log", data, inits = NULL,
                   weights, offset, subset, na.action, contrasts = NULL, 
                   n.chains=3, n.iter=2000, n.burnin=floor(n.iter/2),
                   n.thin=max(1, floor(n.chains * (n.iter - n.burnin) / n.sims)),
                   n.sims=1000, n.report=1000, prior.beta.mean=NULL, 
                   prior.beta.var=NULL, bound.phi=100, bound.p=c(1.01,1.99), 
                   method="dtweedie", tune.iter=4000, n.tune=10, tune.weight=0.25,...) {

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
  n.obs <- NROW(X)
  n.beta <- NCOL(X)
  
  # check arguments
  check.args.bcplm(call,n.beta, n.chains)
  # check initial values
  if (!is.null(inits))
    check.inits.bcpglm(inits,n.beta, n.chains)  
  		  
  # default prior mean and var if missing           
  if (is.null(prior.beta.mean))
    prior.beta.mean <- rep(0, n.beta)			
  if (is.null(prior.beta.var))
  	prior.beta.var <- rep(10000, n.beta)
  
  # default weights and offset
  if (is.null(weights))     
      weights <- rep(1, n.obs)
  if (is.null(offset)) 
        offset <- rep(0, n.obs)   
     
  # dimensions used in simulation  
  n.keep <- floor((n.iter-n.burnin) / n.thin)
  n.sims <- n.chains * n.keep  
  dims <- list(n.obs= as.integer(n.obs),
           n.beta=as.integer(n.beta),
           n.pos= as.integer(sum(Y>0)),     
           n.term = as.integer(0),
           n.u = as.integer(0),
           n.all = as.integer(n.beta+2),
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
 
  C <- 2.38 * 2.38             
  # proposal covariance matrix              
  ebeta.var <- diag(1, nrow=n.beta) * C/n.beta
  ephi.var <- ep.var <- C 
     
  # generate initial values if necessary
  if (is.null(inits)) {
       # generating starting values and scale matrix in metropolis update 
      fit.start <- cpglm_profile(X=X,Y=Y,weights=weights,offset=offset,
                  link.power=link.power, intercept=attr(mt, "intercept") > 0L)
                      
      pstart <- fit.start$p
      betastart <- as.numeric(fit.start$coefficients)
      phistart <- fit.start$phi
		  inits.start <- c(betastart, phistart, pstart)
		  inits <- vector("list",n.chains)
		  inits[[1]] <- c(betastart, phistart, pstart)
	    if (n.chains>1){
		    for (i in 2:n.chains)
			    inits[[i]] <- c(betastart + rnorm(n.beta,0,0.5),
						        runif(1,phistart/2,1.5*phistart),
						        runif(1,(bound.p[1]+pstart)/2,(bound.p[2]+pstart)/2))
	    }
      # update proposal covariance matrix
      ebeta.var <- fit.start$vcov*C/n.beta
      ephi.var <- C * attr(fit.start$vcov, "phi_p")[1,1]
      ep.var <- C * attr(fit.start$vcov, "phi_p")[2,2]
  } else {
    inits <- lapply(inits, function(x) c(x$beta, x$phi, x$p ))
  }

  # run MCMC   
    # input for the C function 	     
    input <- list(X=as.double(X),
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
  
  # get names
  sims.list <- lapply(sims.list, function(x){ 
                  dimnames(x) <- list(NULL, c(dimnames(X)[[2L]],"phi","p"))
                  return(x)})  
  # coerce to mcmc object                  
  sims <- lapply(sims.list, as.mcmc)
  sims <- as.mcmc.list(sims)
   
  # coerce sims.list to mcmc.list from coda
  ans <- new("bcpglm", 
             n.chains=as.integer(n.chains), 
             n.iter=as.integer(n.iter), 
             n.burnin=as.integer(n.burnin),
             n.thin=as.integer(n.thin), 
             n.sims=as.integer(dims["n.sims"]), 
             sims.list=sims,
             link.power=link.power,
             call=call,
             formula=formula,
             model.frame = mf,
             contrasts=contrasts,
             inits = inits)  
  return(ans)
}


