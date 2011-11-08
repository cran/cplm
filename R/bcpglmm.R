#######################################################
##    Bayesian Compound Poisson GLMM                 ##
## Author: Wayne Zhang, actuary_zhang@hotmail.com    ##
#######################################################

bcpglmm <- function(formula, link = "log", data, inits = NULL,
                   weights, offset, subset, na.action, contrasts = NULL, 
                   n.chains=3, n.iter=2000, n.burnin=floor(n.iter/2),
                   n.thin=max(1, floor(n.chains * (n.iter - n.burnin) / n.sims)),
                   n.sims=1000, n.report=1000, prior.beta.mean=NULL, prior.beta.var=NULL, 
                   bound.phi=100, bound.p=c(1.01,1.99),  prior.Sigma = NULL,
                   tune.iter=4000, n.tune=10, tune.weight=0.25,...) {
      
   
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
  
  # get factor list and dims                  
  FL <- lme4:::lmerFactorList(formula, fr, 0L, 0L)
  dm <- lme4:::mkZt(FL, NULL)
  n.obs <- unname(dm$dd['n'])
  n.beta <- unname(dm$dd['p'])
  
  # check initial values
  if (!is.null(inits))
    check.inits.bcpglmm(inits=inits,n.beta=ncol(fr$X),
                       n.term=as.numeric(dm$dd['nt']), n.chains=n.chains)
  M1 <- length(levels(dm$flist[[1]]))
  if (M1 >= n.obs) {
        msg1 <- "Number of levels of a grouping factor for the random effects\n"
        msg3 <- "n, the number of observations"
        if (dm$dd["useSc"]) 
            stop(msg1, "must be less than ", msg3)
        else if (M1 == n.obs) 
            message(msg1, "is *equal* to ", msg3)
    }
        
  # default prior info
  if (is.null(prior.beta.mean))
    prior.beta.mean <- rep(0, n.beta)
  if (is.null(prior.beta.var))
    prior.beta.var <- rep(10000, n.beta)
  if (is.null(prior.Sigma))
    prior.Sigma <- prior.Sigma.default(dm$ST)
  
  # default weights and offsets if NULL    
  if (is.null(wts))     
      wts <- rep.int(1, n.obs)
  if (is.null(offset)) 
        offset <- rep.int(0, n.obs)
  
  # dimensions used in simulation
  nc <- unlist(lapply(dm$ST, ncol))
  nlev <- diff(dm$Gp)/nc  
  n.keep <- floor((n.iter-n.burnin) / n.thin)
  n.sims <- n.chains * n.keep    
  dims <- list(n.obs= n.obs,
           n.beta=as.integer(unname(n.beta)),           
           n.pos= as.integer(sum(fr$Y>0)),
           n.term = as.integer(dm$dd['nt']),
           n.u = as.integer(dm$dd['q']),
           n.all = as.integer(dm$dd['p'] + dm$dd['q'] + sum(nc^2)+2),
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
  ebeta.var <- diag(1, nrow=n.beta) * C /n.beta
  eu.var <- diag(1, nrow=dm$dd['q']) * C/dm$dd['q']
  ephi.var <- ep.var <-  C 
                   
  # generate initial values if necessary
  if (is.null(inits)) {
    fit.start <- cplm:::cpglmm(formula=formula, link = link, data=data)
    pstart <- fit.start@p
    betastart <- as.numeric(fixef(fit.start))
    ustart <- as.double(unlist(lapply(ranef(fit.start), as.vector)))
    phistart <- fit.start@phi
		Sigmastart <- lapply(fit.start@ST,function(x) x%*%t(x))
  	Sigmastartv <- unlist(lapply(Sigmastart, as.numeric))
    inits.start <- c(betastart, phistart, pstart, ustart, Sigmastartv)
    inits <- vector("list",n.chains)
	  inits[[1]] <- c(betastart, phistart, pstart, ustart, Sigmastartv)
	  if (n.chains>1){
		    for (i in 2:n.chains)
			    inits[[i]] <- c(as.numeric(betastart + rnorm(n.beta,0,0.5)),
						        runif(1,phistart/2,1.5*phistart),
						        runif(1,(bound.p[1]+pstart)/2,(bound.p[2]+pstart)/2),
                    ustart+rnorm(length(ustart),0,0.5) ,
                    Sigmastartv)
	   }
    #update proposal
    ebeta.var <- as.matrix(vcov(fit.start)*C/n.beta)
    eu.var <- lapply(1:dm$dd['nt'], function(x) 
              kronecker(Sigmastart[[x]],diag(1,nrow=nlev[x])))
    eu.var <- as.matrix(do.call(bdiag, eu.var))*C/dm$dd['q']
  } else{
    inits <- lapply(inits, function(x) c(x$beta, x$phi, x$p, x$u,
                                      unlist(lapply(x$Sigma, as.numeric))))
  }

  # run MCMC   
    # input for the C function 	     
    input <- list(X=fr$X,
               y=as.double(fr$Y),
               Zt = dm$Zt, 
               ygt0= as.integer(which(fr$Y>0L)-1),
               offset=as.double(offset),
               pWt=as.double(wts),
               mu = double(dims["n.obs"]),
               eta = double(dims["n.obs"]),
               inits = inits,
               beta=as.double(betastart),
               u= as.double(ustart),
               phi=as.double(phistart),
               p=as.double(pstart),
               link.power=as.double(link.power),
               pbeta.mean=as.double(prior.beta.mean),
               pbeta.var=as.double(prior.beta.var),
               pSigma = prior.Sigma,
               bound.p=as.double(bound.p),
               bound.phi=as.double(bound.phi),    
               ebeta.var=as.double(ebeta.var),
               eu.var = as.double(eu.var),    
               ep.var= as.double(ep.var),
               ephi.var = as.double(ephi.var),               
               dims=dims,
               tune.weight=as.double(tune.weight),
               Gp = unname(dm$Gp),
               Sigma = Sigmastart,
               ncol = as.integer(nc), 
               nlev = as.integer(nlev))

  sims.list <- .Call("bcpglmm_gibbs_tw", input) 
  xnames <- names(fr$fixef)
  unames <- paste("u", 1:dm$dd['q'],sep="")
  snames <- sapply(1:length(dm$ST), function(x){
              tm <- apply(expand.grid(1:nc[x],1:nc[x]),1,
                          paste, collapse=",")
              paste("Sigma",x,"[",tm,"]",sep="")
  })
  snames <- unlist(snames)
  sims.list <- lapply(sims.list, function(x){ 
                  dimnames(x) <- list(NULL, 
                      c(xnames,"phi","p",unames,snames))
                  return(x)}) 
  # coerce to mcmc object                  
  sims <- lapply(sims.list, as.mcmc)
  sims <- as.mcmc.list(sims)
 
  # coerce sims.list to mcmc.list from coda
  ans <- new("bcpglmm", 
             n.chains=as.integer(n.chains), 
             n.iter=as.integer(n.iter), 
             n.burnin=as.integer(n.burnin),
             n.thin=as.integer(n.thin), 
             n.sims=as.integer(dims["n.sims"]), 
             sims.list=sims,
             link.power=link.power,
             call=call,
             formula=formula,
             model.frame = fr$mf,
             contrasts=contrasts,
             inits = inits)  
  return(ans)
                  
}        
