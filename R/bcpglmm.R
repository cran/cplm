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
                   tune.iter=4000, n.tune=10, tune.weight=0.25,
                   basisGenerators = c("tp","tpU","bsp","sp2d"),
                   method="dtweedie",...) {
      
  call <- amer:::expand.call()  
  if (missing(data)) 
    data <- environment(formula)   
  link.power <- make.link.power(link)
  # identify smooth terms 
  tf <- terms.formula(formula, specials = eval(call$basisGenerators, 
        parent.frame(2)))
  n.f <- length(unlist(attr(tf, "specials")))
  # create model frame and get factor list
  if (n.f) {
    call2 <- as.list(call)[-1]
    m <- match(c("formula", "data", "weights", "offset", 
                        "contrasts", "basisGenerators"), names(call2),0L)
    call2 <- call2[m]
    #setup <- do.call(amer:::amerSetup, as.list(call2))     
    setup <- do.call(frFL, as.list(call2))
    fr <- setup$m$fr 
    FL <- setup$m$FL
  } else {  
    fr <- lme4:::lmerFrames(call, formula, contrasts)
    FL <- lme4:::lmerFactorList(formula, fr, 0L, 0L)
  }
  
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
 # default offset and prior wts
 if (is.null(fr$wts) || length(fr$wts)==0)  
    wts <- as.double(rep(1,n.obs)) else 
    wts <- fr$wts 
  if (is.null(fr$off) || length(fr$off)==0)
    off <- as.double(rep(0,n.obs)) else 
    off <- fr$off
  
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
    pstart <- 1.5
    fit.start <- glm.fit(fr$X, fr$Y, weights = wts, offset = off, 
        family = tweedie(var.power=pstart,link.power=link.power), 
        intercept = attr(attr(fr$mf, "terms"), "intercept") > 0)
    betastart <- as.numeric(fit.start$coefficients)
    ustart <- rnorm(dm$dd[["q"]])
    phistart <- sum((fr$Y-fit.start$fitted.values)^2/fit.start$fitted.values^pstart)/
                    fit.start$df.residual
		Sigmastart <- lapply(dm$ST,function(x) x%*%t(x))
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
    #ebeta.var <- as.matrix(vcov(fit.start)*C/n.beta)
    eu.var <- lapply(1:dm$dd['nt'], function(x) 
              kronecker(Sigmastart[[x]],diag(1,nrow=nlev[x])))
    eu.var <- as.matrix(do.call(bdiag, eu.var))*C/dm$dd['q']
  } else{
    betastart <- inits[[1]]$beta
    phistart <- inits[[1]]$phi
    pstart <- inits[[1]]$p
    ustart <- inits[[1]]$u
    Sigmastart <- inits[[1]]$Sigma
    inits <- lapply(inits, function(x) c(x$beta, x$phi, x$p, x$u,
                                      unlist(lapply(x$Sigma, as.numeric))))
  }

  # run MCMC   
    # input for the C function 	     
    input <- list(X=fr$X,
               y=as.double(fr$Y),
               Zt = dm$Zt, 
               ygt0= as.integer(which(fr$Y>0L)-1),
               offset=as.double(off),
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
               nlev = as.integer(nlev),
               lambda = as.double(2),
               k = as.integer(1),
               simT = as.integer(rep(1,dims["n.pos"])))

  if (method=="dtweedie")
    sims.list<- .Call("bcpglmm_gibbs_tw",input) 
  if (method=="latent")
    sims.list<- .Call("bcpglmm_gibbs_lat",input)
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
             inits = inits,
             prop.var = list(beta.var = matrix(input$ebeta.var,n.beta,n.beta), 
                             u.var = matrix(input$eu.var,dm$dd['q'],dm$dd['q']),
                             phi.var = input$ephi.var,
                             p.var = input$ep.var))  
  return(ans)
                  
}        
