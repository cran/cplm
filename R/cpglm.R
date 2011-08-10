#######################################################
## MCEM method for fitting compound Poisson GLM
#######################################################

cpglm <- function(formula, link = "log", data, weights, subset, 
    na.action,betastart=NULL, phistart=NULL, pstart=NULL, offset, 
    contrasts = NULL,   control=list(),
    ... ) {
  
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
  if (!is.null(weights)) 
        stop("'weights' is not implemented yet")
  if (!is.null(offset)) {
    if (length(offset) != NROW(Y)) 
      stop(gettextf("number of 'offset' is %d should equal %d (number of observations)", 
                length(offset), NROW(Y)), domain = NA)
    }
  if (!is.null(betastart)){
    if (length(betastart) != ncol(X))
      stop(gettextf("number of 'betastart' is %d should equal %d (number of mean parameters)", 
                length(betastart), ncol(X)), domain = NA)
    }
  if (!is.null(phistart) && phistart<=0)
    stop("value of 'phistart' should be greater than 0")
  if (!is.null(pstart) && (pstart<=1 || pstart>=2))
    stop("value of 'pstart' should be between 1 and 2")   
  link.power <- make.link.power(link)

  cpfit <- .cpglm.fit(X,Y,weights=weights,offset=offset,
                     link.power=link.power,contrasts=contrasts,
                     betastart=betastart,phistart=phistart,pstart=pstart,
                      control=control)
  ans <- new("cpglm", 
             coefficients=cpfit$coefficients, 
             residuals=cpfit$residuals,
             fitted.values=cpfit$fitted.values,
             weights=cpfit$weights,
             df.residual=cpfit$df.residual,
             df.null=cpfit$df.null,
             y=cpfit$y,
             call=call,
             formula=formula,
             data=data,
             offset=cpfit$offset,
             control=cpfit$control,
             contrasts=contrasts,
             theta=cpfit$theta,
             vcov=cpfit$vcov,
             iter=cpfit$iter)  
  return(ans)
}


# internal function to run the MCEM 
.cpglm.fit <- function(X,Y,weights=NULL,offset=NULL,
                      link.power=0,contrasts=NULL,
                      betastart,phistart,pstart,
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
    nobs <- NROW(Y)
    #FIXME: fix weights now. does not handle prior weights
    weights <- rep.int(1, nobs)
    if (is.null(offset)) 
        offset <- rep.int(0, nobs)          
    # generating starting values if necessary
    if (is.null(pstart)) 
      pstart <- 1.5
    if (is.null(betastart) || is.null(phistart)) {
      fit.start <- glm(Y~-1+X,weights=weights,offset=offset,
                  family=tweedie(var.power=pstart,link.power=link.power))
      if (is.null(betastart))
        betastart <- as.numeric(fit.start$coefficients)
      if (is.null(phistart))
        phistart <- sum(residuals(fit.start,"pearson")^2)/df.residual(fit.start)
    }
    
    # index of positive y's     
    ygt0 <- as.numeric(which(Y>0L) )     
    npobs <- length(ygt0)    
    nbeta <- length(betastart)
    sz <- control$init.size 
    sr <- control$sample.iter
    C <- qchisq(1-control$alpha, nbeta+2)
    
    # input list to be used in sampling or likelihood calculation
    input <- list(X=X,Y=Y,ygt0=ygt0,beta=betastart,
                  phi=phistart,p=pstart,k=ygt0[1L],
                  link.power=link.power,offset=offset)         
         
    # start iterations:
    diff <- 10
    iter <- 0
    theta <- as.numeric(unlist(input[c("beta","phi","p")]))                  
    # weight for each sample
    w <- rep.int(1,sz)
      	
    while (diff>control$epsilon2){
      iter <- iter + 1       
    ######################################################### 
    ###########            The E Step             ###########
    #########################################################        
    # sample the t's using rejection sampling if needed
    if (iter <=sr) {
		sim.t <- sample.t.rej(sz, input)
    } else                    
    # use importance sampling if iter > sample.iter
    {   
      w <- import.weight(theta,input,sim.t,sr) 
    }  
              
    ######################################################### 
    ###########            The M Step             ###########
    #########################################################            
    #1. compute beta given p
    fit <- glm.fit(x = X, y = Y, weights = weights,offset = offset, 
		family = tweedie(var.power=input$p,link.power=input$link.power))
    input$beta <- as.numeric(fit$coefficients)
    mu <- as.numeric(fit$fitted.values)
    
    #2. update phi and p
    nb <- optim.phi.p(mu,sim.t,w,bd=control$bound.p,input)
    input$phi <- nb[1] 
    input$p <- nb[2] 
 	
    # store the history of theta
    theta <- rbind(theta,as.numeric(unlist(input[c("beta","phi","p")])))                    

    # print out iteration info if necessary  
    if (control$verbose) {
      cat("Iteration:",iter,"\n")
      parm <- round(theta[iter+1,],4)
      names(parm) <- c(xnames,"phi","p")
      cat(parm,"\n")
    }
    if (iter >= control$maxit) {
      warning("maximum iteration limit reached!\n")
      break
    }
    # update diff
    diff <- max(abs(theta[iter+1,]-theta[iter,])/
    			(abs(theta[iter,])+control$epsilon1))
    
    # determine if sample size needs to be increased
    if (!control$fixed.size & iter <=sr){		
    	T <- tstat.intv.cpglm(theta,input, sim.t,w)  
    	if (iter <=sr &  T>C) 	
	    	sz <- min(sz + round(sz/control$k), control$max.size)
	cat("sample size:", sz, "\n")     
	}
  }

    # fit the glm using the converged parameters  
    fit <- glm(Y~-1+ X, weights = weights,offset = offset,
          family=tweedie(var.power=input$p,link.power=input$link.power))

    if (control$verbose) {
      cat("Algorithm has converged\n")
      cat("Computing hessian matrix\n")      
    }  

    # compute the vcov matrix  
    vcov <- var.theta.cpglm(theta,input, sim.t,w)

    # get output  
    names(fit$coefficients)<- xnames    
    ans <- c(fit[c("coefficients","residuals","fitted.values",
                   "df.residual", "df.null","y")],
                   list(theta=theta,vcov=vcov,iter=iter,
                        control=control,offset=offset,
                        weights=weights))
    return(ans)

}
  
# function to simulate T's from the posterior distribution
sample.t.rej <- function(n,input){
  out <- .Call("sampleTRej",n=as.integer(n),
		Y=as.numeric(input$Y[input$ygt0]),
		phi=input$phi,
		p=input$p )
  return(out)				
}

# function to compute the weight in importance sampling
import.weight <- function(theta,input,sim.t,sr){
  ix <- c(ncol(theta)-1, ncol(theta))      
  w <- .Call("importWeight",nc=as.integer(ncol(sim.t)),
      			Y=input$Y[input$ygt0], 
      			sT=sim.t,parm=list("old"=theta[sr,ix],
                                   "new"=theta[nrow(theta),ix]))  
}  

# function to maximize over phi and p
optim.phi.p <- function(mu,sim.t,w,bd=c(1.01,1.99),input){
  out <- .Call("optimNbGlm",Y=input$Y[input$ygt0], 
	mu=mu, ygt0=as.integer(input$ygt0), p= input$p, phi=input$phi,  
	simT=sim.t, w=w, bd=bd)
  return(out)				
}

# function to take inverse of a matrix using svd 
svd.inv <- function(x){
	sx <- svd(x)
	return(sx$v%*% diag(1/sx$d)%*%t(sx$u))	
}

# function to compute approximate normal confidence interval
# used for determining whether to increase sample size
tstat.intv.cpglm <- function(theta,input,sim.t,w){
	nt <- nrow(theta)
	out <- .Call("hessGlmEst",X=input$X,Y=input$Y[input$ygt0],
            ygt0= as.integer(input$ygt0),theta=theta[nt,],
		simT=sim.t,w=w,linkpower=as.double(input$link.power), 
		n=as.integer(nrow(input$X)),offset=as.numeric(input$offset))
	# compute approximate covariance matrix	
	hI.inv <- svd.inv(out$HQ)
	var.app <- hI.inv %*% out$VQ %*% t(hI.inv)  
	mu.app <- theta[nt-1,]-theta[nt,]
	# compute test stat
	return(t(mu.app)%*%svd.inv(var.app)%*%mu.app)
}

               
# function to compute covariance matrix
var.theta.cpglm <- function(theta,input,sim.t,w){
	out <- .Call("hessGlmEst",X=input$X,Y=input$Y[input$ygt0],
            ygt0= as.integer(input$ygt0),theta=theta[nrow(theta),],
		simT=sim.t,w=w,linkpower=as.double(input$link.power), 
		n=as.integer(nrow(input$X)),offset=as.numeric(input$offset))
	# compute information matrix 
	IM <- -out[[1]]-out[[2]]	
	# invert information matrix
	svdI <- svd(IM)              
	var.theta <- svdI$v%*% diag(1/svdI$d)%*%t(svdI$u)  
	return(var.theta)
}

# function to compute the link.power needed in tweedie
make.link.power <- function(link) {   
  #link.temp<- substitute(link)
  #if (!is.character(link.temp))
  #      link.temp <- deparse(link.temp)
  link.temp<- link
  okLinks <- c("log", "identity", "sqrt","inverse")
  if (link.temp %in% okLinks) 
      switch(link.temp,log=0, identity=1, sqrt=0.5, inverse=2) else
      stop("invalid link function!")
}

# control options intializer
cpglm.control <- function(init.size=100L,
                       sample.iter=50L,
                       max.size=10000L,
                       maxit=200,
                       epsilon1=1e-03,
                       epsilon2=1e-04,
                       alpha =0.25,
                       k=5,                       
                       bound.p=c(1.01,1.99),
                       fixed.size=TRUE,   
                       verbose=TRUE){
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
  if (!is.numeric(maxit) || maxit <= 0) 
        stop("value of 'maxit' must be > 0")
  if (min(bound.p)<1 || max(bound.p)>2)
        stop("value of 'bound.p' must be between 1 and 2")
  if (!is.numeric(fixed.size) && !is.logical(fixed.size))
        stop("'fixed.size' must be logical or numeric")
  if (!is.numeric(verbose) && !is.logical(verbose))
        stop("'verbose' must be logical or numeric")
  
  bound.p <- sort(bound.p)
  fixed.size <- as.logical(fixed.size)
  verbose <- as.logical(verbose)
  
    list(init.size=init.size,
         sample.iter=sample.iter,
         maxit=maxit,
         epsilon1 = epsilon1,
         epsilon2=epsilon2,
         alpha=alpha,
         k=k,
         fixed.size=fixed.size,
         max.size=max.size,
         bound.p=bound.p,
         verbose=verbose)  
}






