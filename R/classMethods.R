
################################################
# classes defined in the cplm package
################################################

# virtual classes used in other class definitions
setClassUnion("NullNum",c("NULL","numeric"))
setClassUnion("NullList",c("NULL","list"))  
setClassUnion("NullFunc",c("NULL","function"))  

# import from package coda
setOldClass(c("mcmc","mcmc.list","summary.mcmc"))

# class defining slots common to all derived classes 
setClass("cplm", 
  representation(
    call="call",
    formula="formula",
    contrasts="NullList",
    link.power="numeric",
    model.frame = "data.frame",
    inits="NullList")
)

# class of "cpglm", returned by a call to "cpglm" 
setClass("cpglm", 
  representation(
    coefficients="numeric",
    residuals="numeric",
    fitted.values="numeric",
    linear.predictors="numeric",
    y = "numeric",
    offset = "NullNum",
    prior.weights = "NullNum",
    weights="numeric",
    df.residual="integer",
    deviance="numeric",
    aic="numeric",
    control="list",
    p="numeric",
    phi="numeric",
    iter = "integer",
    converged = "logical",
    na.action = "NullFunc",
    vcov = "matrix"),
    contains = "cplm"
)

# class "cpglmm" returned from a call of cpglmm
setClass("cpglmm", 
 representation(
  p = "numeric", 
  phi = "numeric",
  bound.p = "numeric",
  vcov = "matrix"),
  contains = c("mer","cplm")
)

# class "summary.cpglmm" 
setClass("summary.cpglmm",                 
  representation(           
    methTitle = "character",
    logLik= "logLik",
    ngrps = "integer",
    sigma = "numeric", # scale, non-negative number
    coefs = "matrix",
    REmat = "matrix",
    AICtab= "data.frame"),
  contains = "cpglmm")

# class of "bcpglm", returned from a call to "bcpglm"
setClass("bcpglm", 
  representation(
    n.chains="integer",
    n.iter="integer", 
    n.burnin="integer",
    n.thin="integer", 
    n.sims="integer",  
    sims.list="mcmc.list",
    prop.var = "list"),
  contains="cplm")

# class of "bcpglmm", returned from a call to "bcpglmm"
setClass("bcpglmm", 
  representation(
    Zt = "dgCMatrix",
    flist = "list"),
  contains="bcpglm")
         
################################################
# methods defined for cplm 
################################################

# extraction of slots using $
setMethod("$",
    signature(x = "cplm"),
    function (x, name) 
        slot(x,name)
)

# names to get slot names
setMethod("names",
    signature(x = "cplm"),
    function (x) 
        return(slotNames(x))
)

# extraction of slots using "[["
setMethod("[[",
    signature(x = "cplm",i="numeric",j="missing"),
    function (x, i, j, ...) 
	    return(slot(x,names(x)[i]))
)

setMethod("[[",
    signature(x = "cplm",i="character",j="missing"),
    function (x, i, j, ...) 
      return(slot(x,i))
)

setMethod("[",
    signature(x = "cplm",i="numeric",j="missing",drop="missing"),
    function (x, i, j, ..., drop) {
      output <- lapply(i, function(y) slot(x,names(x)[y]))
        names(output) <- names(x)[i]
	    return(output)
    }
)

setMethod("[",
    signature(x = "cplm",i="character",j="missing",drop="missing"),
    function (x, i, j, ..., drop) {
      output <- lapply(1:length(i), function(y) slot(x,i[y]))
      names(output) <- i
      return(output)
    }
)

setMethod("terms",
    signature(x = "cplm"),
    function (x,...) 
      attr(x@model.frame,"terms")
)

setMethod("model.matrix",
    signature(object = "cplm"),
    function (object,...) 
      model.matrix(attr(object@model.frame,"terms"), 
            object@model.frame, object@contrasts)
)

setMethod("formula",
    signature(x = "cplm"),
    function (x,...) 
     x@formula
)         

setMethod("show",signature(object = "cplm"),
  function(object)
    summary(object)                                                    
)


setMethod("vcov", signature(object = "cplm"),
    function(object, ...){            
      object@vcov
})

################################################
# methods defined for cpglm 
################################################
         
setMethod("coef",
          signature(object = "cpglm"),
    function (object,...) 
	    return(object@coefficients)
)

setMethod("residuals",
    signature(object = "cpglm"),
    function (object,type = c("deviance", "pearson", "working", 
    "response", "partial"),...) {      
    type <- match.arg(type)
    y <- object@y
    r <- object@residuals
    mu <- object@fitted.values
    wts <- object@prior.weights
    family <- tweedie(var.power=object@p,link.power=object@link.power)
    switch(type, deviance = , pearson = , response = if (is.null(y)) {
        eta <- object@linear.predictors
        y <- mu + r * family$mu.eta(eta)
    })
    res <- switch(type, 
      deviance = if (object@df.residual > 0) {
        d.res <- sqrt(pmax((family$dev.resids)(y, mu, 
            wts), 0))
        ifelse(y > mu, d.res, -d.res)
        } else rep.int(0, length(mu)), 
      pearson = (y - mu) * sqrt(wts)/sqrt(family$variance(mu)), 
      working = r, 
      response = y - mu, 
      partial = r)
    na.action <- attr(object@model.frame,"na.action")
    if (!is.null(na.action)) 
        res <- naresid(na.action, res)
    #if (type == "partial") 
    #    res <- res + predict(object, type = "terms")
    res
    }
)

setMethod("resid",
    signature(object = "cpglm"),
    function (object, type = c("deviance", "pearson", "working", 
    "response", "partial"),...) 
	    return(residuals(object, type=type))
)

# generate fitted values on the original scale
setMethod("fitted",
    signature(object = "cpglm"),
    function (object,...) 
      return(object@fitted.values)
)
	

setMethod("AIC",
    signature(object = "cpglm",k="missing" ),
    function (object,...,k) 
      object@aic
)


setMethod("deviance",
    signature(object = "cpglm"),
    function (object,...) 
      object@deviance
)


setMethod("summary", signature(object="cpglm"),
	function(object,...){
    coef.beta <- coef(object)  
    vc <- vcov(object)
    s.err <- sqrt(diag(vc))    
    err.beta <- s.err
    test.value <- coef.beta/err.beta
    dn <- c("Estimate", "Std. Error")             
    pvalue <- 2 * pt(-abs(test.value), object@df.residual)
    
    coef.table <- cbind(coef.beta, err.beta, test.value, pvalue)  
    dn2 <- c("t value", "Pr(>|t|)")
    dimnames(coef.table) <- list(names(coef.beta), c(dn, dn2))
    keep <- match(c("call", "deviance", "aic", "contrasts", "df.residual",  
        "iter","na.action"), names(object), 0L)  
    ans <- c(object[keep], list(deviance.resid = residuals(object, 
        type = "deviance"), coefficients = coef.table, 
        dispersion = object@phi, vcov=vc, p=object@p))    
    .print.cpglm.summary(ans)    
    }
)

.print.cpglm.summary<-function(x,digits=max(3, getOption("digits") - 3),
                               signif.stars = getOption("show.signif.stars"), ...){
  
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    cat("Deviance Residuals: \n")
    if (x$df.residual > 5) {
        x$deviance.resid <- quantile(x$deviance.resid, na.rm = TRUE)
        names(x$deviance.resid) <- c("Min", "1Q", "Median", "3Q", 
            "Max")
    }
    xx <- zapsmall(x$deviance.resid, digits + 1)
    print.default(xx, digits = digits, na.print = "", print.gap = 2)
    printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars, 
            na.print = "NA",...)
        
    cat("\n(MLE estimate for the dispersion parameter is ",  
        format(x$dispersion,digits = max(5, digits + 1)), ";") 
    cat("\n MLE estimate for the index parameter is ",  
        format(x$p,digits = max(5, digits + 1)),")\n\n") 
    cat("Residual deviance:", format(x$deviance, digits = max(5, digits + 1)), 
        " on", format(x$df.residual), " degrees of freedom\n") 
    if (nzchar(mess <- naprint(x$na.action))) 
        cat("  (", mess, ")\n", sep = "")
    cat("AIC: ", format(x$aic, digits = max(4, digits + 1)), "\n\n")
    cat("Number of Fisher Scoring iterations: ", x$iter, "\n") 
    cat("\n")
    invisible(x)
}

# simple prediction method for cpglm
setMethod("predict", signature(object = "cpglm"),
  function (object, newdata, type = c("response", "link"), 
                  na.action = na.pass, ...) {
    tt <- attr(object@model.frame,"terms")
    if (missing(newdata) || is.null(newdata)) {
        X <- model.matrix(object)
        offset <- object$offset
    }
    else {
        Terms <- delete.response(tt)
        m <- model.frame(Terms, newdata, na.action = na.action)
        X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
        offset <- rep(0, nrow(X))
        if (!is.null(off.num <- attr(tt, "offset"))) 
            for (i in off.num) offset <- offset + eval(attr(tt, 
                "variables")[[i + 1]], newdata)
        if (!is.null(object$call$offset)) 
            offset <- offset + eval(object$call$offset, newdata)
    }
    beta <- object$coefficients
    predictor <- X%*% beta
    if (!is.null(offset)) 
        predictor <- predictor + offset
    mu <- tweedie(link.power=object@link.power)$linkinv(predictor)
    type <- match.arg(type)                                                            
    switch(type,link=predictor, response=mu)                                                            
})


# method for mcmcsamp, based on bcpglm
setMethod("mcmcsamp",
    signature(object = "cpglm"),
    function(object, inits = NULL, n.chains=3, n.iter=2000, 
        n.burnin=floor(n.iter/2),
        n.thin=max(1, floor(n.chains * (n.iter - n.burnin) / n.sims)),
        n.sims=1000, n.report=1000, 
        prior.beta.mean=NULL, prior.beta.var=NULL, 
        bound.phi=100, bound.p=c(1.01,1.99),  
        tune.iter=4000, n.tune=10, tune.weight=0.25,
        method="dtweedie", #prop.var = NULL,
             ...) {
  call <- match.call()
  X <- model.matrix(attr(object@model.frame,"terms"), 
            object@model.frame, object@contrasts)
  n.obs <- nrow(X)
  n.beta <- ncol(X)
  # default prior info
  if (is.null(prior.beta.mean))
    prior.beta.mean <- rep(0, n.beta)
  if (is.null(prior.beta.var))
    prior.beta.var <- rep(10000, n.beta)
 
 #  offset and prior wts
  wts <- object$prior.weights 
  off <- object$offset
  if (is.null(wts))     
      wts <- rep(1, n.obs)
  if (is.null(off)) 
        off <- rep(0, n.obs) 
  
  # dimensions used in simulation  
  n.keep <- floor((n.iter-n.burnin) / n.thin)
  n.sims <- n.chains * n.keep  
  dims <- list(n.obs= as.integer(n.obs),
           n.beta=as.integer(n.beta),
           n.pos= as.integer(sum(object$y>0)),     
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
 
  # proposal covariance matrix
  #if (is.null(prop.var)){
    C <- 2.38 * 2.38
    vc <- vcov(object)
    ebeta.var <- vc*C/n.beta
    ephi.var <- C * attr(vc, "phi_p")[1,1]
    ep.var <- C * attr(vc, "phi_p")[2,2]
  #} else {
  #  ebeta.var <- prop.var$beta.var
  #  ephi.var <- prop.var$phi.var
  #  ep.var <- prop.var$p.var
  #}
  # generate initial values if necessary
  if (is.null(inits)) {
      pstart <- object$p
      betastart <- as.numeric(object$coefficients)
      phistart <- object$phi
		  inits.start <- c(betastart, phistart, pstart)
		  inits <- vector("list",n.chains)
		  inits[[1]] <- c(betastart, phistart, pstart)
	    if (n.chains>1){
		    for (i in 2:n.chains)
			    inits[[i]] <- c(betastart + rnorm(n.beta,0,0.5),
						        runif(1,phistart/2,1.5*phistart),
						        runif(1,(bound.p[1]+pstart)/2,(bound.p[2]+pstart)/2))
	    }
  } else {
    inits <- lapply(inits, function(x) c(x$beta, x$phi, x$p ))
  }

  # run MCMC   
    # input for the C function 	     
    input <- list(X=as.double(X),
               y=as.double(object$y),
               ygt0= as.integer(which(object$y>0L)-1),
               offset=as.double(off),
               pWt=as.double(wts),
               mu = double(dims["n.obs"]),
               eta = double(dims["n.obs"]),
               inits = inits,
               beta=as.double(inits[[1]][1:n.beta]),
               phi=as.double(inits[[1]][n.beta+1]),
               p=as.double(inits[[1]][n.beta+2]),
               link.power=as.double(object$link.power),
               pbeta.mean=as.double(prior.beta.mean),
               pbeta.var=as.double(prior.beta.var),
               bound.phi=as.double(bound.phi),
               bound.p=as.double(bound.p),    
               ebeta.var=as.double(ebeta.var),
               ep.var= as.double(ep.var),
               ephi.var = as.double(ephi.var),               
               dims=dims,
               tune.weight=as.double(tune.weight),               
               lambda = as.double(2),
               k = as.integer(1),
               simT = as.integer(rep(1,dims["n.pos"])))
  
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
             link.power=object$link.power,
             call=call,
             formula=object$formula,
             model.frame = object$model.frame,
             contrasts=object$contrasts,
             inits = inits,
             prop.var = list(beta.var = matrix(input$ebeta.var,n.beta, n.beta), 
                             phi.var = input$ephi.var,
                             p.var = input$ep.var))  
  return(ans)
                  
})        


################################################
# methods defined for bcpglm
################################################

setMethod("summary", signature(object="bcpglm"),
  function(object, quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975),...){
    summary(object@sims.list,quantiles=quantiles)   
  }
)

setMethod("plot", signature(x="bcpglm",y="missing"),
  function(x,y,...)
    plot(x@sims.list)
)

        
################################################
# methods defined for cpglmm
################################################


setMethod("vcov", signature(object = "cpglmm"),
  function(object, ...){
    rr <- object$phi * chol2inv(object@RX, size = object@dims['p'])
    nms <- colnames(object@X)
    dimnames(rr) <- list(nms, nms)
      
    # compute vcov for phi and p numerically 
    cpglmm_dev <- function(x, ...){
      parm <- c(.Call(lme4:::mer_ST_getPars, object), 
                object$fixef, log(x[1]), x[2])
      .Call("cpglmm_update_dev",object, parm)  
    }
    x <- c(object$phi, object$p)
    hs <- hess(x, cpglmm_dev)
    dimnames(hs) <- list(c("phi","p"),c("phi","p"))      
    attr(rr,"phi_p") <- solve(hs)
    rr
})


setMethod("VarCorr", signature(x = "cpglmm"),
    function(x, ...){
    sc <- sqrt(x@phi)
	  ans <- lapply(cc <- .Call(lme4:::mer_ST_chol, x),
                        function(ch) {
                            val <- crossprod(sc * ch) # variance-covariance
                            stddev <- sqrt(diag(val))
                            correl <- t(val / stddev)/stddev
                            diag(correl) <- 1
                            attr(val, "stddev") <- stddev
                            attr(val, "correlation") <- correl
                            val
                        })
          fl <- x@flist
          names(ans) <- names(fl)[attr(fl, "assign")]
          attr(ans, "sc") <- sc
          ans
      })


setMethod("summary", signature(object = "cpglmm"),
    function(object, ...){
      fcoef <- fixef(object)
      vcov <- object@vcov
      dims <- object@dims
      coefs <- cbind("Estimate" = fcoef, "Std. Error" = sqrt(diag(vcov)) )
      llik <- logLik(object, 0)
      dev <- object@deviance
      mType <- "LMM"
      mName <- "Compound Poisson linear"
	    method <- paste("the", if(dims[["nAGQ"]] == 1) "Laplace" else
			  "adaptive Gaussian Hermite","approximation")
	  
      AICframe <- data.frame(AIC = AIC(llik), BIC = BIC(llik),
                                 logLik = as.vector(llik),
                                 deviance = dev[["ML"]],
                                 row.names = "")
      
      varcor <- VarCorr(object)
      REmat <- lme4:::formatVC(varcor)
      if (is.na(attr(varcor, "sc")))
          REmat <- REmat[-nrow(REmat), , drop = FALSE]

      if (nrow(coefs) > 0) {
        if (!dims[["useSc"]]) {
          coefs <- coefs[, 1:2, drop = FALSE]
          stat <- coefs[,1]/coefs[,2]
          pval <- 2*pnorm(abs(stat), lower = FALSE)
          coefs <- cbind(coefs, "z value" = stat, "Pr(>|z|)" = pval)
        } else {
          stat <- coefs[,1]/coefs[,2]
          ##pval <- 2*pt(abs(stat), coefs[,3], lower = FALSE)
          coefs <- cbind(coefs, "t value" = stat) #, "Pr(>|t|)" = pval)
        }
      } 
      new("summary.cpglmm",
              object,
              methTitle = paste(mName, "mixed model fit by", method),
              logLik = llik,
              ngrps = sapply(object@flist, function(x) length(levels(x))),
              sigma = sqrt(object@phi),
              coefs = coefs,
              REmat = REmat,
              AICtab = AICframe)
  }
)

## This is modeled a bit after  print.summary.lm :
print.cpglmm <- function(x, digits = max(3, getOption("digits") - 3),
                     correlation = FALSE, symbolic.cor = FALSE,
                     signif.stars = getOption("show.signif.stars"), ...){
    so <- summary(x)
    llik <- so@logLik
    dev <- so@deviance
    dims <- x@dims
    cat(so@methTitle, "\n")
    if (!is.null(x@call$formula))
        cat("Formula:", deparse(x@call$formula),"\n")
    if (!is.null(x@call$data))
        cat("   Data:", deparse(x@call$data), "\n")
    if (!is.null(x@call$subset))
        cat(" Subset:", deparse(x@call$subset),"\n")
    print(so@AICtab, digits = digits)

    cat("Random effects:\n")
    print(so@REmat, quote = FALSE, digits = digits, ...)

    ngrps <- so@ngrps
    cat(sprintf("Number of obs: %d, groups: ", dims[["n"]]))
    cat(paste(paste(names(ngrps), ngrps, sep = ", "), collapse = "; "))
    cat("\n")
     
    if (nrow(so@coefs) > 0) {
	    cat("\nFixed effects:\n")
	    printCoefmat(so@coefs, zap.ind = 3, #, tst.ind = 4
		     digits = digits, signif.stars = signif.stars)
  
      cat("\nEstimated scale parameter:", round(so@phi, digits=digits))
      cat("\n")
      cat("Estimated index parameter:", round(so@p, digits=digits))
      cat("\n")
     
	    if(correlation) {
	      corF <- so@vcov@factors$correlation
	    if (!is.null(corF)) {
		    p <- ncol(corF)
		    if (p > 1) {
		      rn <- rownames(so@coefs)
		      rns <- abbreviate(rn, minlength=11)
		      cat("\nCorrelation of Fixed Effects:\n")
		      if (is.logical(symbolic.cor) && symbolic.cor) {
			      corf <- as(corF, "matrix")
			      dimnames(corf) <- list(rns,
					       abbreviate(rn, minlength=1, strict=TRUE))
			      print(symnum(corf))
		      } else {
			      corf <- matrix(format(round(corF@x, 3), nsmall = 3),
				       ncol = p,dimnames = list(rns, abbreviate(rn, minlength=6)))
			      corf[!lower.tri(corf)] <- ""
			      print(corf[-1, -p, drop=FALSE], quote = FALSE)
		      }
		    }
	    }
	  }
  }
  invisible(x)
}

setMethod("print", "cpglmm", print.cpglmm)
setMethod("show", "cpglmm", 
  function(object) 
    print.cpglmm(object)
)


# predict method for cpglmm
getZt <- function(formula, oldmf, newmf){
  bars <- lme4:::expandSlash(lme4:::findbars(formula[[3]]))
  names(bars) <- unlist(lapply(bars, function(x) deparse(x[[3]])))
  fl <- lapply(bars, function(x) {
        oldlvl <- eval(substitute(levels(as.factor(fac)[, drop = TRUE]), 
            list(fac = x[[3]])), oldmf)
        ff <- eval(substitute(factor(fac,levels=oldlvl)[, drop = TRUE], 
            list(fac = x[[3]])), newmf)
        # fill columns of 0's if some levels are missing
        im <- as(ff, "sparseMatrix")
        im2 <- Matrix(0,nrow=length(oldlvl),ncol=length(ff),sparse=T)
        # this is awkward as the Matrix package seems to fail
        for (i in 1:nrow(im)){
          ind <- match(rownames(im)[i],oldlvl)
          im2[as.numeric(ind),] <- im[as.numeric(i),]            
        }        
        if (!isTRUE(validObject(im, test = TRUE))) 
            stop("invalid conditioning factor in random effect: ", 
                format(x[[3]]))
        mm <- model.matrix(eval(substitute(~expr, list(expr = x[[2]]))), newmf)
        mm <- mm[!is.na(ff),,drop=F]
        Zt <- do.call(rBind, lapply(seq_len(ncol(mm)), 
            function(j) {
                im2@x <- mm[, j]
                im2
            }))
        ans <- list(f = oldlvl, Zt=Zt)       
        ans
    })
  nlev <- sapply(fl, function(el) length(levels(el$f)))
  if (any(diff(nlev)) > 0) 
        fl <- fl[rev(order(nlev))]        
  Zt <- do.call(rBind,lapply(fl,"[[","Zt"))
  Zt
}         

setMethod("predict", signature(object = "cpglmm"),
    function(object,  newdata, type = c("link","response"), 
          na.action = na.pass, ...) {
    tt <- attr(object@model.frame,"terms")
    if (missing(newdata) || is.null(newdata)) {
        mm <- X <- model.matrix(object)
        Zt <- object@Zt
        offset <- object$offset
    }
    else {
        Terms <- delete.response(tt)
        # design matrix for fixed effects
        X <- model.matrix(Terms, newdata, contrasts.arg = object@contrasts)
        # design matrix for random effects
        formula <- object@formula
        oldmf <- object@model.frame
        Zt <- getZt(formula, oldmf, newdata)        
        # get offset
        offset <- rep(0, nrow(X))
        if (!is.null(off.num <- attr(tt, "offset"))) 
            for (i in off.num) offset <- offset + eval(attr(tt, 
                "variables")[[i + 1]], newdata)
        if (!is.null(object$call$offset)) 
            offset <- offset + eval(object$call$offset, newdata)                
    }
    beta <- object@fixef
    u <- object@ranef
    predictor <- as.numeric(X %*% beta + t(Zt)%*% u)
    if (!is.null(offset)) 
        predictor <- predictor + offset
    mu <- tweedie(link.power=object@link.power)$linkinv(predictor)
    type <- match.arg(type)
    switch(type,link=predictor, response=mu)   
})


# methods for mcmcsamp, based on bcpglmm
setMethod("mcmcsamp",
    signature(object = "cpglmm"),
    function(object, inits = NULL, n.chains=3, n.iter=2000, 
        n.burnin=floor(n.iter/2),
        n.thin=max(1, floor(n.chains * (n.iter - n.burnin) / n.sims)),
        n.sims=1000, n.report=1000, 
        prior.beta.mean=NULL, prior.beta.var=NULL, 
        bound.phi=100, bound.p=c(1.01,1.99),  prior.Sigma = NULL,
        tune.iter=4000, n.tune=10, tune.weight=0.25,
        method="dtweedie", #prop.var = NULL, 
             ...) {
  call <- match.call()
  dd <- object$dims
  n.obs <- as.integer(dd['n'])
  n.beta <- as.integer(dd['p'])
  # default prior info
  if (is.null(prior.beta.mean))
    prior.beta.mean <- rep(0, n.beta)
  if (is.null(prior.beta.var))
    prior.beta.var <- rep(10000, n.beta)
  if (is.null(prior.Sigma))
    prior.Sigma <- prior.Sigma.default(object$ST)
 
 #  offset and prior wts
  wts <- object$pWt 
  off <- object$offset
  
  # dimensions used in simulation
  nc <- unlist(lapply(object$ST, ncol))
  nlev <- diff(object$Gp)/nc  
  n.keep <- floor((n.iter-n.burnin) / n.thin)
  n.sims <- n.chains * n.keep 
  dims <- list(n.obs= n.obs,
           n.beta=as.integer(unname(n.beta)),           
           n.pos= as.integer(sum(object$y>0)),
           n.term = as.integer(dd['nt']),
           n.u = as.integer(dd['q']),
           n.all = as.integer(dd['p'] + dd['q'] + sum(nc^2)+2),
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
  
  # proposal covariance matrix              
  #if (is.null(prop.var)){
    C <- 2.38 * 2.38
    vc <- vcov(object)
    ebeta.var <- vc * C /n.beta
    tmp <- unlist(lapply(ranef(object, postVar=T), 
                       function(x) as.numeric(attr(x, 'postVar'))))
    eu.var <- diag(as.numeric(tmp), nrow=dd['q']) * C/dd['q']
    ephi.var <- C * attr(vc, "phi_p")[1,1]
    ep.var <- C * attr(vc, "phi_p")[2,2]
  #} else {
  #  ebeta.var <- prop.var$beta.var
  #  eu.var <- prop.var$u.var
  #  ephi.var <- prop.var$phi.var
  #  ep.var <- prop.var$p.var
  #}

                   
  # generate initial values 
  if (is.null(inits)) {
    pstart <- object$p
    betastart <- as.numeric(fixef(object))
    ustart <- lapply(ranef(object), function(t) as.numeric(unlist(t)))
    ustart <- as.numeric(unlist(ustart))
    phistart <- object$phi
    Sigmastart <- lapply(object$ST,function(x) x%*%t(x))
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
  input <- list(X=object$X,
               y=as.double(object$y),
               Zt = object$Zt, 
               ygt0= as.integer(which(object$y>0L)-1),
               offset=as.double(off),
               pWt=as.double(wts),
               mu = double(dims["n.obs"]),
               eta = double(dims["n.obs"]),
               inits = inits,
               beta=as.double(betastart),
               u= as.double(ustart),
               phi=as.double(phistart),
               p=as.double(pstart),
               link.power=as.double(object$link.power),
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
               Gp = unname(object$Gp),
               Sigma = Sigmastart,
               ncol = as.integer(nc), 
               nlev = as.integer(nlev),
               lambda = as.double(2),
               k = as.integer(1),
               simT = as.integer(rep(1,dims["n.pos"])),
               mh.var = c(diag(ebeta.var), diag(eu.var), ephi.var, ep.var))

  if (method=="dtweedie")
    sims.list<- .Call("bcpglmm_gibbs_tw",input) 
  if (method=="latent")
    sims.list<- .Call("bcpglmm_gibbs_lat",input)
  xnames <- names(fixef(object))
  unames <- paste("u", 1:dd['q'],sep="")
  snames <- sapply(1:length(object$ST), function(x){
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
             link.power=object$link.power,
             call=call,
             formula= object$formula,
             model.frame = object$model.frame,
             contrasts=object$contrasts,
             inits = inits,
             prop.var = list(beta.var = matrix(input$ebeta.var,n.beta,n.beta), 
                             u.var = matrix(input$eu.var,dd['q'],dd['q']),
                             phi.var = input$ephi.var,
                             p.var = input$ep.var))  
  return(ans)
                  
})