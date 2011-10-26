
################################################
# classes defined in the cplm package
################################################

# virtual classes used in other class definitions
setClassUnion("NullNum",c("NULL","numeric"))
setClassUnion("NullList",c("NULL","list"))  
setClassUnion("NullFunc",c("NULL","function"))

# import from package coda
setOldClass(c("mcmc","mcmc.list","summary.mcmc"))

# class of "cpglm", returned by a call to "cpglm" 
setClass("cpglm", 
 representation(
  coefficients="numeric",
  residuals="numeric",
  fitted.values="numeric",
  linear.predictors="numeric",
  weights="numeric",
  df.residual="integer",
  deviance="numeric",
  aic="numeric",
  offset="NullNum",
  prior.weights="NullNum",
  call="call",
  formula="formula",
  data="data.frame",
  control="list",
  contrasts="NullList",
  theta="numeric",
  theta.all="matrix",
  p="numeric",
  phi="numeric",
  vcov="matrix",
  iter="integer",
  converged="logical",
  method="character",
  y="numeric",
  link.power="numeric",
  na.action="NullFunc",
  model.frame="data.frame")
)

# class "cpglmm" returned from a call of cpglmm
setClass("cpglmm", 
 representation(
  p="numeric", 
  phi="numeric",
  link.power="numeric", 
  bound.p="numeric"),
  contains = "mer"
)

# class "summary.cpglmm" 
setClass("summary.cpglmm",                 
  representation(           
    methTitle = "character",
    logLik= "logLik",
    ngrps = "integer",
    sigma = "numeric", # scale, non-negative number
    coefs = "matrix",
    vcov = "dpoMatrix",
    REmat = "matrix",
    AICtab= "data.frame"),
    contains = "cpglmm")


# class of "bcplm", base class used in all Bayesian models 
setClass("bcplm", 
 representation( 
  n.chains="integer",
  n.iter="integer", 
  n.burnin="integer",
  n.thin="integer", 
  n.sims="integer",  
  sims.list="mcmc.list",
  link.power="numeric",
  model.frame = "data.frame",
  call="call",
  formula="formula",
  data="data.frame",
  contrasts="NullList",
  inits="list")
)             

# class of "bcpglm", returned from a call to "bcpglm"
setClass("bcpglm", contains="bcplm")             

################################################
# methods defined for cpglm 
################################################

# extraction of slots using $
setMethod("$",
    signature(x = "cpglm"),
    function (x, name) 
        slot(x,name)
)

# names to get slot names
setMethod("names",
    signature(x = "cpglm"),
    function (x) 
        return(slotNames(x))
)

# extraction of slots using "[["
setMethod("[[",
    signature(x = "cpglm",i="numeric",j="missing"),
    function (x, i, j, ...) 
	    return(slot(x,names(x)[i]))
)

setMethod("[[",
    signature(x = "cpglm",i="character",j="missing"),
    function (x, i, j, ...) 
      return(slot(x,i))
)

setMethod("[",
    signature(x = "cpglm",i="numeric",j="missing",drop="missing"),
    function (x, i, j, ..., drop) {
      output <- lapply(i, function(y) slot(x,names(x)[y]))
        names(output) <- names(x)[i]
	    return(output)
    }
)

setMethod("[",
    signature(x = "cpglm",i="character",j="missing",drop="missing"),
    function (x, i, j, ..., drop) {
      output <- lapply(1:length(i), function(y) slot(x,i[y]))
      names(output) <- i
      return(output)
    }
)

setMethod("coef",
          signature(object = "cpglm"),
    function (object,...) 
	    return(object@coefficients)
)

setMethod("vcov",
	signature(object = "cpglm"),
    function (object,...) 
    	return(object@vcov)
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
    if (!is.null(object@na.action)) 
        res <- naresid(object@na.action, res)
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
		
setMethod("fitted.values",
    signature(object = "cpglm"),
    function (object,...) 
      fitted(object)
)

  	
setMethod("df.residual",
    signature(object = "cpglm"),
    function (object,...) 
      object@df.residual
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

setMethod("terms",
    signature(x = "cpglm"),
    function (x,...) 
      attr(x@model.frame,"terms")
)

setMethod("model.matrix",
    signature(object = "cpglm"),
    function (object,...) 
      model.matrix(terms(object), 
            object@model.frame, object@contrasts)
)

setMethod("formula",
    signature(x = "cpglm"),
    function (x,...) 
     x@formula
)

setMethod("summary", signature(object="cpglm"),
	function(object,...){
    coef.beta <- coef(object)    
    s.err <- sqrt(diag(object@vcov))    
    err.beta <- switch(object@method, 
                        MCEM=s.err[1:(length(s.err)-2)],
                        profile=s.err)
    test.value <- coef.beta/err.beta
    dn <- c("Estimate", "Std. Error")             
    pvalue <- switch(object@method, 
                        MCEM=2 * pnorm(-abs(test.value)),
                        profile=2 * pt(-abs(test.value), object@df.residual))
    
    coef.table <- cbind(coef.beta, err.beta, test.value, pvalue)  
    dn2 <- switch(object@method, 
                        MCEM=c("z value", "Pr(>|z|)"),
                        profile=c("t value", "Pr(>|t|)"))
    dimnames(coef.table) <- list(names(coef.beta), c(dn, dn2))
    keep <- match(c("call", "deviance", "aic", "contrasts", "df.residual","method",  
        "iter", "na.action"), names(object), 0L)  
    ans <- c(object[keep], list(deviance.resid = residuals(object, 
        type = "deviance"), coefficients = coef.table, 
        dispersion = object@phi, vcov=object@vcov, p=object@p))    
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
    if (x$method=="MCEM")
      cat("Number of Monte Carlo EM iterations: ", x$iter, "\n") 
    if (x$method=="profile")
      cat("Number of Fisher Scoring iterations: ", x$iter, "\n") 
    cat("\n")
    invisible(x)
}
    
setMethod("show",signature(object = "cpglm"),
  function(object)
    summary(object)                                                    
)     



################################################
# methods defined for bcplm
################################################

# extraction of slots using $
setMethod("$",
    signature(x = "bcplm"),
    function (x, name) 
        slot(x,name)
)

# names to get slot names
setMethod("names",
    signature(x = "bcplm"),
    function (x) 
        return(slotNames(x))
)

# extraction of slots using "[["
setMethod("[[",
    signature(x = "bcplm",i="numeric",j="missing"),
    function (x, i, j, ...) 
      return(slot(x,names(x)[i]))
)

setMethod("[[",
    signature(x = "bcplm",i="character",j="missing"),
    function (x, i, j, ...) 
      return(slot(x,i))
)

setMethod("[",
    signature(x = "bcplm",i="numeric",j="missing",drop="missing"),
    function (x, i, j, ..., drop) {
      output <- lapply(i, function(y) slot(x,names(x)[y]))
        names(output) <- names(x)[i]
	    return(output)
    }
)

setMethod("[",
    signature(x = "bcplm",i="character",j="missing",drop="missing"),
    function (x, i, j, ..., drop) {
      output <- lapply(1:length(i), function(y) slot(x,i[y]))
      names(output) <- i
      return(output)
    }
)

setMethod("terms",
    signature(x = "bcplm"),
    function (x,...) 
      attr(x@model.frame,"terms")
)

setMethod("model.matrix",
    signature(object = "bcplm"),
    function (object,...) 
      model.matrix(terms(object), 
            object@model.frame, object@contrasts)
)

setMethod("formula",
    signature(x = "bcplm"),
    function (x,...) 
      x@formula
)


setMethod("summary", signature(object="bcplm"),
  function(object, quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975),...){
    summary(object@sims.list,quantiles=quantiles)   
  }
)


setMethod("plot", signature(x="bcplm",y="missing"),
  function(x,y,...)
    plot(x@sims.list)
)


setMethod("show",signature(object = "bcplm"),
  function(object)
    summary(object)                                                    
)


################################################
# methods defined for cpglmm
################################################


setMethod("vcov", signature(object = "cpglmm"),
    function(object, ...){
      rr <- as(object@phi *
              chol2inv(object@RX, size = object@dims['p']), "dpoMatrix")
      nms <- colnames(object@X)
      dimnames(rr) <- list(nms, nms)
      rr@factors$correlation <- as(rr, "corMatrix")
      rr
})

setMethod("VarCorr", signature(x = "cpglmm"),
    function(x, ...){
### Create the VarCorr object of variances and covariances
      sc <- sqrt(x@phi)
	    ans <- lapply(cc <- .Call(lme4:::mer_ST_chol, x),
                        function(ch) {
                            val <-  ch %*%t(ch) # variance-covariance
                            stddev <- sqrt(diag(val))
                            correl <- t(val / stddev)/stddev
                            diag(correl) <- 1
                            attr(val, "stddev") <- stddev
                            attr(val, "correlation") <- correl
                            val
                        })
      fl <- x@flist
      names(ans) <- names(fl)[attr(fl, "assign")]
      attr(ans, "sc") <- if (x@dims[["useSc"]]) sc else NA
      ans
    }
)

setMethod("summary", signature(object = "cpglmm"),
    function(object, ...){
      REML <- object@dims[["REML"]]
      fcoef <- fixef(object)
      vcov <- vcov(object)
      corF <- vcov@factors$correlation
      dims <- object@dims
      coefs <- cbind("Estimate" = fcoef, "Std. Error" = corF@sd) #, DF = DF)
      llik <- logLik(object, REML)
      dev <- object@deviance
      mType <- if((non <- as.logical(length(object@V)))) "NMM" else "LMM"
      mName <- "Compound Poisson linear"
	    method <- paste("the", if(dims[["nAGQ"]] == 1) "Laplace" else
			  "adaptive Gaussian Hermite","approximation")
	  
      AICframe <- data.frame(AIC = AIC(llik), BIC = BIC(llik),
                                 logLik = as.vector(llik),
                                 deviance = dev[["ML"]],
                                 REMLdev = dev[["REML"]],
                                 row.names = "")
      if (is.na(AICframe$REMLdev)) AICframe$REMLdev <- NULL
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
          } ## else : append columns to 0-row matrix ...
          new("summary.cpglmm",
              object,
              methTitle = paste(mName, "mixed model fit by", method),
              logLik = llik,
              ngrps = sapply(object@flist, function(x) length(levels(x))),
              sigma = sqrt(object@phi),
              coefs = coefs,
              vcov = vcov,
              REmat = REmat,
              AICtab= AICframe)
  }
)

## This is modeled a bit after  print.summary.lm :
print.cpglmm <- function(x, digits = max(3, getOption("digits") - 3),
                     correlation = FALSE, symbolic.cor = FALSE,
                     signif.stars = getOption("show.signif.stars"), ...){
    so <- summary(x)
    REML <- so@dims[["REML"]]
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
		      rns <- abbreviate(rn, minlen=11)
		      cat("\nCorrelation of Fixed Effects:\n")
		      if (is.logical(symbolic.cor) && symbolic.cor) {
			      corf <- as(corF, "matrix")
			      dimnames(corf) <- list(rns,
					       abbreviate(rn, minlength=1, strict=TRUE))
			      print(symnum(corf))
		      } else {
			      corf <- matrix(format(round(corF@x, 3), nsmall = 3),
				       nc = p,dimnames = list(rns, abbreviate(rn, minlen=6)))
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
                  