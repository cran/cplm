
################################################
# classes defined in the cplm package
################################################

setClassUnion("NullNum",c("NULL","numeric"))
setClassUnion("NullList",c("NULL","list"))	

# class of "cpglm" 
setClass("cpglm", 
 representation(
  coefficients="numeric",
  residuals="numeric",
  fitted.values="numeric",
  weights="NullNum",
  df.residual="integer",
  df.null="integer",
  y="numeric",
  call="call",
  formula="formula",
  #terms="terms",
  data="data.frame",
  offset="NullNum",
  control="list",
  contrasts="NullList",
  theta="matrix",
  vcov="matrix",
  iter="numeric"),
 contains="list" 
)

################################################
# methods defined for cpglm
################################################

# extraction of slots using $
setMethod("$",
    signature(x = "cpglm"),
    function (x, name) 
    {
        slot(x,name)
    }
)

# names to get slot names
setMethod("names",
    signature(x = "cpglm"),
    function (x) 
    {
        return(slotNames(x))
    }
)

# extraction of slots using "[["
setMethod("[[",
    signature(x = "cpglm",i="numeric",j="missing"),
    function (x, i, j, ...) 
    {
	output <- lapply(i, function(y) slot(x,names(x)[y]))
        names(output) <- names(x)[i]
	return(output)
    }
)

setMethod("[[",
    signature(x = "cpglm",i="character",j="missing"),
    function (x, i, j, ...) 
    {
      output <- lapply(1:length(i), function(y) slot(x,i[y]))
      names(output) <- i
      return(output)
    }
)


setMethod("coef",
          signature(object = "cpglm"),
    function (object,...) 
    {
	return(object@coefficients)
    }
)

# variance-covariance matrix as returned by systemfit.
setMethod("vcov",
	signature(object = "cpglm"),
    function (object,...) 
    {
	return(object@vcov)
    }
)


setMethod("residuals",
    signature(object = "cpglm"),
    function (object, ...) 
    {
      return(object@residuals)
    }
)

setMethod("resid",
    signature(object = "cpglm"),
    function (object, ...) 
    {
	return(residuals(object))
    }
)

# generate fitted values on the original scale
setMethod("fitted",
    signature(object = "cpglm"),
    function (object,...) 
    {
      return(object@fitted.values)
    }
)
		
setMethod("fitted.values",
    signature(object = "cpglm"),
    function (object,...) 
    {
      fitted(object)
    }
)


setMethod("show",signature(object = "cpglm"),
	function(object){
		print(summary(object)) 
	}
)

setMethod("summary", signature(object="cpglm"),
	function(object,portfolio=NULL,...){
          
          coef.beta <- coef(object)
          s.err <- sqrt(diag(object@vcov))
          err.beta <- s.err[1:(length(s.err)-2)] 
          tvalue <- coef.beta/err.beta
          dn <- c("Estimate", "Std. Error")             
          pvalue <- 2 * pnorm(-abs(tvalue))
          coef.table <- cbind(coef.beta, err.beta, tvalue, pvalue)  
          dimnames(coef.table) <- list(names(coef.beta), c(dn, 
                "z value", "Pr(>|z|)"))
 
          coef.other <- cbind(object$theta[nrow(object$theta),],s.err)[-(1:length(coef.beta)),]
          dimnames(coef.other) <- list(c("phi","p"),dn) 

          return(list(coef.table,coef.other))

        }
)    
