#######################################################
##             compound Poisson GLM                  ##
## Author: Wayne Zhang, actuary_zhang@hotmail.com    ##
#######################################################

cpglm <- function(formula, link = "log", data, weights, offset, 
                  subset, na.action = NULL, contrasts = NULL, 
                  control = list(), chunksize = 0, ...) {

  call <- match.call()  
  if (missing(data)) 
    data <- environment(formula) 
  # use bigglm for big data sets 
  if (chunksize) {
    ans <-  cpglm.profile.bigglm(call, data)
  } else {
    fr <- cpglm.mf(call, contrasts)
    link.power <- make.link.power(link)
    control <- do.call("cpglm.control", control)
    ans <- cpglm.profile(fr, link.power, control)  
  }
  ans@formula <- formula 
  ans@call <- call
  ans@na.action <- na.action
  ans@contrasts <- contrasts
  return(ans)
}


# function to implement the profile likelihood approach 
# Use cpglm.profile.bigglm for big data sets  
cpglm.profile <- function(fr, link.power = 0, control = list()){
  
  # return profile loglikelihood
  # log.phi indicates whether phi (as in parm[1]) is on the log scale
  llik_profile <- function(parm, log.phi = TRUE){
    phi <- ifelse(log.phi, exp(parm[1]), parm[1])
    p <- parm[2]
    fit2 <- cpglm.fit(fr, p = p, link.power) 
    - sum(log(dtweedie.series(fr$Y, p, fit2$fitted.values, phi)))    
  }
  
  # generate starting values for phi
  pstart <- 1.5
  fit <- cpglm.fit(fr, pstart, link.power)
  phistart <- sum((fit$weights * fit$residuals^2)) / fit$df.residual
  parm <- c(log(phistart), pstart)
  
  # optimize the profile loglikelihood
  opt_ans <- nlminb(parm,llik_profile, gradient = NULL, 
                    hessian = NULL, log.phi = TRUE, 
                    lower = c(-Inf, control$bound.p[1]),
                    upper = c(Inf, control$bound.p[2]),
                    control = list(trace = control$trace,
                                   iter.max = control$max.iter))
  if (opt_ans$convergence) warning(opt_ans$message)
  # parameter values at optima
  p.max <- opt_ans$par[2]
  phi.max <- exp(opt_ans$par[1])
  
  # fit glm using the optimized index parameter
  fit <- cpglm.fit(fr, p.max, link.power)                  
  class(fit) <- "glm" 
  
  # compute vcov for p and phi      
  pm <- c(phi.max, p.max) 
  hs <- hess(pm,llik_profile, log.phi = FALSE)
  dimnames(hs) <- list(c("phi","p"),c("phi","p"))  
  vc <- vcov(fit)
  attr(vc,"phi_p") <- svd.inv(hs)
    
  # return results
  out <- new("cpglm", 
             coefficients = fit$coefficients, residuals = fit$residuals,
             fitted.values = fit$fitted.values, weights = fit$weights,
             linear.predictors =  fit$linear.predictors,
             df.residual = as.integer(fit$df.residual),
             deviance = fit$deviance, call = call("foo"),
             aic = dtweedie.nlogl(fr$Y, fitted(fit), phi.max, p.max) + 2 * fit$rank,           
             formula = ~ 1, control = control, contrasts = NULL,
             p = p.max, phi = phi.max, iter = fit$iter, converged = fit$converged,
             link.power= link.power, model.frame = fr$mf, na.action = NULL,
             offset = fr$offset, prior.weights = fit$prior.weights, y = fr$Y,
             inits = NULL, vcov = vc)
  return(out)  
}               


# function to implement the  profile likelihood approach for big data sets
cpglm.profile.bigglm <- function(call, data){
  # restruct the call
  if (is.null(control <- call$control)) 
    control <- list()  
  if (is.null(link <- call$link)) 
    link <- "log"
  control <- do.call("cpglm.control", eval(control))
  link.power <- make.link.power(eval(link))
  mc <- match(c("link", "control", "family"), names(call), 0L)
  call <- call[-c(1, mc[mc>0])]
  pstart <- 1.5
  call$family <- tweedie(var.power = pstart, link.power = link.power)
  # default maxit in bigglm: the default seems not to be enough 
  call$maxit <- ifelse(is.null(call$maxit), 50, call$maxit)
  Y <- eval(call$formula[[2]], envir = data)
  
  # profile loglikelihood 
  llik_profile <- function(parm, log.phi = TRUE){
    phi <- ifelse(log.phi, exp(parm[1]), parm[1])
    p <- parm[2]
    call$family <- tweedie(var.power = p, link.power = link.power)
    fit2 <- do.call("bigglm", as.list(call))
    - sum(log(dtweedie.series(Y, p, fitted(fit2, data)$fitted.values, phi)))    
  }
  
  # generate starting values for phi  
  fit <- do.call("bigglm", as.list(call))  
  phistart <- fit$qr$ss / fit$df.resid 
  parm <- c(log(phistart), pstart)
  
  # optimize the profiled loglikelihood
  opt_ans <- nlminb(parm, llik_profile, gradient = NULL,
                    hessian = NULL, log.phi = TRUE,
                    lower = c(-Inf, control$bound.p[1]),
                    upper = c(Inf, control$bound.p[2]),
                    control = list(trace = control$trace,
                                   iter.max = control$max.iter))
  if (opt_ans$convergence) warning(opt_ans$message)
  p.max <- opt_ans$par[2]
  phi.max <- exp(opt_ans$par[1])

  # fit glm using the optimized index parameter
  call$family <- tweedie(var.power = p.max, link.power = link.power)
  fit <- do.call("bigglm", as.list(call))
  fs <- fitted(fit, data)
  
  # compute vcov for p and phi      
  pm <- c(phi.max, p.max) 
  hs <- hess(pm, llik_profile, log.phi = FALSE)
  dimnames(hs) <- list(c("phi","p"), c("phi","p"))  
  vc <- vcov(fit)
  attr(vc,"phi_p") <- svd.inv(hs)
    
  # return results
  out <- new("cpglm",  coefficients = coef(fit), 
             residuals = fs$residuals, fitted.values = fs$fitted.values,
             linear.predictors = fs$linear.predictors,
             weights = fs$weights, df.residual = as.integer(fit$df.resid),
             aic = dtweedie.nlogl(Y, fs$fitted.values, phi.max, p.max) + 
                    2 * (fit$n - fit$df.resid),           
             deviance = fit$deviance, call = call("foo"),
             formula = ~ 1, control = control,
             contrasts = NULL, p = p.max, phi = phi.max, 
             iter = fit$iterations, converged = fit$converged,
             link.power= link.power, model.frame = as.data.frame(NULL),
             na.action = NULL, offset = fs$offset,
             prior.weights = fs$prior.weights, y = Y,
             inits = NULL, vcov = vc)
  return(out)
}

  

