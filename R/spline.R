###########################################################
# Functions used in fitting additive models in cpglmm     #
###########################################################


#######################
# get model frame and factor list 
# for cpglmm with smoothing terms
#######################
frFL <- function (formula, data, family, control = list(),  
    verbose, weights, offset, contrasts, basisGenerators, bySetToZero = T) 
{
    call <- match.call()
    formula <- eval(call$formula)
    tf <- terms.formula(formula, specials = eval(call$basisGenerators, 
        parent.frame(2)))
    f.ind <- unlist(attr(tf, "specials"))
    n.f <- length(f.ind)
    rhs <- amer:::safeDeparse(formula[[3]])
    fctterm <- fct <- vector(mode = "list", length = n.f)
    for (i in 1:n.f) 
      fctterm[[i]] <- attr(tf, "variables")[[f.ind[i] + 1]]
    fct <- lapply(fctterm, eval, envir = data, enclos = parent.frame(2))
    for (i in seq_along(fct)) 
      fct[[i]] <- amer:::expandBasis(fct[[i]], eval(attr(fct[[i]], "call")$by, data), 
                  eval(attr(fct[[i]], "call")$varying, data), bySetToZero)
    names(fct) <- names(fctterm) <- paste("f.", lapply(fct, 
            function(x) {
                paste(as.character(attr(x, "call")$x), ifelse(!is.null(eval(attr(x, 
                  "call")$varying, data)), paste("X", deparse(attr(x, 
                  "call")$varying), sep = ""), ""), ifelse(eval(attr(x, 
                  "call")$allPen), paste(".", deparse(attr(x, 
                  "call")$by), sep = ""), ""), sep = "")
            }), sep = "")
    rhs <- amer:::subFcts(rhs, fctterm, fct, data)
    data <- amer:::expandMf(data, fct)
    call[[1]] <- as.name("lmer")
    call$doFit <- FALSE
    call$data <- as.name("data")
    call$formula <- as.formula(paste(formula[[2]], "~", rhs))
    call["basisGenerators"] <- NULL
    m <- eval(call, data)
    #    m$fr$mf <- data
    m <- amer:::subAZ(m, fct)
    fctterm <- lapply(fct, function(x) attr(x, "call"))
    return(list(m = m, fct = fct, fctterm = fctterm))
}

############################
# function for 2d splines (from Ngo and Wand)
############################
# compute thin plate spline covariance function  
tps.cov <- function(r) {
  r <- as.matrix(r)
  num.row <- nrow(r)
  num.col <- ncol(r)
  r <- as.vector(r)
  nzi <- (1:length(r))[r!=0]
  ans <- rep(0,length(r))
  ans[nzi] <- r[nzi]*r[nzi]*log(abs(r[nzi]))
  if (num.col>1) 
    ans <- matrix(ans,num.row,num.col)
  return(ans)
}

sp2d <- function(x1, x2, k = max(20,min(length(x1)/4,150)), 
                 by = NULL, allPen = FALSE, varying = NULL, 
                 diag = FALSE, knots1 = quantile(x1, probs = 1:k/(k+1)),
                 knots2 = quantile(x1, probs = 1:k/(k+1))) {
  call <- as.list(amer:::expand.call())
  knots1 <- eval(knots1)
  knots2 <- eval(knots2)
  k <- eval(k)
  # design matrix for fixed effects
  X <- cbind(x1,x2)
  knots <- cbind(knots1,knots2)
  dist.mat <- matrix(0,k,k)
  dist.mat[lower.tri(dist.mat)] <- dist(knots)
  dist.mat <- dist.mat + t(dist.mat)
  Omega <- tps.cov(dist.mat)
  diffs.1 <- outer(x1,knots1,"-")
  diffs.2 <- outer(x2,knots2,"-")
  dists <- sqrt(diffs.1*diffs.1+diffs.2*diffs.2)
  svd.Omega <- svd(Omega)
  sqrt.Omega <- t(svd.Omega$v %*% (t(svd.Omega$u) * sqrt(svd.Omega$d)))
  # design matrix for random effects
  Z <- t(solve(sqrt.Omega,t(tps.cov(dists))))
  res <- list(X=X, Z=Z, knots=knots)
  attr(res, "call") <- as.call(call)  
  return(res)
}


