\name{cplm-class-methods}
\docType{class}
\alias{cplm-class}
\alias{$,cplm-method}
\alias{[,cplm,character,missing,missing-method}
\alias{[,cplm,numeric,missing,missing-method}
\alias{[[,cplm,character,missing-method}
\alias{[[,cplm,numeric,missing-method}
\alias{model.matrix,cplm-method}
\alias{names,cplm-method}
\alias{show,cplm-method}
\alias{terms,cplm-method}
\alias{formula,cplm-method}
\alias{vcov,cplm-method}

\alias{cpglm-class}
\alias{coef,cpglm-method}
\alias{fitted,cpglm-method}
\alias{residuals,cpglm-method}
\alias{resid,cpglm-method}
\alias{AIC,cpglm,missing-method}
\alias{deviance,cpglm-method}
\alias{summary,cpglm-method}
\alias{predict,cpglm-method}

\alias{cpglmm-class}
\alias{print,cpglmm-method}
\alias{summary,cpglmm-method}
\alias{show,cpglmm-method}
\alias{VarCorr,cpglmm-method}
\alias{vcov,cpglmm-method}
\alias{predict,cpglmm-method}
\alias{summary.cpglmm-class}

\alias{bcpglm-class}
\alias{plot,bcpglm,missing-method}
\alias{summary,bcpglm-method}

\alias{bcpglmm-class}

\alias{NullNum-class}
\alias{NullList-class}
\alias{NullFunc-class}

\title{Classes and Methods for a Compound Poisson Linear Model Object}

\description{
Documented here are the class \code{"cplm"} used to represent a compound Poisson linear model object and its derived classes \code{"cpglm"}, \code{"cpglmm"}, \code{"bcpglm"} and \code{"bcpglmm"}. Several primitive methods and statistical methods are created to facilitate the extraction of specific slots and further statistical analysis. \code{"NullNum"},  \code{"NullList"} and \code{"NullFunc"} are virtual classes for \code{c("NULL", "numeric")}, \code{c("NULL","list")} and \code{c("NULL","function")}, respectively. 
}

\section{Objects from the Class}{

\describe{
  \item{\code{"cplm"}}{Objects can be created by calls of the form \code{new("cplm", ...)}.} 
  \item{\code{"cpglm"}}{Objects can be created by calls from \code{new("cpglm", ...)} or \code{cpglm}.}
  \item{\code{"cpglmm"}}{Objects can be created by calls of the form \code{new("cpglmm", ...)}, or a call to \code{cpglmm}.}
  \item{\code{"summary.cpglmm"}}{Objects can be created by calls of the form \code{new("summary.cpglmm", ...)}, or a call to \code{summary} on a \code{cpglmm} object. }
  \item{\code{"bcpglm"}}{Objects can be created by calls from \code{new("bcpglm", ...)} or \code{bcpglm}.}
  \item{\code{"bcpglmm"}}{Objects can be created by calls from \code{new("bcpglmm", ...)} or \code{bcpglmm}.}
  \item{\code{"NullNum"}, \code{"NullList"}, \code{"NullFunc"}}{These are  virtual classes and  no objects may be created from them.}
 } 
}

\section{Slots}{
  The \code{"cplm"} class defines the slots common in all the classes in the \code{cplm} package, and thus the utility methods defined on the \code{"cplm"} class such as \code{[} and \code{names} are also applicable to all of the derived classes. 
  \describe{
    \item{\code{call}:}{the matched call. }
    \item{\code{formula}:}{the formula supplied, class \code{"formula"}}
    \item{\code{contrasts}:}{the contrasts used, class \code{"NullList"} }
    \item{\code{link.power}:}{index of power link function, class \code{"numeric"}. See \code{\link[statmod]{tweedie}}.}
    \item{\code{model.frame}:}{the data frame used. class \code{"data.frame"}. }
    \item{\code{inits}:}{initial values used, class \code{"NullList"}.}
  }
  

 The \code{"cpglm"} class extends \code{"cplm"} directly and is a variant of \code{\link{glm}} with additional functionality to estimate the index parameter. So most of the slots have the same definition as those in \code{\link{glm}}. But for the Monte Carlo EM algorithm (\code{method="MCEM"}), some of these slots have different meanings. The following slots are in addition to those in \code{"cplm"}:
\describe{
    \item{\code{coefficients}:}{estimated mean parameters, class \code{"numeric"}. }
    \item{\code{residuals}:}{the working residuals, that is the residuals in the final iteration of the IWLS fit, class \code{"numeric"}}
    \item{\code{fitted.values}:}{the fitted mean values, obtained by transforming the linear predictors by the inverse of the link function, class \code{"numeric"} }
     \item{\code{linear.predictors}:}{the fitted linear predictors, class \code{"numeric"}}    
    \item{\code{weights}:}{working weights from the last iteration of the iterative least square, class \code{"numeric"}}
    \item{\code{df.residual}:}{residual degrees of freedom, class \code{"integer"}}
    \item{\code{deviance}:}{up to a constant, minus twice the maximized log-likelihood. Where sensible, the constant is chosen so that a saturated model has deviance zero. This is computed using \code{\link[tweedie]{tweedie.dev}}.}
    \item{\code{aic}:}{a version of Akaike's Information Criterion, minus twice the maximized log-likelihood plus twice the number of mean parameters. This is computed using the tweedie density approximation as in \code{\link[tweedie]{dtweedie}}. }
    \item{\code{offset}:}{the offset vector used, class \code{"NullNum"},} 
    \item{\code{prior.weights}:}{the weights initially supplied, a vector of \code{1}s if none were, class \code{"NullNum"}}        
    \item{\code{y}:}{the response vector used.}
    \item{\code{control}:}{the value of the control argument used, class \code{"list"} }
    \item{\code{p}:}{the maximum likelihood estimate of the index parameter.}    
    \item{\code{phi}:}{the maximum likelihood estimate of the dispersion parameter.}    
    \item{\code{vcov}:}{estimated variance-covariance matrix, class \code{"matrix"}}    
    \item{\code{iter}:}{the number of Fisher's scoring iterations in the GLM.}
    \item{\code{converged}:}{indicating whether the algorithm has converged, class \code{"logical"}.}
    \item{\code{na.action}:}{method of handling \code{NA}, class \code{"NullFunc"}.}
}    
    

  The class \code{"cpglmm"} extends \code{"\linkS4class{mer}"} and \code{"cplm"}, and it has the following additional slots (see \code{"\linkS4class{mer}"} for details of the slots not documented here):
  \describe{
    \item{\code{p}:}{estimated value of the index parameter, class \code{"numeric"}  }
    \item{\code{phi}:}{estimated value of the dispersion parameter, class \code{"numeric"} }
    \item{\code{bound.p}:}{the specified bound of the index parameter, class \code{"numeric"}   }
    \item{\code{vcov}:}{estimated variance-covariance matrix, class \code{"matrix"}}
    }

  The \code{"summary.cpglmm"} class \emph{contains} the \code{"cpglmm"}
  class and has the following additional slots:
  \describe{
    \item{\code{methTitle}:}{character string specifying a method title}
    \item{\code{logLik}:}{the same as \code{logLik(object)}.}
    \item{\code{ngrps}:}{the number of levels per grouping factor in the
      \code{flist} slot.}
    \item{\code{sigma}:}{the scale factor for the variance-covariance estimates}
    \item{\code{coefs}:}{the matrix of estimates, standard errors,
      etc. for the fixed-effects coefficients}
    \item{\code{REmat}:}{the formatted Random-Effects matrix}
    \item{\code{AICtab}:}{A named vector of values of AIC, BIC, log-likelihood
      and deviance}
  }    
  
  The class \code{"bcpglm"} extends the class \code{"cplm"}, and it has the following additional slots: 
  \describe{
    \item{\code{n.chains}:}{number of Markov chains, class \code{"integer"}}
    \item{\code{n.iter}:}{number of total iterations per chain, class \code{"integer"}}
    \item{\code{n.burnin}:}{length of burn in, i.e. number of iterations to discard at the beginning, class \code{"integer"}}
    \item{\code{n.thin}:}{thinning rate. class \code{"integer"}}
    \item{\code{n.sims}:}{the approximate number of simulations to keep after thinning, class \code{"integer"}}
    \item{\code{sims.list}:}{Object of class \code{"mcmc.list"}, a list of \code{n.chains} \code{mcmc} objects, each \code{mcmc} object storing the simulation result from a Markov chain. See \code{\link[coda]{mcmc}} and \code{\link[coda]{mcmc.convert}}. Since this is an \code{"mcmc.list"} object, most methods defined in the package \code{coda} can be directly applied to it. }
    \item{\code{summary}:}{Object of class \code{"summary.mcmc"}, two sets of summary statistics for each variable: Mean, standard deviation, naive standard error of the mean (ignoring autocorrelation of the chain) and time-series standard error based on an estimate of the spectral density at 0. See \code{\link[coda]{summary.mcmc}}.}
    \item{\code{prop.var}:}{a named list of proposal variance-covariance matrix used in the Metropolis-Hasting update}
}

  The class \code{"bcpglmm"} extends the class \code{"bcpglm"}, and it has the following additional slots: 
  \describe{
    \item{\code{Zt}:}{The transpose of model matrix for the random effects, stored as a compressed column-oriented sparse matrix (class \code{"dgCMatrix"}).}
    \item{\code{flist}:}{The list of grouping factors for the random effects.}
    \item{\code{prop.var}:}{a named list of proposal variance-covariance matrix used in the Metropolis-Hasting update}
 }
  
}


\section{Extends}{

Class \code{"cpglm"} extends class \code{"\linkS4class{cplm}"}, directly.

Class \code{"cpglmm"} extends class \code{"\linkS4class{mer}"} and class \code{"\linkS4class{cplm}"}, directly.

Class \code{"summary.cpglmm"}  extends class \code{"\linkS4class{cpglmm}"}, directly; 
class \code{"\linkS4class{mer}"}, by class "cpglmm", distance 2;
class \code{"\linkS4class{cplm}"}, by class "cpglmm", distance 2.

Class \code{"bcpglm"} extends class \code{"\linkS4class{cplm}"}, directly.

Class \code{"bcpglmm"} extends class \code{"\linkS4class{bcpglm}"}, directly;
Class \code{"\linkS4class{cplm}"}, by class "bcpglm", distance 2.

}


\section{Methods}{
  The following methods are defined for the class \code{"cplm"}, which are also applicable to all of the derived classes:
  \describe{
    \item{$}{\code{signature(x = "cplm")}:  extract a slot of \code{x} with a specified slot name, just as in list. }
    \item{[[}{\code{signature(x = "cplm", i = "numeric", j = "missing")}:  extract the i-th slot of a \code{"cpglm"} object, just as in list. }
    \item{[[}{\code{signature(x = "cplm", i = "character", j = "missing")}: extract the slots of a \code{"cpglm"} object with names in \code{i}, just as in list.}
    \item{[}{\code{signature(x = "cplm", i = "numeric", j = "missing", drop="missing")}:  extract the i-th slot of a \code{"cpglm"} object, just as in list. \code{i} could be a vetor. }
    \item{[}{\code{signature(x = "cplm", i = "character", j = "missing", drop="missing")}:  extract the slots of a \code{"cpglm"} object with names in \code{i}, just as in list. \code{i} could be a vetor. }    
    \item{names}{\code{signature(x = "cplm")}:  return the slot names. }
    \item{terms}{\code{signature(x = "cplm")}: extract the \code{terms} object from the model frame.  See \code{\link[stats]{terms}}.}
    \item{formula}{\code{signature(x = "cplm")}: extract the \code{formula} slot.  See \code{\link[stats]{formula}}.}    
    \item{model.matrix}{\code{signature(object = "cplm")}: extract the design matrix.  }
    \item{show}{\code{signature(object = "cplm")}: method for \code{show}. }
    \item{vcov}{\code{signature(object = "cplm")}: extract the variance-covariance matrix of a \code{"cplm"} object.}
 }
 
 The following methods are defined for the \code{"cpglm"} class:
  \describe{    
    \item{coef}{\code{signature(object = "cpglm")}:  extract the estimated coefficients.}
    \item{fitted}{\code{signature(object = "cpglm")}: return the fitted values. }
    \item{residuals}{\code{signature(object = "cpglm")}: extract residuals  from a \code{cpglm} object. You can also specify a \code{type} argument to indicate the type of residuals to be computed. See \code{\link[stats]{glm.summaries}}.}
    \item{resid}{\code{signature(object = "cpglm")}: same as \code{residuals}.}   
    \item{AIC}{\code{signature(object = "cpglm",k="missing")}: extract the AIC information from the \code{"cpglm"} object.  See \code{\link[stats]{AIC}}.}    
    \item{deviance}{\code{signature(object = "cpglm")}: extract the deviance from the \code{"cpglm"} object.  See \code{\link[stats]{deviance}}.}
    \item{summary}{\code{signature(object = "cpglm")}: for the profiled likelihood approach, this is the same as \code{\link[stats]{glm.summaries}} except that both the dispersion and the index parameter are estimated using maximum likelihood estimation. For \code{method="MCEM"}, a Z-test is performed that relies on the asymptotic normality of the MLE estimates.}    
    \item{predict}{\code{signature(object = "cpglm")}: generate predictions for new data sets}

   }
   
   
  For the class \code{"cpglmm"}, most methods defined for \code{"\linkS4class{mer}"} can be applied directly. The following are written specifically for \code{"cpglmm"}: 
  \describe{
    \item{print}{\code{signature(x = "cpglmm")}:  print the object }
    \item{summary}{\code{signature(object = "cpglmm")}: summary results}
    \item{predict}{\code{signature(object = "cpglmm")}: generate predictions for new data sets}
    \item{VarCorr}{\code{signature(x = "cpglmm")}: estimation for the variance components }
    \item{vcov}{\code{signature(object = "cpglmm")}: variance-covariance matrix for fixed effects }
   }
   
  The following methods are available for the class \code{"bcpglm"}: 
  \describe{
    \item{plot}{\code{signature(x = "bcpglm", y = "missing")}: summarize the \code{"bcpglm"} object with a trace of the sampled output and a density estimate for each variable in the chain. See \code{\link[coda]{plot.mcmc}}. }
    \item{summary}{\code{signature(object = "bcpglm")}:  produce two sets of summary statistics. See \code{\link[coda]{summary.mcmc}}. }
   }
    
}

\author{ Wayne  Zhang \email{actuary_zhang@hotmail.com} }
\seealso{
  See also \code{\link{cpglm}}, \code{\link{cpglmm}}, \code{\link{bcpglm}}, \code{\link{bcpglmm}}, \code{\link[stats]{glm}} and \code{\link[lme4]{mer-class}}.   
}

\keyword{classes}
