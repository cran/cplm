\name{cpglm-class-method}
\docType{class}
\alias{cpglm-class}
\alias{$,cpglm-method}
\alias{[[,cpglm,numeric,missing-method}
\alias{[[,cpglm,character,missing-method}
\alias{[,cpglm,numeric,missing,missing-method}
\alias{[,cpglm,character,missing,missing-method}
\alias{coef,cpglm-method}
\alias{fitted.values,cpglm-method}
\alias{fitted,cpglm-method}
\alias{names,cpglm-method}
\alias{residuals,cpglm-method}
\alias{resid,cpglm-method}
\alias{show,cpglm-method}
\alias{summary,cpglm-method}
\alias{AIC,cpglm,missing-method}
\alias{deviance,cpglm-method}
\alias{model.matrix,cpglm-method}
\alias{terms,cpglm-method}
\alias{formula,cpglm-method}
\alias{df.residual,cpglm-method}
\alias{vcov,cpglm-method}


\title{Representation of a compound Poisson GLM object}
\description{This class is used to repesent a compound Poisson GLM object, usually a result from calling the \code{\link{cpglm}} function. Several primitive methods and statistical methods are created to facilitate the extraction of specific slots and further statistical analysis. }

\section{Objects from the Class}{
Objects can be created by calls from \code{new("cpglm", ...)} or \code{cpglm}.
}

\section{Slots}{
Class \code{"cpglm"} is a variant of \code{\link{glm}} with additional functionality to estimate the index parameter, and thus most of the slots have the same definition as those in \code{\link{glm}}. But for the Monte Carlo EM algorithm (\code{method="MCEM"}), some of these slots have different meanings:
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
    \item{\code{call}:}{the matched call. }
    \item{\code{formula}:}{the formula supplied, class \code{"formula"} }
    \item{\code{data}:}{the supplied data, class \code{"data.frame"} }
    \item{\code{control}:}{the value of the control argument used, class \code{"list"} }
    \item{\code{contrasts}:}{the contrasts used, class \code{"NullList"} }
    \item{\code{theta}:}{a vecotr that records the maximized values of all parameters, class \code{"numeric"}}    
    \item{\code{theta.all}:}{a matrix (class \code{"matrix"})that records the iteration history of all the parameters. This is meaningless if \code{method="profile"}.}
    \item{\code{p}:}{the maximum likelihood estimate of the index parameter.}    
    \item{\code{phi}:}{the maximum likelihood estimate of the dispersion parameter.}    
    \item{\code{vcov}:}{estimated variance-covariance matrix, class \code{"matrix"}}    
    \item{\code{iter}:}{Object of class \code{"integer"}. For \code{method="profile"}, this is the number of Fisher's scoring iterations in the GLM, and for \code{method="MCEM"}, this is the number of iterations used in the MCEM algorithm.}
    \item{\code{converged}:}{indicating whether the algorithm has converged, class \code{"logical"}.}
    \item{\code{method}:}{the method used in \code{cplm}, class \code{"character"}}
    \item{\code{y}:}{the response vector used.}    
    \item{\code{link.power}:}{index of power link function, class \code{"numeric"}. See \code{\link[statmod]{tweedie}}.}
    \item{\code{na.action}:}{Object of class \code{"NullFunc"}, a function which indicates what should happen when the data contain \code{NA}s. See \code{\link[stats]{glm}}.}
    \item{\code{model.frame}}{the data frame used in \code{cpglm}.}
  }
}


\section{Methods}{
  \describe{
    \item{$}{\code{signature(x = "cpglm")}:  extract a slot of \code{x} with a specified slot name, just as in list. }
    \item{[[}{\code{signature(x = "cpglm", i = "numeric", j = "missing")}:  extract the i-th slot of a \code{"cpglm"} object, just as in list. }
    \item{[[}{\code{signature(x = "cpglm", i = "character", j = "missing")}: extract the slots of a \code{"cpglm"} object with names in \code{i}, just as in list.}
    \item{[}{\code{signature(x = "cpglm", i = "numeric", j = "missing", drop="missing")}:  extract the i-th slot of a \code{"cpglm"} object, just as in list. \code{i} could be a vetor. }
    \item{[}{\code{signature(x = "cpglm", i = "character", j = "missing", drop="missing")}:  extract the slots of a \code{"cpglm"} object with names in \code{i}, just as in list. \code{i} could be a vetor. }    
    \item{coef}{\code{signature(object = "cpglm")}:  extract the estimated coefficients.}
    \item{fitted.values}{\code{signature(object = "cpglm")}:  return the fitted values.  }
    \item{fitted}{\code{signature(object = "cpglm")}: same as \code{fitted.values} in the above. }
    \item{names}{\code{signature(x = "cpglm")}:  return the slot names of a \code{"cpglm"} object. }
    \item{residuals}{\code{signature(object = "cpglm")}: extract residuals  from a \code{cpglm} object. You can also specify a \code{type} argument to indicate the type of residuals to be computed. See \code{\link[stats]{glm.summaries}}.}
    \item{resid}{\code{signature(object = "cpglm")}: same as \code{residuals}.}
    \item{df.residual}{\code{signature(object = "cpglm")}: extract the residual degree of freedome from the \code{"cpglm"} object. See \code{\link[stats]{df.residual}}.  } 
    \item{AIC}{\code{signature(object = "cpglm",k="missing")}: extract the AIC information from the \code{"cpglm"} object.  See \code{\link[stats]{AIC}}.}    
    \item{deviance}{\code{signature(object = "cpglm")}: extract the deviance from the \code{"cpglm"} object.  See \code{\link[stats]{deviance}}.}
    \item{terms}{\code{signature(x = "cpglm")}: extract the \code{terms} object from the \code{"cpglm"} object.  See \code{\link[stats]{terms}}.}
    \item{formula}{\code{signature(x = "cpglm")}: extract the \code{formula} object from the \code{"cpglm"} object.  See \code{\link[stats]{formula}}.}    
    \item{model.matrix}{\code{signature(object = "cpglm")}: extract the design matrix from the \code{"cpglm"} object.  }
    \item{show}{\code{signature(object = "cpglm")}: method for \code{show}. }
    \item{summary}{\code{signature(object = "cpglm")}: for the profiled likelihood approach, this is the same as \code{\link[stats]{glm.summaries}} except that both the dispersion and the index parameter are estimated using maximum likelihood estimation. For \code{method="MCEM"}, a Z-test is performed that relies on the asymptotic normality of the MLE estimates.}
    \item{vcov}{\code{signature(object = "cpglm")}: extract the variance-covariance matrix of a \code{"cpglm"} object.}
	 }
}

\author{ Wayne  Zhang \email{actuary_zhang@hotmail.com} }
\seealso{
	See also \code{\link{cpglm}}, \code{\link[stats]{glm}} and \code{\link[stats]{glm.summaries}}.   
}

\keyword{classes}
