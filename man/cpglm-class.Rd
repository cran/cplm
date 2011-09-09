\name{cpglm-class}
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


\title{Class "cpglm" of the result from fitting compound Poisson GLM }
\description{This class includes the result from calling the \code{\link{cpglm}} function to fit the compound Poisson generalized linear model. Several primitive methods and statistical methods are created to facilitate the extraction of specific slots and further statistical analysis. }

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("cpglm", ...)}, or they could also be a result of calls from \code{cpglm}.
}

\section{Slots}{
The function \code{cpglm} is a variant of the generalized linear model that allows the estimation of the index parameter from the data, and thus most of the slots in the result have the same definition as those in \code{\link{glm}}. But for the Monte Carlo EM algorithm (\code{method="MCEM"}), some of these slots have different meanings:
\describe{
    \item{\code{.Data}:}{Object of class \code{"list"}.}
    \item{\code{coefficients}:}{Object of class \code{"numeric"}, estimated mean parameters. }
    \item{\code{residuals}:}{Object of class \code{"numeric"}, the working residuals, that is the residuals in the final iteration of the IWLS fit.}
    \item{\code{fitted.values}:}{Object of class \code{"numeric"}, the fitted mean values, obtained by transforming the linear predictors by the inverse of the link function. }
     \item{\code{linear.predictors}:}{Object of class \code{"numeric"}, the fitted linear predictors}    
    \item{\code{weights}:}{Object of class \code{"numeric"}, working weights from the last iteration of the iterative least square.}
    \item{\code{df.residual}:}{Object of class \code{"integer"}, residual degrees of freedom.}
    \item{\code{deviance}:}{Object of class \code{"numeric"}, up to a constant, minus twice the maximized log-likelihood. Where sensible, the constant is chosen so that a saturated model has deviance zero. This is computed using \code{\link[tweedie]{tweedie.dev}}.}
    \item{\code{aic}:}{Object of class \code{"numeric"}, a version of Akaike's Information Criterion, minus twice the maximized log-likelihood plus twice the number of parameters. The number of parameters is the number of mean parameters plus two, as both the dispersion parameter and the index parameter are estimated. This is computed using the tweedie density approximation as in \code{\link[tweedie]{dtweedie}}. }
    \item{\code{offset}:}{Object of class \code{"NullNum"}, the offset vector used.} 
    \item{\code{prior.weights}:}{Object of class \code{"NullNum"}, the weights initially supplied, a vector of \code{1}s if none were.}    
    \item{\code{call}:}{Object of class \code{"call"}, the matched call. }
    \item{\code{formula}:}{Object of class \code{"formula"}, the formula supplied. }
    \item{\code{data}:}{Object of class \code{"data.frame"}, the supplied data. }
    \item{\code{control}:}{Object of class \code{"list"}, the value of the control argument used. }
    \item{\code{contrasts}:}{Object of class \code{"NullList"}, the contrasts used. }
    \item{\code{theta}:}{Object of class \code{"numeric"}, a vecotr that records the maximized values of all parameters.}    
    \item{\code{theta.all}:}{Object of class \code{"matrix"}, a matrix that records the iteration history of all the parameters. This is meaningless if \code{method="profile"}.}
    \item{\code{p}:}{Object of class \code{"numeric"}, the maximum likelihood estimate of the index parameter.}    
    \item{\code{phi}:}{Object of class \code{"numeric"}, the maximum likelihood estimate of the dispersion parameter.}    
    \item{\code{vcov}:}{Object of class \code{"matrix"}, estimated variance-covariance matrix.}    
    \item{\code{iter}:}{Object of class \code{"integer"}. For \code{method="profile"}, this is the number of Fisher's scoring iterations in the GLM, and for \code{method="MCEM"}, this is the number of iterations used in the MCEM algorithm.}
    \item{\code{converged}:}{Object of class \code{"logical"}, indicating whether the algorithm has converged.}
    \item{\code{method}:}{Object of class \code{"character"}, the method used in \code{cplm}.}
    \item{\code{y}:}{Object of class \code{"numeric"}, the response vector used.}    
    \item{\code{link.power}:}{Object of class \code{"numeric"}, index of power link function. See \code{\link[statmod]{tweedie}}.}
    \item{\code{na.action}:}{Object of class \code{"NullFunc"}, a function which indicates what should happen when the data contain \code{NA}s. See \code{\link[stats]{glm}}.}
    \item{\code{model.frame}}{Object of class \code{"data.frame"}, the data frame used in \code{cpglm}.}
  }
}

\section{Extends}{
Class \code{"\linkS4class{list}"}, from data part.
Class \code{"\linkS4class{vector}"}, by class "list", distance 2.
Class \code{"\linkS4class{NullList}"}, by class "list", distance 2.
}

\section{Methods}{
  \describe{
    \item{$}{\code{signature(x = "cpglm")}: Method for primitive function \code{"$"}. It extracts a slot of \code{x} with a specified slot name, just as in list. }
    \item{[[}{\code{signature(x = "cpglm", i = "numeric", j = "missing")}:  Method for primitive function \code{"[["}. It extracts the i-th slot of a \code{"cpglm"} object, just as in list. }
    \item{[[}{\code{signature(x = "cpglm", i = "character", j = "missing")}: Method for primitive function \code{"[["}. It extracts the slots of a \code{"cpglm"} object with names in \code{i}, just as in list.}
    \item{[}{\code{signature(x = "cpglm", i = "numeric", j = "missing", drop="missing")}:  Method for primitive function \code{"["}. It extracts the i-th slot of a \code{"cpglm"} object, just as in list. \code{i} could be a vetor. }
    \item{[}{\code{signature(x = "cpglm", i = "character", j = "missing", drop="missing")}: Method for primitive function \code{"["}. It extracts the slots of a \code{"cpglm"} object with names in \code{i}, just as in list. \code{i} could be a vetor. }    
    \item{coef}{\code{signature(object = "cpglm")}: Method for function \code{coef},  to extract the estimated development matrix. The output is a list.  }
    \item{fitted.values}{\code{signature(object = "cpglm")}:  Method for function \code{fitted.values},  to calculate the fitted values in the orignal triangles. Note that the return value is a list of fitted valued based on the orignal scale, not the model scale which is first divided by \eqn{Y_{i,k}^{\delta/2}}.  }
    \item{fitted}{\code{signature(object = "cpglm")}: Same as \code{fitted.values} in the above. }
    \item{names}{\code{signature(x = "cpglm")}: Method for function \code{names}, which returns the slot names of a \code{"cpglm"} object. }
    \item{residuals}{\code{signature(object = "cpglm")}: Method for function \code{residuals}, to extract residuals  from a \code{cpglm} object. You can also specify a \code{type} argument to indicate the type of residuals to be computed. See \code{\link[stats]{glm.summaries}}.}
    \item{resid}{\code{signature(object = "cpglm")}: Same as \code{residuals}.}
    \item{df.residual}{\code{signature(object = "cpglm")}: Extract the residual degree of freedome from the \code{"cpglm"} object. See \code{\link[stats]{df.residual}}.  } 
    \item{AIC}{\code{signature(object = "cpglm",k="missing")}: Extract the AIC information from the \code{"cpglm"} object.  See \code{\link[stats]{AIC}}.}    
    \item{deviance}{\code{signature(object = "cpglm")}: Extract the deviance from the \code{"cpglm"} object.  See \code{\link[stats]{deviance}}.}
    \item{terms}{\code{signature(x = "cpglm")}: Extract the \code{terms} object from the \code{"cpglm"} object.  See \code{\link[stats]{terms}}.}
    \item{formula}{\code{signature(x = "cpglm")}: Extract the \code{formula} object from the \code{"cpglm"} object.  See \code{\link[stats]{formula}}.}    
    \item{model.matrix}{\code{signature(object = "cpglm")}: Extract the design matrix from the \code{"cpglm"} object.  }
    \item{show}{\code{signature(object = "cpglm")}: Method for \code{show}. }
    \item{summary}{\code{signature(object = "cpglm")}: Method for function \code{summary}. For the profiled likelihood approach, this is the same as \code{\link[stats]{glm.summaries}} excepte that both the dispersion and the index parameter are estimated using maximum likelihood estimation. For \code{method="MCEM"}, a Z-test is performed that relies on the asymptotic normality of the MLE estimates.}
    \item{vcov}{\code{signature(object = "cpglm")}: Method for function \code{vcov}, to extract the variance-covariance matrix of a \code{"cpglm"} object.}
	 }
}

\author{ Wayne  Zhang \email{actuary_zhang@hotmail.com} }
\seealso{
	See also \code{\link{cpglm}}, \code{\link[stats]{glm}} and \code{\link[stats]{glm.summaries}}.   
}

\keyword{classes}
