\name{cpglm-class}
\docType{class}
\alias{cpglm-class}
\alias{$,cpglm-method}
\alias{[[,cpglm,numeric,missing-method}
\alias{[[,cpglm,character,missing-method}
\alias{coef,cpglm-method}
\alias{fitted.values,cpglm-method}
\alias{fitted,cpglm-method}
\alias{names,cpglm-method}
\alias{residuals,cpglm-method}
\alias{resid,cpglm-method}
\alias{show,cpglm-method}
\alias{summary,cpglm-method}
\alias{vcov,cpglm-method}

\title{Class "cpglm" of the result from fitting compound Poisson GLM }
\description{ This class includes the results from calling the \code{\link{cpglm}} function to fit the compound Poisson generalized linear model using the Monte Carlo EM algorithm. Several primitive methods and statistical methods are also created to facilitate the extraction of specific slots and further statistical analysis. }

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("cpglm", ...)}, or they could also be a result of calls from \code{cpglm}.

}
\section{Slots}{
\describe{
    \item{\code{.Data}:}{Object of class \code{"list"}}
    \item{\code{coefficients}:}{Object of class \code{"numeric"}, estimated mean parameters }
    \item{\code{residuals}:}{Object of class \code{"numeric"}, the working residuals, that is the residuals in the final iteration of the IWLS fit }
    \item{\code{fitted.values}:}{Object of class \code{"numeric"}, the fitted mean values, obtained by transforming the linear predictors by the inverse of the link function }
    \item{\code{weights}:}{Object of class \code{"NullNum"}, prior weights }
    \item{\code{df.residual}:}{Object of class \code{"integer"}, residual degrees of freedom: \code{df.null} minus number of parameters }
    \item{\code{df.null}:}{Object of class \code{"integer"}, null degrees of freedom: number of observations }
    \item{\code{y}:}{Object of class \code{"numeric"}, the response variable used }
    \item{\code{call}:}{Object of class \code{"call"}, the matched call }
    \item{\code{formula}:}{Object of class \code{"formula"}, the formula supplied }
    \item{\code{data}:}{Object of class \code{"data.frame"}, the supplied data }
    \item{\code{offset}:}{Object of class \code{"NullNum"}, the offset vector used}
    \item{\code{control}:}{Object of class \code{"list"}, the value of the control argument used }
    \item{\code{contrasts}:}{Object of class \code{"NullList"}, the contrasts used }
    \item{\code{theta}:}{Object of class \code{"matrix"}, a matrix that records the iteration history of all the parameters}
    \item{\code{vcov}:}{Object of class \code{"matrix"}, estimated variance-covariance matrix at convergence }
    \item{\code{iter}:}{Object of class \code{"numeric"}, number of iterations used }
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
    \item{[[}{\code{signature(x = "cpglm", i = "numeric", j = "missing")}:  Method for primitive function \code{"[["}. It extracts the i-th slot of a \code{"cpglm"} object, just as in list. \code{i} could be a vetor. }
    \item{[[}{\code{signature(x = "cpglm", i = "character", j = "missing")}: Method for primitive function \code{"[["}. It extracts the slots of a \code{"cpglm"} object with names in \code{i}, just as in list. \code{i} could be a vetor. }
    \item{coef}{\code{signature(object = "cpglm")}: Method for function \code{coef},  to extract the estimated development matrix. The output is a list.  }
    \item{fitted.values}{\code{signature(object = "cpglm")}:  Method for function \code{fitted.values},  to calculate the fitted values in the orignal triangles. Note that the return value is a list of fitted valued based on the orignal scale, not the model scale which is first divided by \eqn{Y_{i,k}^{\delta/2}}.  }
    \item{fitted}{\code{signature(object = "cpglm")}: Same as \code{fitted.values} in the above. }
    \item{names}{\code{signature(x = "cpglm")}: Method for function \code{names}, which returns the slot names of a \code{"cpglm"} object. }
    \item{residuals}{\code{signature(object = "cpglm")}: Method for function \code{residuals}, to extract residuals  from a \code{cpglm} object.}
    \item{resid}{\code{signature(object = "cpglm")}: Same as \code{residuals}.  }
    \item{show}{\code{signature(object = "cpglm")}: Method for \code{show}. }
    \item{summary}{\code{signature(object = "cpglm")}: Method for function \code{summary}.  }
    \item{vcov}{\code{signature(object = "cpglm")}: Method for function \code{vcov}, to extract the variance-covariance matrix of a \code{"cpglm"} object.}
	 }
}

\author{ Wayne  Zhang \email{actuary_zhang@hotmail.com} }
\seealso{
	See also \code{\link{cpglm}}.   
}

\keyword{classes}
