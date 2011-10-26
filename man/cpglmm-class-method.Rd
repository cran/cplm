\name{cpglmm-class-method}
\docType{class}
\alias{cpglmm-class}
\alias{print,cpglmm-method}
\alias{show,cpglmm-method}
\alias{sigma,cpglmm-method}
\alias{summary,cpglmm-method}
\alias{VarCorr,cpglmm-method}
\alias{vcov,cpglmm-method}
\alias{summary.cpglmm-class}

\title{Representaion of the Compound Poisson Mixed Models}
\description{
Class \code{"cpglmm"}  represents compound Poisson generalized
  linear mixed-effects models, and the \code{"summary.cpglmm"} class represents the summary
  of the object.
}

\section{Objects from the Class}{
Objects of class \code{"cpglmm"} can be created by calls of the form \code{new("cpglmm", ...)}, or a call to \code{cpglmm}.  Objects of class \code{"summary.cpglmm"} can be created by calls of the form \code{new("summary.cpglmm", ...)}, or a call to \code{summary}. 
}

\section{Slots}{
  This class extends \code{"\linkS4class{mer}"} directly, and it has the following addtional slots:
  \describe{
    \item{\code{p}:}{estimated value of the index parameter, class \code{"numeric"}  }
    \item{\code{phi}:}{estimated value of the dispersion parameter, class \code{"numeric"} }
    \item{\code{link.power}:}{the power of the link function, class \code{"numeric"}. See \code{\link[statmod]{tweedie}}.  }
    \item{\code{bound.p}:}{the specified bound of the index parameter, class \code{"numeric"}   }
    }

  The \code{"summary.cpglmm"} class \emph{contains} the \code{"cpglm"}
  class and has the following additional slots:
  \describe{

    \item{\code{methTitle}:}{character string specifying a method title}
    \item{\code{logLik}:}{the same as \code{logLik(object)}.}
    \item{\code{ngrps}:}{the number of levels per grouping factor in the
      \code{flist} slot.}
    \item{\code{sigma}:}{the scale factor for the variance-covariance estimates}
    \item{\code{coefs}:}{the matrix of estimates, standard errors,
      etc. for the fixed-effects coefficients}
    \item{\code{vcov}:}{the same as \code{\link{vcov}(object)}.}
    \item{\code{REmat}:}{the formatted Random-Effects matrix}
    \item{\code{AICtab}:}{A named vector of values of AIC, BIC, log-likelihood
      and deviance}
  }% summary slots
}

\section{Extends}{
The \code{"cpglmm"} code extends class \code{"\linkS4class{mer}"}, directly.

The \code{"summary.cpglmm"} code extends: class \code{"\linkS4class{cpglmm}"}, directly; 
class \code{"\linkS4class{mer}"}, by class "cpglmm", distance 2.
}


\section{Methods}{
  Most methods defined for \code{"\linkS4class{mer}"} can be applied directly, and the following is written specifically for \code{"cpglmm"}. 
  \describe{
    \item{print}{\code{signature(x = "cpglmm")}:  print the object }
    \item{show}{\code{signature(object = "cpglmm")}: show the object }
    \item{sigma}{\code{signature(object = "cpglmm")}: return square root of the dispersion parameter }
    \item{summary}{\code{signature(object = "cpglmm")}: summary results}
    \item{VarCorr}{\code{signature(x = "cpglmm")}: estimation for the variance components }
    \item{vcov}{\code{signature(object = "cpglmm")}: variance-covariance matrix for fixed effects }
	 }
}


\author{
Wayne (Yanwei) Zhang \email{actuary_zhang@hotmail.com}
}

\seealso{
\code{\link[lme4]{mer-class}} and \code{\link{cpglmm}}
}

\keyword{classes}
