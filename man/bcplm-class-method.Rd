\name{bcplm-class-method}
\docType{class}
\alias{bcplm-class}
\alias{bcpglm-class}
\alias{$,bcplm-method}
\alias{[,bcplm,character,missing,missing-method}
\alias{[,bcplm,numeric,missing,missing-method}
\alias{[[,bcplm,character,missing-method}
\alias{[[,bcplm,numeric,missing-method}
\alias{formula,bcplm-method}
\alias{model.matrix,bcplm-method}
\alias{names,bcplm-method}
\alias{plot,bcplm,missing-method}
\alias{show,bcplm-method}
\alias{summary,bcplm-method}
\alias{terms,bcplm-method}

\title{Representaion of a Bayesian compound Poisson Linear Model object}
\description{Class \code{"bcplm"} is a representation of a Bayesian compound Poisson Linear Model object resulted from implementations of   
Markov Chain Monte Calo algorithms.  Class \code{"bcpglm"} extends class \code{"bcplm"} directly. 

}
\section{Objects from the Class}{
Objects of class \code{"bcplm"} can be created by calls of the form \code{new("bcplm", ...)}. Objects of class \code{"bcpglm"} can be created by calls from \code{new("bcplm", ...)} or \code{bcpglm}. 
}


\section{Slots}{
  \describe{
    \item{\code{n.chains}:}{number of Markov chains, class \code{"integer"}}
    \item{\code{n.iter}:}{number of total iterations per chain, class \code{"integer"}}
    \item{\code{n.burnin}:}{length of burn in, i.e. number of iterations to discard at the beginning, class \code{"integer"}}
    \item{\code{n.thin}:}{thinning rate. class \code{"integer"}}
    \item{\code{n.sims}:}{the approximate number of simulations to keep after thinning, class \code{"integer"}}
    \item{\code{sims.list}:}{Object of class \code{"mcmc.list"}, a list of \code{n.chains} \code{mcmc} objects, each \code{mcmc} object storing the simulation result from a Markov chain. See \code{\link[coda]{mcmc}} and \code{\link[coda]{mcmc.convert}}. Since this is an \code{"mcmc.list"} object, most methods defined in the package \code{coda} can be directly applied to it. }
    \item{\code{summary}:}{Object of class \code{"summary.mcmc"}, two sets of summary statistics for each variable: Mean, standard deviation, naive standard error of the mean (ignoring autocorrelation of the chain) and time-series standard error based on an estimate of the spectral density at 0. See \code{\link[coda]{summary.mcmc}}.}
    \item{\code{link.power}:}{index of power link function, class \code{"numeric"}. See \code{\link[statmod]{tweedie}}.}
    \item{\code{model.frame}:}{the model frame used, class \code{"data.frame"}}
    \item{\code{call}:}{the matched call, class \code{"call"} }
    \item{\code{formula}:}{the formula supplied, class \code{"formula"} }
    \item{\code{data}:}{the supplied data, class \code{"data.frame"}}
    \item{\code{contrasts}:}{the supplied contrasts, class \code{"NullList"} }
    \item{\code{inits}:}{initial values used for each chain, class \code{"list"} }
  }
}
\section{Methods}{
  \describe{
    \item{$}{\code{signature(x = "bcplm")}: extract a slot of \code{x} with a specified slot name, see \code{\link{cpglm-class}}. }
    \item{[}{\code{signature(x = "bcplm", i = "character", j = "missing", drop = "missing")}:  see \code{\link{cpglm-class}}. }
    \item{[}{\code{signature(x = "bcplm", i = "numeric", j = "missing", drop = "missing")}:  see \code{\link{cpglm-class}}. }
    \item{[[}{\code{signature(x = "bcplm", i = "character", j = "missing")}: extract the slots with names in \code{i}, see \code{\link{cpglm-class}}. }
    \item{[[}{\code{signature(x = "bcplm", i = "numeric", j = "missing")}: extract the i-th slot, see \code{\link{cpglm-class}}. }
    \item{formula}{\code{signature(x = "bcplm")}: extract the formula }
    \item{model.matrix}{\code{signature(object = "bcplm")}: extract the model matrix }
    \item{names}{\code{signature(x = "bcplm")}: extract the slot names }
    \item{plot}{\code{signature(x = "bcplm", y = "missing")}: summarize the \code{"bcplm"} object with a trace of the sampled output and a density estimate for each variable in the chain. See \code{\link[coda]{plot.mcmc}}. }
     \item{show}{\code{signature(object = "bcplm")}: method for \code{show}. }
    \item{summary}{\code{signature(object = "bcplm")}:  produce two sets of summary statistics. See \code{\link[coda]{summary.mcmc}}. }
    \item{terms}{\code{signature(x = "bcplm")}: extract the model terms. }
	 }
}


\author{ Wayne  Zhang \email{actuary_zhang@hotmail.com} }
\seealso{
  See also \code{\link{bcpglm}}, \code{\link[coda]{mcmc}} and \code{\link[coda]{summary.mcmc}}.   
}

\keyword{classes}