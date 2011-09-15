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

\title{Class "bcplm" of MCMC compound Poisson linear models}
\description{This class corresponds to the result from a 
Markov Chain Monte Calo algorithm for the compound Poisson linear models
defined in this package. It extends class "cpglm" directly. 

}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("bcplm", ...)}.
}


\section{Slots}{
  \describe{
    \item{\code{n.chains}:}{Object of class \code{"integer"}, number of Markov chains (default: \code{3}).}
    \item{\code{n.iter}:}{Object of class \code{"integer"}, number of total iterations per chain (including burn in; default: \code{2000}).}
    \item{\code{n.burnin}:}{Object of class \code{"integer"}, length of burn in, i.e. number of iterations to discard at the beginning. Default is \code{n.iter/2}, that is, discarding the first half of the simulations.}
    \item{\code{n.thin}:}{Object of class \code{"integer"}, thinning rate. Must be a positive integer. Set \code{n.thin > 1} to save memory and computation time if \code{n.iter} is large. The default value is \code{max(1, floor(n.chains*(n.iter-n.burnin)/1000))} which will only thin if there are at least 2000 simulations.}
    \item{\code{n.sims}:}{Object of class \code{"integer"}, the approximate number of simulations to keep after thinning.}
    \item{\code{sims.list}:}{Object of class \code{"mcmc.list"}, a list of \code{n.chains} \code{mcmc} objects, each \code{mcmc} object storing the simulation result from a Markov chain. See \code{\link[coda]{mcmc}} and \code{\link[coda]{mcmc.convert}}. Since this is an \code{"mcmc.list"} object, most methods defined in the package \code{coda} can be directly applied to it. }
    \item{\code{summary}:}{Object of class \code{"summary.mcmc"}, two sets of summary statistics for each variable: Mean, standard deviation, naive standard error of the mean (ignoring autocorrelation of the chain) and time-series standard error based on an estimate of the spectral density at 0. See \code{\link[coda]{summary.mcmc}}.}
    \item{\code{link.power}:}{Object of class \code{"numeric"}, index of power link function. See \code{\link[statmod]{tweedie}}.}
    \item{\code{model.frame}:}{Object of class \code{"data.frame"}, the data frame used.}
    \item{\code{call}:}{Object of class \code{"call"}, the matched call. }
    \item{\code{formula}:}{Object of class \code{"formula"}, the formula supplied. }
    \item{\code{data}:}{Object of class \code{"data.frame"}, the supplied data.}
    \item{\code{contrasts}:}{Object of class \code{"NullList"}, the supplied contrasts. }
    \item{\code{inits}:}{Object of class \code{"list"}, initial values used for each chain. }
  }
}
\section{Methods}{
  \describe{
    \item{$}{\code{signature(x = "bcplm")}: see \code{\link{cpglm-class}}. }
    \item{[}{\code{signature(x = "bcplm", i = "character", j = "missing", drop = "missing")}: see \code{\link{cpglm-class}}. }
    \item{[}{\code{signature(x = "bcplm", i = "numeric", j = "missing", drop = "missing")}: see \code{\link{cpglm-class}}. }
    \item{[[}{\code{signature(x = "bcplm", i = "character", j = "missing")}: see \code{\link{cpglm-class}}. }
    \item{[[}{\code{signature(x = "bcplm", i = "numeric", j = "missing")}: see \code{\link{cpglm-class}}. }
    \item{formula}{\code{signature(x = "bcplm")}: see \code{\link{cpglm-class}}. }
    \item{model.matrix}{\code{signature(object = "bcplm")}: see \code{\link{cpglm-class}}. }
    \item{names}{\code{signature(x = "bcplm")}: see \code{\link{cpglm-class}}. }
    \item{plot}{\code{signature(x = "bcplm", y = "missing")}: see \code{\link[coda]{plot.mcmc}}. }
     \item{show}{\code{signature(object = "bcplm")}: Method for \code{show}. }
    \item{summary}{\code{signature(object = "bcplm")}: see \code{\link[coda]{summary.mcmc}}. }
    \item{terms}{\code{signature(x = "bcplm")}: see \code{\link{cpglm-class}}. }
	 }
}


\author{ Wayne  Zhang \email{actuary_zhang@hotmail.com} }
\seealso{
  See also \code{\link{bcpglm}}, \code{\link[coda]{mcmc}} and \code{\link[coda]{summary.mcmc}}.   
}

\keyword{classes}